/**
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **/

#include "celeres_util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rng.h"
#include "sha3.h"

#define DOMAIN_BITS "CELERES-BITS-v1"
#define DOMAIN_DUMMY "CELERES-DUMMY-v1"

int celeres_randombytes(uint8_t *out, size_t len) {
    FILE *f = fopen("/dev/urandom", "rb");
    if (f == NULL) {
        return -1;
    }
    const size_t got = fread(out, 1, len, f);
    fclose(f);
    return got == len ? 0 : -1;
}

int celeres_log_round_up_u32(uint32_t n) {
    int logn = 0;
    uint32_t v = 1;
    while (v < n) {
        v <<= 1;
        logn++;
    }
    return logn;
}

static void xof_bytes(uint8_t *out, size_t out_len,
                      const uint8_t *in, size_t in_len) {
    SHAKE_STATE_STRUCT st;
    xof_shake_init(&st, SEED_LENGTH_BYTES * 8);
    xof_shake_update(&st, in, in_len);
    xof_shake_final(&st);
    xof_shake_extract(&st, out, (unsigned int)out_len);
}

void celeres_derive_bits(uint8_t out[SEED_LENGTH_BYTES],
                         const uint8_t seed[SEED_LENGTH_BYTES],
                         const uint8_t salt[HASH_DIGEST_LENGTH],
                         uint16_t round,
                         uint32_t ring_index) {
    SHAKE_STATE_STRUCT st;
    xof_shake_init(&st, SEED_LENGTH_BYTES * 8);
    xof_shake_update(&st, (const uint8_t *)DOMAIN_BITS, sizeof(DOMAIN_BITS) - 1);
    xof_shake_update(&st, seed, SEED_LENGTH_BYTES);
    xof_shake_update(&st, salt, HASH_DIGEST_LENGTH);
    xof_shake_update(&st, (const uint8_t *)&round, sizeof(round));
    xof_shake_update(&st, (const uint8_t *)&ring_index, sizeof(ring_index));
    xof_shake_final(&st);
    xof_shake_extract(&st, out, SEED_LENGTH_BYTES);
}

void celeres_dummy_commitment(uint8_t out[HASH_DIGEST_LENGTH],
                              const uint8_t salt[HASH_DIGEST_LENGTH],
                              uint16_t round,
                              uint32_t dummy_index) {
    SHAKE_STATE_STRUCT st;
    xof_shake_init(&st, SEED_LENGTH_BYTES * 8);
    xof_shake_update(&st, (const uint8_t *)DOMAIN_DUMMY, sizeof(DOMAIN_DUMMY) - 1);
    xof_shake_update(&st, salt, HASH_DIGEST_LENGTH);
    xof_shake_update(&st, (const uint8_t *)&round, sizeof(round));
    xof_shake_update(&st, (const uint8_t *)&dummy_index, sizeof(dummy_index));
    xof_shake_final(&st);
    xof_shake_extract(&st, out, HASH_DIGEST_LENGTH);
}

static uint8_t ct_eq_u32_mask(uint32_t a, uint32_t b) {
    uint32_t x = a ^ b;
    x |= (uint32_t)(0u - x);
    return (uint8_t)((((x >> 31) ^ 1u) & 1u) * 0xffu);
}

static uint32_t ct_lt_u8(uint8_t a, uint8_t b) {
    return (((uint32_t)a - (uint32_t)b) >> 31) & 1u;
}

static uint32_t ct_hash_gt(const uint8_t a[HASH_DIGEST_LENGTH],
                           const uint8_t b[HASH_DIGEST_LENGTH]) {
    uint32_t gt = 0;
    uint32_t lt = 0;
    for (size_t i = 0; i < HASH_DIGEST_LENGTH; i++) {
        const uint32_t undecided = (gt | lt) ^ 1u;
        gt |= undecided & ct_lt_u8(b[i], a[i]);
        lt |= undecided & ct_lt_u8(a[i], b[i]);
    }
    return gt;
}

static void ct_cond_swap_hash(uint8_t a[HASH_DIGEST_LENGTH],
                              uint8_t b[HASH_DIGEST_LENGTH],
                              uint32_t swap) {
    const uint8_t mask = (uint8_t)(0u - (swap & 1u));
    for (size_t i = 0; i < HASH_DIGEST_LENGTH; i++) {
        const uint8_t diff = (uint8_t)((a[i] ^ b[i]) & mask);
        a[i] ^= diff;
        b[i] ^= diff;
    }
}

static void ct_select_hash(uint8_t out[HASH_DIGEST_LENGTH],
                           const uint8_t *level,
                           uint32_t width,
                           uint32_t index) {
    memset(out, 0, HASH_DIGEST_LENGTH);
    for (uint32_t i = 0; i < width; i++) {
        const uint8_t mask = ct_eq_u32_mask(i, index);
        const uint8_t *candidate = level + (size_t)i * HASH_DIGEST_LENGTH;
        for (size_t j = 0; j < HASH_DIGEST_LENGTH; j++) {
            out[j] |= (uint8_t)(candidate[j] & mask);
        }
    }
}

void celeres_ihmt_build_root_and_path(const uint8_t *commitments_in,
                                      int logn,
                                      int32_t index,
                                      uint8_t root[HASH_DIGEST_LENGTH],
                                      uint8_t *path) {
    const uint32_t leaves = (uint32_t)1 << logn;
    uint8_t *level = malloc((size_t)leaves * HASH_DIGEST_LENGTH);
    uint8_t pair[2 * HASH_DIGEST_LENGTH];
    memcpy(level, commitments_in, (size_t)leaves * HASH_DIGEST_LENGTH);

    if (logn == 0) {
        memcpy(root, level, HASH_DIGEST_LENGTH);
        free(level);
        return;
    }

    for (int depth = logn - 1; depth >= 0; depth--) {
        const uint32_t width = (uint32_t)1 << (depth + 1);
        if (index >= 0 && path != NULL) {
            ct_select_hash(path + (size_t)depth * HASH_DIGEST_LENGTH,
                           level, width, (uint32_t)index ^ 1u);
            index >>= 1;
        }
        for (uint32_t i = 0; i < width / 2; i++) {
            uint8_t *left = level + (2u * i) * HASH_DIGEST_LENGTH;
            uint8_t *right = level + (2u * i + 1u) * HASH_DIGEST_LENGTH;
            ct_cond_swap_hash(left, right, ct_hash_gt(left, right));
            memcpy(pair, left, HASH_DIGEST_LENGTH);
            memcpy(pair + HASH_DIGEST_LENGTH, right, HASH_DIGEST_LENGTH);
            xof_bytes(level + i * HASH_DIGEST_LENGTH, HASH_DIGEST_LENGTH,
                      pair, sizeof(pair));
        }
    }
    memcpy(root, level, HASH_DIGEST_LENGTH);
    free(level);
}

void celeres_ihmt_reconstruct_root(const uint8_t leaf[HASH_DIGEST_LENGTH],
                                   const uint8_t *path,
                                   int logn,
                                   uint8_t root[HASH_DIGEST_LENGTH]) {
    uint8_t pair[2 * HASH_DIGEST_LENGTH];
    memcpy(pair, leaf, HASH_DIGEST_LENGTH);
    for (int depth = logn - 1; depth >= 0; depth--) {
        memcpy(pair + HASH_DIGEST_LENGTH,
               path + (size_t)depth * HASH_DIGEST_LENGTH,
               HASH_DIGEST_LENGTH);
        ct_cond_swap_hash(pair, pair + HASH_DIGEST_LENGTH,
                          ct_hash_gt(pair, pair + HASH_DIGEST_LENGTH));
        xof_bytes(pair, HASH_DIGEST_LENGTH, pair, sizeof(pair));
    }
    memcpy(root, pair, HASH_DIGEST_LENGTH);
}

void celeres_sample_challenge(uint8_t challenge[T],
                              const uint8_t digest[HASH_DIGEST_LENGTH]) {
    memset(challenge, 0, T);
    SHAKE_STATE_STRUCT st;
    initialize_csprng(&st, digest, HASH_DIGEST_LENGTH);

    uint32_t filled = 0;
    while (filled < W) {
        uint16_t candidate;
        csprng_randombytes((uint8_t *)&candidate, sizeof(candidate), &st);
        candidate = (uint16_t)(candidate % T);
        if (challenge[candidate] == 0) {
            challenge[candidate] = 1;
            filled++;
        }
    }
}
