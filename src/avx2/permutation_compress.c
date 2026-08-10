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

#include "permutation_compress.h"

#include <immintrin.h>
#include <string.h>

#include "permutation_compress_tables.h"

static unsigned bit_length_u64(uint64_t x) {
    unsigned bits = 0;
    do {
        bits++;
        x >>= 1;
    } while (x != 0);
    return bits;
}

static uint64_t block_bound(unsigned block_idx) {
    uint64_t bound = 1;
    for (unsigned i = celeres_perm_j[block_idx - 1u] + 1u;
         i <= celeres_perm_j[block_idx]; i++) {
        bound *= i;
    }
    return bound;
}

static unsigned block_bits(unsigned block_idx) {
    return bit_length_u64(block_bound(block_idx) - 1u);
}

static unsigned store_bits(uint8_t out[CELERES_PERM_COMPRESSED_BYTES],
                           unsigned bit_index,
                           uint64_t value,
                           unsigned bits) {
    for (unsigned i = 0; i < bits; i++) {
        out[(bit_index + i) / 8u] |=
            (uint8_t)(((value >> i) & 1u) << ((bit_index + i) % 8u));
    }
    return bit_index + bits;
}

static unsigned load_bits(uint64_t *value,
                          const uint8_t in[CELERES_PERM_COMPRESSED_BYTES],
                          unsigned bit_index,
                          unsigned bits) {
    *value = 0;
    for (unsigned i = 0; i < bits; i++) {
        *value |= (uint64_t)((in[(bit_index + i) / 8u] >> ((bit_index + i) % 8u)) & 1u) << i;
    }
    return bit_index + bits;
}

static uint32_t lehmer_digit_avx2(const POSITION_T permutation[N], unsigned digit_idx) {
    const unsigned pos = N - 1u - digit_idx;
    const __m256i value = _mm256_set1_epi16((int16_t)permutation[pos]);
    uint32_t count = 0;
    unsigned j = pos + 1u;

    for (; j + 16u <= N; j += 16u) {
        const __m256i chunk = _mm256_loadu_si256((const __m256i *)(const void *)(permutation + j));
        const __m256i cmp = _mm256_cmpgt_epi16(value, chunk);
        count += (uint32_t)__builtin_popcount((unsigned)_mm256_movemask_epi8(cmp)) / 2u;
    }
    for (; j < N; j++) {
        if (permutation[j] < permutation[pos]) {
            count++;
        }
    }
    return count;
}

static unsigned encode_block(uint8_t out[CELERES_PERM_COMPRESSED_BYTES],
                             unsigned bit_index,
                             unsigned block_idx,
                             const POSITION_T permutation[N]) {
    uint64_t block = 0;
    const unsigned start = celeres_perm_j[block_idx - 1u];
    const unsigned end = celeres_perm_j[block_idx];

    for (unsigned i = end - 1u; i > start - 1u; i--) {
        block = lehmer_digit_avx2(permutation, i) + (uint64_t)(i + 1u) * block;
    }
    return store_bits(out, bit_index, block, block_bits(block_idx));
}

static int decode_block(POSITION_T digits[N],
                        const uint8_t in[CELERES_PERM_COMPRESSED_BYTES],
                        unsigned *bit_index,
                        unsigned block_idx) {
    uint64_t block = 0;
    const uint64_t bound = block_bound(block_idx);
    const unsigned start = celeres_perm_j[block_idx - 1u];
    const unsigned end = celeres_perm_j[block_idx];
    unsigned i;

    *bit_index = load_bits(&block, in, *bit_index, block_bits(block_idx));
    if (block >= bound) {
        return -1;
    }

    for (i = start + 1u; i < end; i++) {
        const uint64_t quotient = block / i;
        digits[N - i] = (POSITION_T)(block - quotient * i);
        block = quotient;
    }
    digits[N - i] = (POSITION_T)block;
    return 0;
}

static int padding_is_zero(const uint8_t in[CELERES_PERM_COMPRESSED_BYTES]) {
    for (unsigned bit = CELERES_PERM_COMPRESSED_BITS;
         bit < CELERES_PERM_COMPRESSED_BYTES * 8u; bit++) {
        if (((in[bit / 8u] >> (bit % 8u)) & 1u) != 0) {
            return 0;
        }
    }
    return 1;
}

void celeres_perm_compress(uint8_t out[CELERES_PERM_COMPRESSED_BYTES],
                            const POSITION_T permutation[N]) {
    unsigned bit_index = 0;

    memset(out, 0, CELERES_PERM_COMPRESSED_BYTES);
    for (unsigned block_idx = 1; block_idx <= CELERES_PERM_BLOCK_COUNT; block_idx++) {
        bit_index = encode_block(out, bit_index, block_idx, permutation);
    }
}

int celeres_perm_decompress(POSITION_T permutation[N],
                             const uint8_t in[CELERES_PERM_COMPRESSED_BYTES]) {
    uint8_t used[N] = {0};
    unsigned bit_index = 0;

    if (!padding_is_zero(in)) {
        return -1;
    }

    permutation[N - 1u] = 0;
    for (unsigned block_idx = 1; block_idx <= CELERES_PERM_BLOCK_COUNT; block_idx++) {
        if (decode_block(permutation, in, &bit_index, block_idx) != 0) {
            return -1;
        }
    }

    for (unsigned i = 0; i < N; i++) {
        POSITION_T count = 0;
        int found = 0;

        for (POSITION_T j = 0; j < N; j++) {
            if (used[j] == 0) {
                count++;
            }
            if (count == permutation[i] + 1u) {
                permutation[i] = j;
                used[j] = 1;
                found = 1;
                break;
            }
        }
        if (!found) {
            return -1;
        }
    }
    return 0;
}
