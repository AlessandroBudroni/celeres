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

#include "celeres.h"

#include <stdlib.h>
#include <string.h>

#include "codes.h"
#include "fq_arith.h"
#include "celeres_internal.h"
#include "celeres_util.h"
#include "monomial_mat.h"
#include "permutation_compress.h"
#include "rng.h"
#include "seedtree.h"
#include "sha3.h"
#include "utils.h"

#define DOMAIN_CELERES_COMMIT "CELERES-COMMIT-v1"
#define DOMAIN_CELERES_TAU "CELERES-TAU-v1"
#define DOMAIN_CELERES_FS "CELERES-FS-v1"

static size_t celeres_base_bytes(void) {
    return CELERES_TAG_BYTES + (2u * HASH_DIGEST_LENGTH) + 1u;
}

static size_t celeres_response_bytes(void) {
    return CELERES_PERM_BYTES + SEED_LENGTH_BYTES +
           ((size_t)celeres_log_round_up_u32(CELERES_RING_SIZE) * HASH_DIGEST_LENGTH);
}

static size_t celeres_compact_bytes(uint32_t seed_count) {
    if (seed_count > MAX_PUBLISHED_SEEDS) {
        return 0;
    }
    return celeres_base_bytes() +
           ((size_t)seed_count * SEED_LENGTH_BYTES) +
           ((size_t)W * celeres_response_bytes());
}

size_t celeres_signature_bytes(void) {
    return celeres_compact_bytes(MAX_PUBLISHED_SEEDS);
}

static uint8_t *sig_tag(uint8_t *sig) {
    return sig;
}

static uint8_t *sig_digest(uint8_t *sig) {
    return sig + CELERES_TAG_BYTES;
}

static uint8_t *sig_salt(uint8_t *sig) {
    return sig + CELERES_TAG_BYTES + HASH_DIGEST_LENGTH;
}

static uint8_t *sig_seed_count(uint8_t *sig) {
    return sig + CELERES_TAG_BYTES + 2u * HASH_DIGEST_LENGTH;
}

static uint8_t *sig_seeds(uint8_t *sig) {
    return sig + celeres_base_bytes();
}

static uint8_t *sig_responses(uint8_t *sig, uint32_t seed_count) {
    return sig_seeds(sig) + (size_t)seed_count * SEED_LENGTH_BYTES;
}

static const uint8_t *csig_tag(const uint8_t *sig) {
    return sig;
}

static const uint8_t *csig_digest(const uint8_t *sig) {
    return sig + CELERES_TAG_BYTES;
}

static const uint8_t *csig_salt(const uint8_t *sig) {
    return sig + CELERES_TAG_BYTES + HASH_DIGEST_LENGTH;
}

static const uint8_t *csig_seed_count(const uint8_t *sig) {
    return sig + CELERES_TAG_BYTES + 2u * HASH_DIGEST_LENGTH;
}

static const uint8_t *csig_seeds(const uint8_t *sig) {
    return sig + celeres_base_bytes();
}

static const uint8_t *csig_responses(const uint8_t *sig, uint32_t seed_count) {
    return csig_seeds(sig) + (size_t)seed_count * SEED_LENGTH_BYTES;
}

static int tag_to_canonical_form(generator_mat_t *cf,
                                 const uint8_t tag[CELERES_TAG_BYTES]) {
    uint8_t tag_rref[RREF_MAT_PACKEDBYTES];
    uint8_t pivots[N_pad] = {0};

#if CELERES_TAG_SELF_ORTHOGONAL
    if (celeres_tag_expand_to_rref(tag_rref, tag) != 0) {
        return -1;
    }
#else
    memcpy(tag_rref, tag, RREF_MAT_PACKEDBYTES);
#endif

    memset(cf, 0, sizeof(*cf));
    expand_to_rref(cf, tag_rref, pivots);
    celeres_diagonal_scaling_canonical_form(cf);
    return 0;
}

static int canonical_forms_equal(const generator_mat_t *a,
                                 const generator_mat_t *b) {
    unsigned char diff = 0;
    for (uint32_t row = 0; row < K; row++) {
        diff |= (unsigned char)verify(a->values[row], b->values[row], N);
    }
    return diff == 0;
}

static int challenge_is_valid(const uint8_t challenge[T]) {
    uint32_t weight = 0;
    for (uint32_t i = 0; i < T; i++) {
        if (challenge[i] > 1) {
            return 0;
        }
        weight += challenge[i];
    }
    return weight == W;
}

static uint32_t expected_published_seed_count(const uint8_t challenge[T]) {
    uint8_t dummy_tree[NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES] = {0};
    uint8_t dummy_storage[SEED_TREE_MAX_PUBLISHED_BYTES] = {0};
    return GGMPath(dummy_tree, challenge, dummy_storage);
}

static void compose_monomial(monomial_t *out,
                             const monomial_t *left,
                             const monomial_t *right) {
    for (uint32_t i = 0; i < N; i++) {
        const POSITION_T mid = left->permutation[i];
        out->permutation[i] = right->permutation[mid];
        out->coefficients[i] = fq_mul(left->coefficients[i], right->coefficients[mid]);
    }
}

#if !CELERES_TAG_SELF_ORTHOGONAL
static void sample_q_tilde(monomial_t *out,
                           const monomial_t *q,
                           const uint8_t coeff_seed[CELERES_SECRET_SEED_BYTES]) {
    SHAKE_STATE_STRUCT st;
    memcpy(out->permutation, q->permutation, sizeof(out->permutation));
    initialize_csprng(&st, coeff_seed, CELERES_SECRET_SEED_BYTES);
    fq_star_rnd_state_elements(&st, out->coefficients, N);
}
#endif

static void rref_compact_pivot_flags(uint8_t pivots[N_pad],
                                     const rref_generator_mat_t *g) {
    memset(pivots, 0, N_pad);
    for (uint32_t col = 0; col < N; col++) {
        pivots[col] = 1;
    }
    for (uint32_t i = 0; i < N - K; i++) {
        pivots[g->column_pos[i]] = 0;
    }
}

static void permutation_to_bytes(uint8_t out[CELERES_PERM_BYTES],
                                 const POSITION_T permutation[N]) {
    celeres_perm_compress(out, permutation);
}

static int bytes_to_permutation(POSITION_T permutation[N],
                                const uint8_t in[CELERES_PERM_BYTES]) {
    return celeres_perm_decompress(permutation, in);
}

static void monomial_from_permutation(monomial_t *out,
                                      const POSITION_T permutation[N]) {
    for (uint32_t i = 0; i < N; i++) {
        out->permutation[i] = permutation[i];
        out->coefficients[i] = 1;
    }
}

static void public_seed_rref_bytes(uint8_t out[RREF_MAT_PACKEDBYTES],
                                   const uint8_t seed[SEED_LENGTH_BYTES]) {
    rref_generator_mat_t compact;
    generator_mat_t full = {0};
    uint8_t pivots[N_pad] = {0};

    generator_sample(&compact, seed);
    generator_rref_expand(&full, &compact);
    rref_compact_pivot_flags(pivots, &compact);
    compress_rref(out, &full, pivots);
}

static int sf_and_compress_from_rref_monomial(
    uint8_t out[RREF_MAT_PACKEDBYTES],
    const uint8_t rref[RREF_MAT_PACKEDBYTES],
    const monomial_t *q) {
    generator_mat_t systematic = {0};
    uint8_t pivots[N_pad] = {0};

    if (celeres_sf_from_rref_monomial(&systematic, pivots, rref, q) != 0) {
        return -1;
    }
    compress_rref(out, &systematic, pivots);
    return 0;
}

static int sf_and_compress_from_public_seed_monomial(
    uint8_t out[RREF_MAT_PACKEDBYTES],
    const uint8_t seed[SEED_LENGTH_BYTES],
    const monomial_t *q) {
    uint8_t base_rref[RREF_MAT_PACKEDBYTES];

    public_seed_rref_bytes(base_rref, seed);
    return sf_and_compress_from_rref_monomial(out, base_rref, q);
}

#if CELERES_TAG_SELF_ORTHOGONAL
static int sf_and_compress_tag_from_seed_permutation(
    uint8_t out[CELERES_TAG_BYTES],
    const uint8_t seed[SEED_LENGTH_BYTES],
    const POSITION_T permutation[N]) {
    uint8_t base_rref[RREF_MAT_PACKEDBYTES];
    uint8_t tag_rref[RREF_MAT_PACKEDBYTES];
    monomial_t p;

    celeres_tag_sample_base_rref(base_rref, seed);
    monomial_from_permutation(&p, permutation);
    if (sf_and_compress_from_rref_monomial(tag_rref, base_rref, &p) != 0) {
        return -1;
    }
    return celeres_tag_compress_from_rref(out, tag_rref);
}
#endif

static void absorb_cf_matrix(LESS_SHA3_INC_CTX *st,
                             const generator_mat_t *cf) {
    for (uint32_t row = 0; row < K; row++) {
        LESS_SHA3_INC_ABSORB(st, cf->values[row], N);
    }
}

static void hash_commitment(uint8_t out[HASH_DIGEST_LENGTH],
                            const uint8_t bits[SEED_LENGTH_BYTES],
                            const generator_mat_t *cf) {
    LESS_SHA3_INC_CTX st;
    LESS_SHA3_INC_INIT(&st);
    LESS_SHA3_INC_ABSORB(&st, (const uint8_t *)DOMAIN_CELERES_COMMIT,
                         sizeof(DOMAIN_CELERES_COMMIT) - 1);
    LESS_SHA3_INC_ABSORB(&st, bits, SEED_LENGTH_BYTES);
    absorb_cf_matrix(&st, cf);
    LESS_SHA3_INC_FINALIZE(out, &st);
}

static void hash_tau(uint8_t out[HASH_DIGEST_LENGTH],
                     const uint8_t root[HASH_DIGEST_LENGTH],
                     const generator_mat_t *tau) {
    LESS_SHA3_INC_CTX st;
    LESS_SHA3_INC_INIT(&st);
    LESS_SHA3_INC_ABSORB(&st, (const uint8_t *)DOMAIN_CELERES_TAU,
                         sizeof(DOMAIN_CELERES_TAU) - 1);
    LESS_SHA3_INC_ABSORB(&st, root, HASH_DIGEST_LENGTH);
    absorb_cf_matrix(&st, tau);
    LESS_SHA3_INC_FINALIZE(out, &st);
}

static int compute_cf_from_rref_and_monomial(generator_mat_t *cf,
                                             const uint8_t rref[RREF_MAT_PACKEDBYTES],
                                             const monomial_t *q) {
    const int use_ct_scaling = 1;
    return celeres_cf_from_rref_monomial_is(cf, rref, q, use_ct_scaling);
}

static int compute_cf_from_rref_and_monomial_with_reuse(
    generator_mat_t *cf,
    const uint8_t rref[RREF_MAT_PACKEDBYTES],
    const monomial_t *q) {
    return celeres_cf_from_rref_monomial_is(cf, rref, q, 0);
}

static int compute_cf_from_rref_and_perm_with_reuse(
    generator_mat_t *cf,
    const uint8_t base_rref[RREF_MAT_PACKEDBYTES],
    const POSITION_T permutation[N]) {
    monomial_t p;

    monomial_from_permutation(&p, permutation);
    return celeres_cf_from_rref_monomial_is(cf, base_rref, &p, 0);
}

static int compute_commitment_from_pk(uint8_t out[HASH_DIGEST_LENGTH],
                                      const celeres_public_key_t *pk,
                                      const monomial_t *qbar,
                                      const uint8_t bits[SEED_LENGTH_BYTES]) {
    generator_mat_t cf = {0};
    if (compute_cf_from_rref_and_monomial(&cf, pk->rref, qbar) != 0) {
        return -1;
    }
    hash_commitment(out, bits, &cf);
    return 0;
}

static int compute_commitment_from_pk_with_reuse(
    uint8_t out[HASH_DIGEST_LENGTH],
    const celeres_public_key_t *pk,
    const monomial_t *qbar,
    const uint8_t bits[SEED_LENGTH_BYTES]) {
    generator_mat_t cf = {0};
    if (compute_cf_from_rref_and_monomial_with_reuse(&cf, pk->rref, qbar) != 0) {
        return -1;
    }
    hash_commitment(out, bits, &cf);
    return 0;
}

static int compute_commitment_from_perm(uint8_t out[HASH_DIGEST_LENGTH],
                                        const uint8_t base_g_rref[RREF_MAT_PACKEDBYTES],
                                        const POSITION_T permutation[N],
                                        const uint8_t bits[SEED_LENGTH_BYTES]) {
    generator_mat_t cf = {0};
    if (compute_cf_from_rref_and_perm_with_reuse(&cf, base_g_rref,
                                                 permutation) != 0) {
        return -1;
    }
    hash_commitment(out, bits, &cf);
    return 0;
}

static int compute_tau_from_perm(generator_mat_t *tau,
                                 const uint8_t base_t_rref[RREF_MAT_PACKEDBYTES],
                                 const POSITION_T permutation[N]) {
    return compute_cf_from_rref_and_perm_with_reuse(tau, base_t_rref,
                                                    permutation);
}

static void fs_digest(uint8_t digest[HASH_DIGEST_LENGTH],
                      const celeres_public_parameters_t *params,
                      const celeres_public_key_t *ring,
                      const uint8_t tag[CELERES_TAG_BYTES],
                      const uint8_t *message,
                      size_t message_len,
                      const uint8_t salt[HASH_DIGEST_LENGTH],
                      const uint8_t *round_hashes) {
    LESS_SHA3_INC_CTX st;
    LESS_SHA3_INC_INIT(&st);
    LESS_SHA3_INC_ABSORB(&st, (const uint8_t *)DOMAIN_CELERES_FS, sizeof(DOMAIN_CELERES_FS) - 1);
    LESS_SHA3_INC_ABSORB(&st, salt, HASH_DIGEST_LENGTH);
    LESS_SHA3_INC_ABSORB(&st, params->g_seed, SEED_LENGTH_BYTES);
    LESS_SHA3_INC_ABSORB(&st, params->t_seed, SEED_LENGTH_BYTES);
    const uint32_t ring_size = CELERES_RING_SIZE;
    LESS_SHA3_INC_ABSORB(&st, (const uint8_t *)&ring_size, sizeof(ring_size));
    LESS_SHA3_INC_ABSORB(&st, (const uint8_t *)ring,
                         (size_t)CELERES_RING_SIZE * sizeof(*ring));
    LESS_SHA3_INC_ABSORB(&st, tag, CELERES_TAG_BYTES);
    LESS_SHA3_INC_ABSORB(&st, message, message_len);
    LESS_SHA3_INC_ABSORB(&st, round_hashes, (size_t)T * HASH_DIGEST_LENGTH);
    LESS_SHA3_INC_FINALIZE(digest, &st);
}

int celeres_generate_parameters(celeres_public_parameters_t *params) {
    if (params == NULL) {
        return -1;
    }
    return celeres_randombytes(params->g_seed, sizeof(params->g_seed)) != 0 ||
           celeres_randombytes(params->t_seed, sizeof(params->t_seed)) != 0
               ? -1
               : 0;
}

int celeres_keygen(const celeres_public_parameters_t *params,
                        celeres_public_key_t *pk,
                        celeres_tag_t *tag,
                        celeres_secret_key_t *sk) {
    if (params == NULL || pk == NULL || tag == NULL || sk == NULL) {
        return -1;
    }
    /* KEYGEN.1: G and T are public matrices, expanded from public seeds. */

    /* KEYGEN.2.i-ii: sample Q and Qtilde in the same PermRep class. */
    if (celeres_randombytes(sk->q_seed, sizeof(sk->q_seed)) != 0 ||
        celeres_randombytes(sk->q_tilde_seed, sizeof(sk->q_tilde_seed)) != 0) {
        return -1;
    }

    monomial_t q;

    monomial_sample_prikey(&q, sk->q_seed);
#if !CELERES_TAG_SELF_ORTHOGONAL
    monomial_t q_tilde;
    sample_q_tilde(&q_tilde, &q, sk->q_tilde_seed);
#endif

    /* KEYGEN.2.iii-iv: pk = SF_{Q([k])}(G Q). */
    if (sf_and_compress_from_public_seed_monomial(pk->rref,
                                                  params->g_seed, &q) != 0) {
        return -1;
    }

#if CELERES_TAG_SELF_ORTHOGONAL
    /* KEYGEN.2.v: tag = SF_{Q([k])}(T P_Q), using only Q's permutation. */
    return sf_and_compress_tag_from_seed_permutation(tag->tag,
                                                     params->t_seed,
                                                     q.permutation);
#else
    /* KEYGEN.2.v: old full-RREF tag path, tag = SF(T Qtilde). */
    return sf_and_compress_from_public_seed_monomial(tag->tag,
                                                     params->t_seed,
                                                     &q_tilde);
#endif
}

int celeres_sign(const celeres_public_parameters_t *params,
                      const celeres_public_key_t *ring,
                      const celeres_tag_t *tags,
                      uint32_t signer_index,
                      const celeres_secret_key_t *sk,
                      const uint8_t *message,
                      size_t message_len,
                      uint8_t *signature,
                      size_t *signature_len) {
    if (params == NULL || ring == NULL || tags == NULL || sk == NULL ||
        signature == NULL || signature_len == NULL ||
        signer_index >= CELERES_RING_SIZE) {
        return -1;
    }

    const int logn = celeres_log_round_up_u32(CELERES_RING_SIZE);
    const uint32_t rounded_ring = (uint32_t)1 << logn;
    const size_t path_bytes = (size_t)logn * HASH_DIGEST_LENGTH;
    const size_t sig_cap = celeres_signature_bytes();
    memset(signature, 0, sig_cap);

    uint8_t *tag_out = sig_tag(signature);
    uint8_t *stored_digest = sig_digest(signature);
    uint8_t *salt = sig_salt(signature);
    uint8_t *seed_count_ptr = sig_seed_count(signature);
    uint8_t *seed_storage = sig_seeds(signature);
    uint8_t challenge[T];

    /* SIGN output prefix: include T_{i*} in the signature. */
    memcpy(tag_out, tags[signer_index].tag, CELERES_TAG_BYTES);
    uint8_t signer_tag_rref[RREF_MAT_PACKEDBYTES];
#if CELERES_TAG_SELF_ORTHOGONAL
    if (celeres_tag_expand_to_rref(signer_tag_rref,
                                   tags[signer_index].tag) != 0) {
        return -1;
    }
#else
    memcpy(signer_tag_rref, tags[signer_index].tag, RREF_MAT_PACKEDBYTES);
#endif

    /* SIGN.2-3: sample root seed and Fiat-Shamir salt. */
    uint8_t root_seed[SEED_LENGTH_BYTES];
    if (celeres_randombytes(root_seed, sizeof(root_seed)) != 0 ||
        celeres_randombytes(salt, HASH_DIGEST_LENGTH) != 0) {
        return -1;
    }

    /* SIGN.4: expand the seed tree to one seed per round j. */
    uint8_t seed_tree[NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES] = {0};
    uint8_t round_seeds[T * SEED_LENGTH_BYTES] = {0};
    BuildGGM(seed_tree, root_seed, salt);
    seed_leaves(round_seeds, seed_tree);

    uint8_t *round_hashes = malloc((size_t)T * HASH_DIGEST_LENGTH);
    uint8_t *commitments = malloc((size_t)rounded_ring * HASH_DIGEST_LENGTH);
    uint8_t *paths = calloc(T == 0 ? 1 : T, path_bytes == 0 ? 1 : path_bytes);
    uint8_t (*signer_bits)[SEED_LENGTH_BYTES] = calloc(T, SEED_LENGTH_BYTES);
    POSITION_T (*signer_perms)[N] = calloc(T, sizeof(*signer_perms));
    if (round_hashes == NULL || commitments == NULL || paths == NULL ||
        signer_bits == NULL || signer_perms == NULL) {
        free(round_hashes);
        free(commitments);
        free(paths);
        free(signer_bits);
        free(signer_perms);
        return -1;
    }

    monomial_t private_q;
    monomial_t qbar;
    monomial_t product;
    /* SIGN.1: parse sk_{i*}; the private monomial Q is seed-derived. */
    monomial_sample_prikey(&private_q, sk->q_seed);

    /* SIGN.5: build each round transcript h_j. */
    for (uint32_t j = 0; j < T; j++) {
        const uint8_t *seed_j = round_seeds + (size_t)j * SEED_LENGTH_BYTES;

        /* SIGN.5.i: ExpandRand gives Qbar_j and r_{j,i} values. */
        monomial_sample_salt(&qbar, seed_j, salt, (uint16_t)j);

        for (uint32_t i = 0; i < CELERES_RING_SIZE; i++) {
            uint8_t bits[SEED_LENGTH_BYTES];
            celeres_derive_bits(bits, seed_j, salt, (uint16_t)j, i);

            /* SIGN.5.ii.a: h_{j,i}=com(r_{j,i}, CF(G_i Qbar_j)). */
            if (compute_commitment_from_pk(commitments + (size_t)i * HASH_DIGEST_LENGTH,
                                           &ring[i], &qbar, bits) != 0) {
                free(round_hashes);
                free(commitments);
                free(paths);
                free(signer_bits);
                free(signer_perms);
                return -1;
            }
            if (i == signer_index) {
                memcpy(signer_bits[j], bits, SEED_LENGTH_BYTES);
            }
        }

        /* Pad non-power-of-two rings with deterministic dummy leaves. */
        for (uint32_t i = CELERES_RING_SIZE; i < rounded_ring; i++) {
            celeres_dummy_commitment(commitments + (size_t)i * HASH_DIGEST_LENGTH,
                                      salt, (uint16_t)j, i);
        }

        /* SIGN.5.iii and SIGN.7.ii: build rt_j and save Path_j for signer. */
        uint8_t root[HASH_DIGEST_LENGTH];
        celeres_ihmt_build_root_and_path(commitments, logn, (int32_t)signer_index,
                                          root, paths + (size_t)j * path_bytes);

        /* SIGN.5.iv-v: tau_j=CF(T_{i*} Qbar_j), h_j=O(rt_j,tau_j). */
        generator_mat_t tau = {0};
        if (compute_cf_from_rref_and_monomial(&tau, signer_tag_rref, &qbar) != 0) {
            free(round_hashes);
            free(commitments);
            free(paths);
            free(signer_bits);
            free(signer_perms);
            return -1;
        }
        hash_tau(round_hashes + (size_t)j * HASH_DIGEST_LENGTH, root, &tau);

        /* SIGN.7.i: P_j=PermRep(Q_{i*} Qbar_j), stored if e_j=1. */
        compose_monomial(&product, &private_q, &qbar);
        memcpy(signer_perms[j], product.permutation, sizeof(product.permutation));
    }

    /* SIGN.6: derive the fixed-weight challenge from the transcript. */
    uint8_t digest[HASH_DIGEST_LENGTH];
    fs_digest(digest, params, ring, tags[signer_index].tag,
              message, message_len, salt, round_hashes);
    memcpy(stored_digest, digest, HASH_DIGEST_LENGTH);
    celeres_sample_challenge(challenge, digest);

    /* SIGN.8: release seeds for all rounds with e_j=0. */
    const uint32_t seed_count = GGMPath(seed_tree, challenge, seed_storage);
    if (seed_count > MAX_PUBLISHED_SEEDS) {
        free(round_hashes);
        free(commitments);
        free(paths);
        free(signer_bits);
        free(signer_perms);
        return -1;
    }
    *seed_count_ptr = (uint8_t)seed_count;
    uint8_t *responses = sig_responses(signature, seed_count);

    /* SIGN.7: emit rsp_j=(P_j, Path_j, bits_j) for every e_j=1. */
    uint32_t emitted = 0;
    for (uint32_t j = 0; j < T; j++) {
        if (challenge[j] == 0) {
            continue;
        }
        uint8_t *rsp = responses + (size_t)emitted * celeres_response_bytes();
        permutation_to_bytes(rsp, signer_perms[j]);
        memcpy(rsp + CELERES_PERM_BYTES, signer_bits[j], SEED_LENGTH_BYTES);
        if (path_bytes != 0) {
            memcpy(rsp + CELERES_PERM_BYTES + SEED_LENGTH_BYTES,
                   paths + (size_t)j * path_bytes, path_bytes);
        }
        emitted++;
    }

    *signature_len = celeres_compact_bytes(seed_count);
    free(round_hashes);
    free(commitments);
    free(paths);
    free(signer_bits);
    free(signer_perms);
    return emitted == W ? 0 : -1;
}

static int reconstruct_hashes(uint8_t *round_hashes,
                              uint32_t *opened_count,
                              const celeres_public_parameters_t *params,
                              const celeres_public_key_t *ring,
                              const uint8_t *signature,
                              size_t signature_len) {
    const int logn = celeres_log_round_up_u32(CELERES_RING_SIZE);
    const uint32_t rounded_ring = (uint32_t)1 << logn;

    /* VERIFY.1: parse sigma=(T', salt, e, seeds, responses). */
    if (signature_len < celeres_base_bytes()) {
        return -1;
    }

    const uint8_t *stored_digest = csig_digest(signature);
    const uint8_t *salt = csig_salt(signature);
    const uint32_t seed_count = *csig_seed_count(signature);
    if (seed_count > MAX_PUBLISHED_SEEDS ||
        signature_len != celeres_compact_bytes(seed_count)) {
        return -1;
    }

    uint8_t challenge[T];
    celeres_sample_challenge(challenge, stored_digest);

    /* VERIFY.1 continued: reject if e is not in C_{M,K}. */
    if (!challenge_is_valid(challenge) ||
        seed_count != expected_published_seed_count(challenge)) {
        return -1;
    }

    /* VERIFY.2: recover seed_j for all rounds with e_j=0. */
    uint8_t seed_tree[NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES] = {0};
    uint8_t round_seeds[T * SEED_LENGTH_BYTES] = {0};
    if (!RebuildGGM(seed_tree, challenge, csig_seeds(signature), salt)) {
        return -1;
    }
    seed_leaves(round_seeds, seed_tree);

    uint8_t *commitments = malloc((size_t)rounded_ring * HASH_DIGEST_LENGTH);
    if (commitments == NULL) {
        return -1;
    }

    const uint8_t *tag = csig_tag(signature);
    const uint8_t *responses = csig_responses(signature, seed_count);
    uint8_t base_g_rref[RREF_MAT_PACKEDBYTES];
    uint8_t base_t_rref[RREF_MAT_PACKEDBYTES];
    uint8_t tag_rref[RREF_MAT_PACKEDBYTES];
    public_seed_rref_bytes(base_g_rref, params->g_seed);
#if CELERES_TAG_SELF_ORTHOGONAL
    celeres_tag_sample_base_rref(base_t_rref, params->t_seed);
    if (celeres_tag_expand_to_rref(tag_rref, tag) != 0) {
        free(commitments);
        return -1;
    }
#else
    public_seed_rref_bytes(base_t_rref, params->t_seed);
    memcpy(tag_rref, tag, RREF_MAT_PACKEDBYTES);
#endif
    uint32_t employed = 0;
    monomial_t qbar;

    for (uint32_t j = 0; j < T; j++) {
        if (challenge[j] == 0) {
            /* VERIFY.3.i.a: recompute Qbar_j and r_{j,i}. */
            const uint8_t *seed_j = round_seeds + (size_t)j * SEED_LENGTH_BYTES;
            monomial_sample_salt(&qbar, seed_j, salt, (uint16_t)j);
            for (uint32_t i = 0; i < CELERES_RING_SIZE; i++) {
                uint8_t bits[SEED_LENGTH_BYTES];
                celeres_derive_bits(bits, seed_j, salt, (uint16_t)j, i);

                /* VERIFY.3.i.b: h'_{j,i}=com(r_{j,i}, CF(G_i Qbar_j)). */
                if (compute_commitment_from_pk_with_reuse(
                        commitments + (size_t)i * HASH_DIGEST_LENGTH,
                        &ring[i], &qbar, bits) != 0) {
                    free(commitments);
                    return -1;
                }
            }
            for (uint32_t i = CELERES_RING_SIZE; i < rounded_ring; i++) {
                celeres_dummy_commitment(commitments + (size_t)i * HASH_DIGEST_LENGTH,
                                          salt, (uint16_t)j, i);
            }

            /* VERIFY.3.i.c: construct rt'_j from all recomputed leaves. */
            uint8_t root[HASH_DIGEST_LENGTH];
            celeres_ihmt_build_root_and_path(commitments, logn, -1, root, NULL);

            /* VERIFY.3.i.d-e: tau'_j=CF(T' Qbar_j), h'_j=O(rt'_j,tau'_j). */
            generator_mat_t tau = {0};
            if (compute_cf_from_rref_and_monomial_with_reuse(&tau, tag_rref, &qbar) != 0) {
                free(commitments);
                return -1;
            }
            hash_tau(round_hashes + (size_t)j * HASH_DIGEST_LENGTH, root, &tau);
        } else {
            /* VERIFY.3.ii.a: parse rsp_j=(P_j, Path_j, bits_j). */
            const uint8_t *rsp = responses + (size_t)employed * celeres_response_bytes();
            POSITION_T permutation[N];
            uint8_t leaf[HASH_DIGEST_LENGTH];
            uint8_t root[HASH_DIGEST_LENGTH];
            generator_mat_t tau = {0};
            if (bytes_to_permutation(permutation, rsp) != 0) {
                free(commitments);
                return -1;
            }

            /* VERIFY.3.ii.b-c: h''_j=com(bits_j, CF(G P_j)). */
            if (compute_commitment_from_perm(leaf, base_g_rref, permutation,
                                             rsp + CELERES_PERM_BYTES) != 0) {
                free(commitments);
                return -1;
            }

            /* VERIFY.3.ii.d: reconstruct rt'_j from h''_j and Path_j. */
            celeres_ihmt_reconstruct_root(
                leaf,
                rsp + CELERES_PERM_BYTES + SEED_LENGTH_BYTES,
                logn, root);

            /* VERIFY.3.ii.e-f: tau'_j=CF(T P_j), h'_j=O(rt'_j,tau'_j). */
            if (compute_tau_from_perm(&tau, base_t_rref, permutation) != 0) {
                free(commitments);
                return -1;
            }
            hash_tau(round_hashes + (size_t)j * HASH_DIGEST_LENGTH, root, &tau);
            employed++;
        }
    }

    free(commitments);
    if (opened_count != NULL) {
        *opened_count = employed;
    }
    return employed == W ? 0 : -1;
}

int celeres_verify(const celeres_public_parameters_t *params,
                        const celeres_public_key_t *ring,
                        const uint8_t *message,
                        size_t message_len,
                        const uint8_t *signature,
                        size_t signature_len) {
    if (params == NULL || ring == NULL || signature == NULL) {
        return -1;
    }

    uint8_t *round_hashes = malloc((size_t)T * HASH_DIGEST_LENGTH);
    if (round_hashes == NULL) {
        return -1;
    }

    uint32_t opened = 0;
    if (reconstruct_hashes(round_hashes, &opened, params, ring,
                           signature, signature_len) != 0) {
        free(round_hashes);
        return -1;
    }

    /* VERIFY.4: recompute e'=H_FS(salt,m,pk_1,...,pk_r,T',h'_1,...,h'_M). */
    uint8_t digest[HASH_DIGEST_LENGTH];
    fs_digest(digest, params, ring, csig_tag(signature),
              message, message_len, csig_salt(signature), round_hashes);
    free(round_hashes);

    if (opened != W) {
        return -1;
    }

    /* VERIFY.5: accept exactly when e'=e.  The digest encodes e in sigma. */
    return verify(digest, csig_digest(signature), HASH_DIGEST_LENGTH) == 0 ? 0 : -1;
}

int celeres_link(const uint8_t *signature,
                      size_t signature_len,
                      const uint8_t *other_signature,
                      size_t other_signature_len) {
    if (signature == NULL || other_signature == NULL ||
        signature_len < celeres_base_bytes() ||
        other_signature_len < celeres_base_bytes()) {
        return 0;
    }

    generator_mat_t tag_cf = {0};
    generator_mat_t other_tag_cf = {0};
    if (tag_to_canonical_form(&tag_cf, csig_tag(signature)) != 0 ||
        tag_to_canonical_form(&other_tag_cf, csig_tag(other_signature)) != 0) {
        return 0;
    }

    return canonical_forms_equal(&tag_cf, &other_tag_cf) ? 1 : 0;
}
