#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "codes.h"
#include "celeres.h"
#include "fq_arith.h"

#if !CELERES_TAG_SELF_ORTHOGONAL
static int scale_serialized_rref_tag(uint8_t *signature) {
    generator_mat_t tag = {0};
    uint8_t pivots[N_pad] = {0};

    expand_to_rref(&tag, signature, pivots);
    for (uint32_t col = 0; col < N; col++) {
        if (pivots[col]) {
            continue;
        }

        uint8_t nonzero = 0;
        for (uint32_t row = 0; row < K; row++) {
            nonzero |= (uint8_t)(tag.values[row][col] != 0);
        }
        if (!nonzero) {
            continue;
        }

        for (uint32_t row = 0; row < K; row++) {
            tag.values[row][col] = fq_mul(tag.values[row][col], 2);
        }
        compress_rref(signature, &tag, pivots);
        return 0;
    }
    return -1;
}
#endif

int main(void) {
    const uint32_t ring_size = CELERES_RING_SIZE;
    if (ring_size < 2) {
        puts("test_celeres_link: passed");
        return 0;
    }

    const uint8_t msg1[] = "celeres link test one";
    const uint8_t msg2[] = "celeres link test two";

    celeres_public_parameters_t params;
    celeres_public_key_t *ring = calloc(ring_size, sizeof(*ring));
    celeres_tag_t *tags = calloc(ring_size, sizeof(*tags));
    celeres_secret_key_t *sk = calloc(ring_size, sizeof(*sk));
    if (ring == NULL || tags == NULL || sk == NULL) {
        return 1;
    }

    if (celeres_generate_parameters(&params) != 0) {
        return 1;
    }
    for (uint32_t i = 0; i < ring_size; i++) {
        if (celeres_keygen(&params, &ring[i], &tags[i], &sk[i]) != 0) {
            return 1;
        }
    }

    const size_t cap = celeres_signature_bytes();
    uint8_t *sig0 = malloc(cap);
    uint8_t *sig1 = malloc(cap);
    uint8_t *sig_other = malloc(cap);
    if (sig0 == NULL || sig1 == NULL || sig_other == NULL) {
        return 1;
    }

    size_t len0 = 0;
    size_t len1 = 0;
    size_t len_other = 0;
    if (celeres_sign(&params, ring, tags, 0, &sk[0],
                          msg1, sizeof(msg1) - 1, sig0, &len0) != 0 ||
        celeres_sign(&params, ring, tags, 0, &sk[0],
                          msg2, sizeof(msg2) - 1, sig1, &len1) != 0 ||
        celeres_sign(&params, ring, tags, 1, &sk[1],
                          msg1, sizeof(msg1) - 1, sig_other, &len_other) != 0) {
        fprintf(stderr, "sign failed\n");
        return 1;
    }

    if (!celeres_link(sig0, len0, sig1, len1)) {
        fprintf(stderr, "same signer signatures did not link\n");
        return 1;
    }
#if !CELERES_TAG_SELF_ORTHOGONAL
    uint8_t *sig_scaled = malloc(cap);
    if (sig_scaled == NULL) {
        return 1;
    }
    memcpy(sig_scaled, sig0, len0);
    if (scale_serialized_rref_tag(sig_scaled) != 0) {
        fprintf(stderr, "could not build equivalent non-identical tag\n");
        return 1;
    }
    if (memcmp(sig0, sig_scaled, CELERES_TAG_BYTES) == 0) {
        fprintf(stderr, "scaled tag did not change serialization\n");
        return 1;
    }
    if (!celeres_link(sig0, len0, sig_scaled, len0)) {
        fprintf(stderr, "canonically equivalent tags did not link\n");
        return 1;
    }
    free(sig_scaled);
#endif
    if (celeres_link(sig0, len0, sig_other, len_other)) {
        fprintf(stderr, "different signer signatures linked\n");
        return 1;
    }

    free(sig0);
    free(sig1);
    free(sig_other);
    free(ring);
    free(tags);
    free(sk);
    puts("test_celeres_link: passed");
    return 0;
}
