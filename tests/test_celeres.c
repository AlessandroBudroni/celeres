#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "celeres.h"

int main(void) {
    const uint32_t ring_size = CELERES_RING_SIZE;
    const uint32_t signer = ring_size > 1 ? 1u : 0u;
    const uint8_t msg[] = "celeres protocol test";

    celeres_public_parameters_t params;
    celeres_public_key_t *ring = calloc(ring_size, sizeof(*ring));
    celeres_tag_t *tags = calloc(ring_size, sizeof(*tags));
    celeres_secret_key_t *sk = calloc(ring_size, sizeof(*sk));
    if (ring == NULL || tags == NULL || sk == NULL) {
        free(ring);
        free(tags);
        free(sk);
        return 1;
    }

    if (celeres_generate_parameters(&params) != 0) {
        return 1;
    }
    for (uint32_t i = 0; i < ring_size; i++) {
        if (celeres_keygen(&params, &ring[i], &tags[i], &sk[i]) != 0) {
            fprintf(stderr, "keygen failed at %u\n", i);
            return 1;
        }
    }

    const size_t sig_cap = celeres_signature_bytes();
    uint8_t *sig = malloc(sig_cap);
    if (sig == NULL) {
        return 1;
    }

    size_t sig_len = 0;
    if (celeres_sign(&params, ring, tags, signer, &sk[signer],
                          msg, sizeof(msg) - 1, sig, &sig_len) != 0) {
        fprintf(stderr, "sign failed\n");
        return 1;
    }
    if (sig_len == 0 || sig_len > sig_cap) {
        fprintf(stderr, "bad signature size: %zu > %zu\n", sig_len, sig_cap);
        return 1;
    }
    if (celeres_verify(&params, ring, msg, sizeof(msg) - 1,
                            sig, sig_len) != 0) {
        fprintf(stderr, "verify failed\n");
        return 1;
    }

    uint8_t bad_msg[] = "celeres protocol tesu";
    if (celeres_verify(&params, ring, bad_msg, sizeof(bad_msg) - 1,
                            sig, sig_len) == 0) {
        fprintf(stderr, "tampered message verified\n");
        return 1;
    }

    free(sig);
    free(ring);
    free(tags);
    free(sk);
    puts("test_celeres: passed");
    return 0;
}
