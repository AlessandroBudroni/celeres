#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "celeres.h"

#if defined(__has_builtin)
#if __has_builtin(__builtin_readcyclecounter)
#define HAS_BUILTIN_CYCLECOUNTER 1
#endif
#endif

#if !defined(HAS_BUILTIN_CYCLECOUNTER) && \
    (defined(__i386__) || defined(__x86_64__))
#include <x86intrin.h>
#define HAS_RDTSC 1
#endif

#ifndef IMPL_NAME
#define IMPL_NAME "unknown"
#endif

#ifndef BENCH_ITERS
#define BENCH_ITERS 100
#endif

static double now_ms(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec * 1000.0 + (double)tv.tv_usec / 1000.0;
}

static uint64_t now_cycles(void) {
#if defined(HAS_BUILTIN_CYCLECOUNTER)
    return __builtin_readcyclecounter();
#elif defined(HAS_RDTSC)
    return __rdtsc();
#else
    return 0;
#endif
}

static int cycles_available(void) {
#if defined(HAS_BUILTIN_CYCLECOUNTER) || defined(HAS_RDTSC)
    return 1;
#else
    return 0;
#endif
}

static void print_result(const char *label,
                         uint32_t iters,
                         double elapsed_ms,
                         uint64_t elapsed_cycles) {
    printf("%s_avg_ms=%.3f\n", label, elapsed_ms / (double)iters);
    if (cycles_available()) {
        printf("%s_avg_cycles=%.2f M\n", label,
               ((double)elapsed_cycles / (double)iters) / 1000000.0);
    }
}

static void cleanup(uint8_t *signatures,
                    size_t *signature_lengths,
                    celeres_public_parameters_t *keygen_params,
                    celeres_public_key_t *ring,
                    celeres_tag_t *tags,
                    celeres_secret_key_t *sk) {
    free(signatures);
    free(signature_lengths);
    free(keygen_params);
    free(ring);
    free(tags);
    free(sk);
}

int main(void) {
    const uint32_t ring_size = CELERES_RING_SIZE;
    const uint32_t signer = ring_size > 1 ? 1u : 0u;
    const uint8_t msg[] = "celeres benchmark message";

    celeres_public_parameters_t params;
    celeres_public_key_t *ring = calloc(ring_size, sizeof(*ring));
    celeres_tag_t *tags = calloc(ring_size, sizeof(*tags));
    celeres_secret_key_t *sk = calloc(ring_size, sizeof(*sk));
    if (ring == NULL || tags == NULL || sk == NULL || ring_size == 0) {
        cleanup(NULL, NULL, NULL, ring, tags, sk);
        return 1;
    }

    if (celeres_generate_parameters(&params) != 0) {
        cleanup(NULL, NULL, NULL, ring, tags, sk);
        return 1;
    }

    const size_t sig_cap = celeres_signature_bytes();
    uint8_t *signatures = calloc(BENCH_ITERS, sig_cap);
    size_t *signature_lengths = calloc(BENCH_ITERS,
                                       sizeof(*signature_lengths));
    celeres_public_parameters_t *keygen_params =
        calloc(BENCH_ITERS, sizeof(*keygen_params));
    if (sig_cap == 0 || signatures == NULL ||
        signature_lengths == NULL || keygen_params == NULL) {
        cleanup(signatures, signature_lengths, keygen_params, ring, tags, sk);
        return 1;
    }

    for (uint32_t i = 0; i < BENCH_ITERS; i++) {
        if (celeres_generate_parameters(&keygen_params[i]) != 0) {
            cleanup(signatures, signature_lengths, keygen_params, ring, tags, sk);
            return 1;
        }
    }

    double start_ms = now_ms();
    uint64_t start_cycles = now_cycles();
    for (uint32_t i = 0; i < BENCH_ITERS; i++) {
        celeres_public_key_t pk_tmp;
        celeres_tag_t tag_tmp;
        celeres_secret_key_t sk_tmp;
        if (celeres_keygen(&keygen_params[i], &pk_tmp, &tag_tmp, &sk_tmp) != 0) {
            cleanup(signatures, signature_lengths, keygen_params, ring, tags, sk);
            return 1;
        }
    }
    const uint64_t keygen_cycles = now_cycles() - start_cycles;
    const double keygen_ms = now_ms() - start_ms;

    for (uint32_t i = 0; i < ring_size; i++) {
        if (celeres_keygen(&params, &ring[i], &tags[i], &sk[i]) != 0) {
            cleanup(signatures, signature_lengths, keygen_params, ring, tags, sk);
            return 1;
        }
    }

    size_t total_sig_bytes = 0;
    start_ms = now_ms();
    start_cycles = now_cycles();
    for (uint32_t i = 0; i < BENCH_ITERS; i++) {
        uint8_t *sig = signatures + (size_t)i * sig_cap;
        if (celeres_sign(&params, ring, tags, signer, &sk[signer],
                              msg, sizeof(msg) - 1, sig,
                              &signature_lengths[i]) != 0) {
            cleanup(signatures, signature_lengths, keygen_params, ring, tags, sk);
            return 1;
        }
        total_sig_bytes += signature_lengths[i];
    }
    const uint64_t sign_cycles = now_cycles() - start_cycles;
    const double sign_ms = now_ms() - start_ms;

    for (uint32_t i = 0; i < BENCH_ITERS; i++) {
        const uint8_t *sig = signatures + (size_t)i * sig_cap;
        if (celeres_verify(&params, ring, msg, sizeof(msg) - 1,
                                sig, signature_lengths[i]) != 0) {
            cleanup(signatures, signature_lengths, keygen_params, ring, tags, sk);
            return 1;
        }
    }

    start_ms = now_ms();
    start_cycles = now_cycles();
    for (uint32_t i = 0; i < BENCH_ITERS; i++) {
        const uint8_t *sig = signatures + (size_t)i * sig_cap;
        if (celeres_verify(&params, ring, msg, sizeof(msg) - 1,
                                sig, signature_lengths[i]) != 0) {
            cleanup(signatures, signature_lengths, keygen_params, ring, tags, sk);
            return 1;
        }
    }
    const uint64_t verify_cycles = now_cycles() - start_cycles;
    const double verify_ms = now_ms() - start_ms;

    start_ms = now_ms();
    start_cycles = now_cycles();
    for (uint32_t i = 1; i < BENCH_ITERS; i++) {
        const uint8_t *a = signatures + (size_t)(i - 1u) * sig_cap;
        const uint8_t *b = signatures + (size_t)i * sig_cap;
        (void)celeres_link(a, signature_lengths[i - 1u],
                                b, signature_lengths[i]);
    }
    const uint64_t link_cycles = now_cycles() - start_cycles;
    const double link_ms = now_ms() - start_ms;
    const uint32_t link_iters =
        BENCH_ITERS > 1u ? BENCH_ITERS - 1u : 1u;

    printf("impl=%s\n", IMPL_NAME);
    printf("CATEGORY=%d\n", CATEGORY);
    printf("TARGET=%d\n", TARGET);
    printf("N=%d\n", N);
    printf("K=%d\n", K);
    printf("M/T=%d\n", T);
    printf("challenge_weight/W=%d\n", W);
    printf("ring_size=%u\n", ring_size);
    printf("celeres_tag_mode=%s\n", CELERES_TAG_MODE_NAME);
    printf("celeres_tag_compression=%s\n", CELERES_TAG_COMPRESSION_NAME);
    printf("benchmark_iterations=%u\n", BENCH_ITERS);
    printf("public_key_bytes=%zu\n", sizeof(celeres_public_key_t));
    printf("tag_bytes=%zu\n", sizeof(celeres_tag_t));
    printf("secret_key_bytes=%zu\n", sizeof(celeres_secret_key_t));
    printf("permutation_response_bytes=%zu\n", (size_t)CELERES_PERM_BYTES);
    printf("avg_signature_bytes=%.3f\n",
           (double)total_sig_bytes / (double)BENCH_ITERS);
    print_result("keygen", BENCH_ITERS, keygen_ms, keygen_cycles);
    print_result("sign", BENCH_ITERS, sign_ms, sign_cycles);
    print_result("verify", BENCH_ITERS, verify_ms, verify_cycles);
    print_result("link", link_iters, link_ms, link_cycles);

    cleanup(signatures, signature_lengths, keygen_params, ring, tags, sk);
    return 0;
}
