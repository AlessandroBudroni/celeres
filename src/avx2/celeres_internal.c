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

#include "celeres_internal.h"

#include <string.h>

#include "fq_arith.h"
#include "macro.h"

#ifndef CELERES_CT_PROPAGATION_PASSES
#define CELERES_CT_PROPAGATION_PASSES (K + N)
#endif

static FQ_ELEM ct_select_fq(FQ_ELEM old_value,
                            FQ_ELEM new_value,
                            uint8_t take) {
    const uint8_t mask = (uint8_t)(0u - (uint8_t)(take & 1u));
    return (FQ_ELEM)((old_value & (uint8_t)~mask) | (new_value & mask));
}

static vec256_t vec_mul_scalar_mod_q(vec256_t a, FQ_ELEM scalar) {
    vec256_t shuffle;
    vec256_t t;
    vec256_t c01;
    vec256_t c7f;
    vec256_t b;
    vec256_t a_lo;
    vec256_t a_hi;
    vec256_t b_lo;
    vec256_t b_hi;
    vec128_t tmp;

    vload256(shuffle, (const vec256_t *)shuff_low_half);
    vset8(c01, 1);
    vset8(c7f, 127);
    vset8(b, scalar);

    vget_lo(tmp, a);
    vextend8_16(a_lo, tmp);
    vget_hi(tmp, a);
    vextend8_16(a_hi, tmp);
    vget_lo(tmp, b);
    vextend8_16(b_lo, tmp);
    vget_hi(tmp, b);
    vextend8_16(b_hi, tmp);

    barrett_mul_u16(a_lo, a_lo, b_lo, t);
    barrett_mul_u16(a_hi, a_hi, b_hi, t);
    vshuffle8(a_lo, a_lo, shuffle);
    vshuffle8(a_hi, a_hi, shuffle);
    vpermute_4x64(a_lo, a_lo, 0xd8);
    vpermute_4x64(a_hi, a_hi, 0xd8);
    vpermute2(t, a_lo, a_hi, 0x20);
    W_RED127_(t);
    return t;
}

static vec256_t vec_mul_pair_mod_q(vec256_t a, vec256_t b) {
    vec256_t shuffle;
    vec256_t t;
    vec256_t c01;
    vec256_t c7f;
    vec256_t a_lo;
    vec256_t a_hi;
    vec256_t b_lo;
    vec256_t b_hi;
    vec128_t tmp;

    vload256(shuffle, (const vec256_t *)shuff_low_half);
    vset8(c01, 1);
    vset8(c7f, 127);

    vget_lo(tmp, a);
    vextend8_16(a_lo, tmp);
    vget_hi(tmp, a);
    vextend8_16(a_hi, tmp);
    vget_lo(tmp, b);
    vextend8_16(b_lo, tmp);
    vget_hi(tmp, b);
    vextend8_16(b_hi, tmp);

    barrett_mul_u16(a_lo, a_lo, b_lo, t);
    barrett_mul_u16(a_hi, a_hi, b_hi, t);
    vshuffle8(a_lo, a_lo, shuffle);
    vshuffle8(a_hi, a_hi, shuffle);
    vpermute_4x64(a_lo, a_lo, 0xd8);
    vpermute_4x64(a_hi, a_hi, 0xd8);
    vpermute2(t, a_lo, a_hi, 0x20);
    W_RED127_(t);
    return t;
}

static void build_scaled_beta_table(FQ_ELEM scaled_beta[Q][N_pad],
                                    const FQ_ELEM beta[N_pad]) {
    for (uint16_t a = 0; a < Q; a++) {
        for (uint16_t col = 0; col < N_pad; col += 32) {
            vec256_t beta_v;
            vload256(beta_v, (const vec256_t *)(beta + col));
            beta_v = vec_mul_scalar_mod_q(beta_v, (FQ_ELEM)a);
            vstore256((vec256_t *)(scaled_beta[a] + col), beta_v);
        }
    }
}

static void row_scale_by_beta_row(FQ_ELEM row[N_pad],
                                  const FQ_ELEM beta_row[N_pad]) {
    for (uint16_t col = 0; col < N_pad; col += 32) {
        vec256_t row_v;
        vec256_t beta_v;
        vload256(row_v, (const vec256_t *)(row + col));
        vload256(beta_v, (const vec256_t *)(beta_row + col));
        row_v = vec_mul_pair_mod_q(row_v, beta_v);
        vstore256((vec256_t *)(row + col), row_v);
    }
}

int celeres_rref(generator_mat_t *g,
                      uint8_t is_pivot_column[N_pad]) {
    memset(is_pivot_column, 0, N_pad);
    return generator_RREF(g, is_pivot_column) == 1 ? 0 : -1;
}

int celeres_rref_with_hint(generator_mat_t *g,
                                uint8_t is_pivot_column[N_pad],
                                const uint8_t pivot_hint[N_pad]) {
    uint8_t hint[N_pad] = {0};
    memset(is_pivot_column, 0, N_pad);
    memcpy(hint, pivot_hint, N);
    return generator_RREF_pivot_reuse(g, is_pivot_column, hint, K) == 1 ? 0 : -1;
}

void celeres_diagonal_scaling_canonical_form(generator_mat_t *g) {
    FQ_ELEM alpha[K];
    FQ_ELEM beta[N_pad];
    FQ_ELEM scaled_beta[Q][N_pad];
    uint8_t touched_rows[K] = {0};
    uint8_t touched_cols[N] = {0};
    uint16_t touched_row_count = 0;
    uint16_t touched_col_count = 0;

    for (uint16_t row = 0; row < K; row++) {
        alpha[row] = 1;
    }
    for (uint16_t col = 0; col < N; col++) {
        beta[col] = 1;
    }
    for (uint16_t col = N; col < N_pad; col++) {
        beta[col] = 0;
    }

    for (uint16_t root = 0; root < K; root++) {
        if (touched_row_count == K) {
            break;
        }
        if (touched_rows[root] != 0) {
            continue;
        }

        uint8_t row_nonzero = 0;
        for (uint16_t col = 0; col < N; col++) {
            row_nonzero |= (uint8_t)(g->values[root][col] != 0);
        }
        if (row_nonzero == 0) {
            continue;
        }

        uint16_t row_queue[K];
        uint16_t col_queue[N];
        uint16_t row_head = 0;
        uint16_t row_tail = 0;
        uint16_t col_head = 0;
        uint16_t col_tail = 0;

        alpha[root] = 1;
        touched_rows[root] = 1;
        touched_row_count++;
        row_queue[row_tail++] = root;

        while (row_head < row_tail || col_head < col_tail) {
            while (row_head < row_tail) {
                if (touched_col_count == N) {
                    row_head = row_tail;
                    break;
                }
                const uint16_t row = row_queue[row_head++];
                for (uint16_t col = 0; col < N; col++) {
                    if (g->values[row][col] != 0 && touched_cols[col] == 0) {
                        beta[col] = fq_inv(fq_mul(alpha[row], g->values[row][col]));
                        touched_cols[col] = 1;
                        touched_col_count++;
                        col_queue[col_tail++] = col;
                    }
                }
            }

            while (col_head < col_tail) {
                if (touched_row_count == K) {
                    col_head = col_tail;
                    break;
                }
                const uint16_t col = col_queue[col_head++];
                for (uint16_t row = 0; row < K; row++) {
                    if (g->values[row][col] != 0 && touched_rows[row] == 0) {
                        alpha[row] = fq_inv(fq_mul(g->values[row][col], beta[col]));
                        touched_rows[row] = 1;
                        touched_row_count++;
                        row_queue[row_tail++] = row;
                    }
                }
            }
        }
    }

    build_scaled_beta_table(scaled_beta, beta);
    for (uint16_t row = 0; row < K; row++) {
        row_scale_by_beta_row(g->values[row], scaled_beta[alpha[row]]);
    }
}

void celeres_diagonal_scaling_canonical_form_ct(generator_mat_t *g) {
    FQ_ELEM alpha[K];
    FQ_ELEM beta[N_pad];
    FQ_ELEM scaled_beta[Q][N_pad];
    uint8_t touched_rows[K] = {0};
    uint8_t touched_cols[N] = {0};
    uint16_t touched_row_count = 0;
    uint16_t touched_col_count = 0;

    for (uint16_t row = 0; row < K; row++) {
        alpha[row] = 1;
    }
    for (uint16_t col = 0; col < N; col++) {
        beta[col] = 1;
    }
    for (uint16_t col = N; col < N_pad; col++) {
        beta[col] = 0;
    }

    for (uint16_t root = 0; root < K; root++) {
        uint8_t row_nonzero = 0;
        for (uint16_t col = 0; col < N; col++) {
            row_nonzero |= (uint8_t)(g->values[root][col] != 0);
        }

        const uint8_t start_component =
            (uint8_t)((touched_rows[root] == 0) & row_nonzero);
        if (start_component == 0) {
            continue;
        }
        touched_rows[root] |= start_component;
        touched_row_count++;
        alpha[root] = ct_select_fq(alpha[root], 1, start_component);

        uint16_t row_queue[K];
        uint16_t col_queue[N];
        uint16_t row_head = 0;
        uint16_t row_tail = 0;
        uint16_t col_head = 0;
        uint16_t col_tail = 0;

        row_queue[row_tail++] = root;
        for (uint16_t pass = 0; pass < CELERES_CT_PROPAGATION_PASSES; pass++) {
            for (uint16_t step = 0; step < K; step++) {
                if (row_head >= row_tail || touched_col_count == N) {
                    continue;
                }
                const uint16_t row = row_queue[row_head++];
                for (uint16_t col = 0; col < N; col++) {
                    if (g->values[row][col] != 0 && touched_cols[col] == 0) {
                        beta[col] = fq_inv(fq_mul(alpha[row], g->values[row][col]));
                        touched_cols[col] = 1;
                        touched_col_count++;
                        col_queue[col_tail++] = col;
                    }
                }
            }

            for (uint16_t step = 0; step < N; step++) {
                if (col_head >= col_tail || touched_row_count == K) {
                    continue;
                }
                const uint16_t col = col_queue[col_head++];
                for (uint16_t row = 0; row < K; row++) {
                    if (g->values[row][col] != 0 && touched_rows[row] == 0) {
                        alpha[row] = fq_inv(fq_mul(g->values[row][col], beta[col]));
                        touched_rows[row] = 1;
                        touched_row_count++;
                        row_queue[row_tail++] = row;
                    }
                }
            }
        }
    }

    build_scaled_beta_table(scaled_beta, beta);
    for (uint16_t row = 0; row < K; row++) {
        row_scale_by_beta_row(g->values[row], scaled_beta[alpha[row]]);
    }
}

int celeres_cf(generator_mat_t *g) {
    uint8_t pivots[N_pad] = {0};
    const int rc = celeres_rref(g, pivots);
    celeres_diagonal_scaling_canonical_form(g);
    return rc;
}

int celeres_cf_sign_ct(generator_mat_t *g) {
    uint8_t pivots[N_pad] = {0};
    const int rc = celeres_rref(g, pivots);
    celeres_diagonal_scaling_canonical_form_ct(g);
    return rc;
}

int celeres_cf_with_hint(generator_mat_t *g,
                              const uint8_t pivot_hint[N_pad]) {
    uint8_t pivots[N_pad] = {0};
    const int rc = celeres_rref_with_hint(g, pivots, pivot_hint);
    celeres_diagonal_scaling_canonical_form(g);
    return rc;
}

void celeres_perm_rep(POSITION_T out[N],
                           const monomial_t *m) {
    for (uint16_t i = 0; i < N; i++) {
        out[i] = m->permutation[i];
    }
}
