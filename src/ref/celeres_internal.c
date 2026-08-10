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

#ifndef CELERES_CT_PROPAGATION_PASSES
#define CELERES_CT_PROPAGATION_PASSES (K + N)
#endif

static FQ_ELEM ct_select_fq(FQ_ELEM old_value,
                            FQ_ELEM new_value,
                            uint8_t take) {
    const uint8_t mask = (uint8_t)(0u - (uint8_t)(take & 1u));
    return (FQ_ELEM)((old_value & (uint8_t)~mask) | (new_value & mask));
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
    FQ_ELEM beta[N];
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

    for (uint16_t a = 0; a < Q; a++) {
        for (uint16_t col = 0; col < N; col++) {
            scaled_beta[a][col] = fq_mul((FQ_ELEM)a, beta[col]);
        }
        for (uint16_t col = N; col < N_pad; col++) {
            scaled_beta[a][col] = 0;
        }
    }

    for (uint16_t row = 0; row < K; row++) {
        const FQ_ELEM *row_scale = scaled_beta[alpha[row]];
        for (uint16_t col = 0; col < N_pad; col++) {
            g->values[row][col] = fq_mul(g->values[row][col], row_scale[col]);
        }
    }
}

void celeres_diagonal_scaling_canonical_form_ct(generator_mat_t *g) {
    FQ_ELEM alpha[K];
    FQ_ELEM beta[N];
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

    for (uint16_t a = 0; a < Q; a++) {
        for (uint16_t col = 0; col < N; col++) {
            scaled_beta[a][col] = fq_mul((FQ_ELEM)a, beta[col]);
        }
        for (uint16_t col = N; col < N_pad; col++) {
            scaled_beta[a][col] = 0;
        }
    }

    for (uint16_t row = 0; row < K; row++) {
        const FQ_ELEM *row_scale = scaled_beta[alpha[row]];
        for (uint16_t col = 0; col < N_pad; col++) {
            g->values[row][col] = fq_mul(g->values[row][col], row_scale[col]);
        }
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
