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

int celeres_sf_from_rref_monomial(
    generator_mat_t *out,
    uint8_t out_pivots[N_pad],
    const uint8_t rref[RREF_MAT_PACKEDBYTES],
    const monomial_t *q) {
    generator_mat_t base = {0};
    generator_mat_t product = {0};
    uint8_t base_pivots[N_pad] = {0};
    POSITION_T pivot_rank_by_source_col[N];
    POSITION_T source_by_destination_col[N];
    POSITION_T pivot_rows[K];
    POSITION_T pivot_cols[K];
    uint16_t pivot_count = 0;
    uint16_t ordered_count = 0;

    expand_to_rref(&base, rref, base_pivots);
    generator_monomial_mul(&product, &base, q);

    memset(out, 0, sizeof(*out));
    memset(out_pivots, 0, N_pad);
    for (uint16_t col = 0; col < N; col++) {
        pivot_rank_by_source_col[col] = 0;
        source_by_destination_col[col] = 0;
    }

    for (uint16_t src_col = 0; src_col < N; src_col++) {
        const POSITION_T dst_col = q->permutation[src_col];
        if (dst_col >= N) {
            return -1;
        }
        source_by_destination_col[dst_col] = src_col;
        if (base_pivots[src_col] != 0) {
            if (pivot_count >= K) {
                return -1;
            }
            pivot_rank_by_source_col[src_col] = pivot_count++;
        }
    }
    if (pivot_count != K) {
        return -1;
    }

    for (uint16_t dst_col = 0; dst_col < N; dst_col++) {
        const POSITION_T src_col = source_by_destination_col[dst_col];
        if (base_pivots[src_col] != 0) {
            if (ordered_count >= K) {
                return -1;
            }
            pivot_rows[ordered_count] = pivot_rank_by_source_col[src_col];
            pivot_cols[ordered_count] = dst_col;
            out_pivots[dst_col] = 1;
            ordered_count++;
        }
    }
    if (ordered_count != K) {
        return -1;
    }

    for (uint16_t row = 0; row < K; row++) {
        const POSITION_T src_row = pivot_rows[row];
        const POSITION_T pivot_col = pivot_cols[row];
        const FQ_ELEM pivot = product.values[src_row][pivot_col];
        if (pivot == 0) {
            return -1;
        }
        const FQ_ELEM scale = fq_inv(pivot);
        for (uint16_t col = 0; col < N; col++) {
            out->values[row][col] = fq_mul(product.values[src_row][col], scale);
        }
        for (uint16_t col = N; col < N_pad; col++) {
            out->values[row][col] = 0;
        }
    }

    return 0;
}

int celeres_cf_from_rref_monomial_is(
    generator_mat_t *cf,
    const uint8_t rref[RREF_MAT_PACKEDBYTES],
    const monomial_t *q,
    int use_ct_scaling) {
    generator_mat_t systematic = {0};
    uint8_t pivots[N_pad] = {0};
    if (celeres_sf_from_rref_monomial(&systematic, pivots, rref, q) != 0) {
        return -1;
    }
    if (use_ct_scaling) {
        celeres_diagonal_scaling_canonical_form_ct(&systematic);
    } else {
        celeres_diagonal_scaling_canonical_form(&systematic);
    }
    *cf = systematic;
    return 0;
}
