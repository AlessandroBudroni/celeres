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

#include "celeres_tag.h"

#include <string.h>

#include "fq_arith.h"
#include "rng.h"

#if CELERES_TAG_SELF_ORTHOGONAL

static const uint8_t fq_sqrt_table[128] = {
    0, 1, 16, 0, 2, 0, 0, 0, 32, 3, 0, 30, 0, 34, 0, 53,
    4, 12, 48, 20, 0, 23, 28, 0, 0, 5, 36, 0, 0, 0, 41, 44,
    63, 0, 62, 17, 6, 52, 61, 0, 0, 26, 13, 0, 60, 0, 0, 38,
    0, 7, 47, 0, 59, 0, 0, 0, 0, 0, 0, 0, 21, 51, 58, 0,
    8, 0, 0, 0, 24, 14, 18, 43, 31, 33, 57, 0, 40, 0, 0, 29,
    0, 9, 35, 0, 46, 0, 0, 50, 56, 0, 0, 0, 0, 0, 27, 0,
    0, 0, 15, 37, 10, 0, 0, 22, 55, 0, 0, 19, 0, 0, 0, 0,
    0, 42, 0, 49, 0, 25, 0, 0, 45, 11, 54, 0, 39, 0, 0
};

static FQ_ELEM fq_sqrt(FQ_ELEM x) {
    return fq_sqrt_table[x];
}

static void bit_write(uint8_t *out, size_t bit_pos, uint32_t value, uint32_t bits) {
    for (uint32_t i = 0; i < bits; i++) {
        const size_t pos = bit_pos + i;
        const uint8_t bit = (uint8_t)((value >> i) & 1u);
        out[pos >> 3] = (uint8_t)((out[pos >> 3] & ~(uint8_t)(1u << (pos & 7u))) |
                                  (uint8_t)(bit << (pos & 7u)));
    }
}

static uint32_t bit_read(const uint8_t *in, size_t bit_pos, uint32_t bits) {
    uint32_t value = 0;
    for (uint32_t i = 0; i < bits; i++) {
        const size_t pos = bit_pos + i;
        value |= (uint32_t)((in[pos >> 3] >> (pos & 7u)) & 1u) << i;
    }
    return value;
}

static FQ_ELEM inner_square(const FQ_ELEM *row, uint32_t len) {
    FQ_ELEM out = 0;
    for (uint32_t i = 0; i < len; i++) {
        out = fq_add(out, fq_mul(row[i], row[i]));
    }
    return out;
}

static FQ_ELEM scalar_product(const FQ_ELEM *a, const FQ_ELEM *b, uint32_t len) {
    FQ_ELEM out = 0;
    for (uint32_t i = 0; i < len; i++) {
        out = fq_add(out, fq_mul(a[i], b[i]));
    }
    return out;
}

static void row_mul_len(FQ_ELEM *row, FQ_ELEM scale, uint32_t len) {
    for (uint32_t i = 0; i < len; i++) {
        row[i] = fq_mul(row[i], scale);
    }
}

static void row_sum_len(FQ_ELEM *out,
                        const FQ_ELEM *in,
                        FQ_ELEM scale,
                        uint32_t len) {
    for (uint32_t i = 0; i < len; i++) {
        out[i] = fq_add(out[i], fq_mul(in[i], scale));
    }
}

static int anti_normalize(FQ_ELEM row[K]) {
    const FQ_ELEM norm = inner_square(row, K);
    const FQ_ELEM root = fq_sqrt(fq_inv(fq_sub(0, norm)));
    if (root == 0) {
        return 0;
    }
    row_mul_len(row, root, K);
    return 1;
}

static void row_mat_mult(FQ_ELEM *out,
                         const FQ_ELEM *row,
                         FQ_ELEM matrix[K][K],
                         uint32_t rows,
                         uint32_t cols) {
    for (uint32_t col = 0; col < cols; col++) {
        out[col] = 0;
        for (uint32_t r = 0; r < rows; r++) {
            out[col] = fq_add(out[col], fq_mul(row[r], matrix[r][col]));
        }
    }
}

static void swap_columns(FQ_ELEM matrix[K][K],
                         uint32_t c1,
                         uint32_t c2,
                         uint32_t rows) {
    if (c1 == c2) {
        return;
    }
    for (uint32_t row = 0; row < rows; row++) {
        const FQ_ELEM tmp = matrix[row][c1];
        matrix[row][c1] = matrix[row][c2];
        matrix[row][c2] = tmp;
    }
}

static void sample_antiorthogonal(FQ_ELEM out[K][K],
                                  const uint8_t seed[SEED_LENGTH_BYTES]) {
    SHAKE_STATE_STRUCT csprng_state;
    FQ_ELEM echelon[K][K] = {{0}};
    FQ_ELEM row[K] = {0};
    FQ_ELEM constraints[K][K] = {{0}};

    initialize_csprng(&csprng_state, seed, SEED_LENGTH_BYTES);

    do {
        rand_range_q_state_elements(&csprng_state, row, K);
    } while (row[0] == 0 || !anti_normalize(row));

    memcpy(out[0], row, K * sizeof(FQ_ELEM));
    row_mul_len(row, fq_inv(row[0]), K);
    memcpy(echelon[0], row, K * sizeof(FQ_ELEM));

    for (uint32_t pivot = 1; pivot < K; pivot++) {
        for (uint32_t i = 0; i < K - pivot; i++) {
            for (uint32_t j = 0; j < pivot; j++) {
                constraints[i][j] = fq_sub(0, echelon[j][i + pivot]);
            }
        }

        FQ_ELEM free_part[K] = {0};
        do {
            rand_range_q_state_elements(&csprng_state, free_part, K - pivot);
            row_mat_mult(row, free_part, constraints, K - pivot, pivot);
            memcpy(row + pivot, free_part, (K - pivot) * sizeof(FQ_ELEM));
        } while (!anti_normalize(row));

        memcpy(out[pivot], row, K * sizeof(FQ_ELEM));
        memcpy(echelon[pivot], row, K * sizeof(FQ_ELEM));

        for (uint32_t i = 0; i < pivot; i++) {
            if (echelon[pivot][i] != 0) {
                row_sum_len(echelon[pivot], echelon[i],
                            fq_sub(0, echelon[pivot][i]), K);
            }
        }

        if (echelon[pivot][pivot] == 0) {
            for (uint32_t i = pivot + 1; i < K; i++) {
                if (echelon[pivot][i] != 0) {
                    swap_columns(echelon, pivot, i, pivot + 1);
                    swap_columns(out, pivot, i, pivot + 1);
                    break;
                }
            }
        }

        row_mul_len(echelon[pivot], fq_inv(echelon[pivot][pivot]), K);

        for (uint32_t i = 0; i < pivot; i++) {
            if (echelon[i][pivot] != 0) {
                row_sum_len(echelon[i], echelon[pivot],
                            fq_sub(0, echelon[i][pivot]), K);
            }
        }
    }
}

void celeres_tag_sample_base_rref(uint8_t out[RREF_MAT_PACKEDBYTES],
                                  const uint8_t seed[SEED_LENGTH_BYTES]) {
    rref_generator_mat_t compact = {0};
    generator_mat_t full = {0};
    uint8_t pivots[N_pad] = {0};
    FQ_ELEM anti[K][K];

    sample_antiorthogonal(anti, seed);
    for (uint32_t row = 0; row < K; row++) {
        memcpy(compact.values[row], anti[row], K * sizeof(FQ_ELEM));
    }
    for (uint32_t i = 0; i < N - K; i++) {
        compact.column_pos[i] = (POSITION_T)(i + K);
    }
    generator_rref_expand(&full, &compact);
    for (uint32_t i = 0; i < K; i++) {
        pivots[i] = 1;
    }
    compress_rref(out, &full, pivots);
}

static int choose_column_order(POSITION_T order[K], const FQ_ELEM b[K][K]) {
    FQ_ELEM work[K][K];
    uint8_t used[K] = {0};

    memcpy(work, b, sizeof(work));
    for (uint32_t row = 0; row < K; row++) {
        uint32_t pivot_col = K;
        for (uint32_t col = 0; col < K; col++) {
            if (!used[col] && work[row][col] != 0) {
                pivot_col = col;
                break;
            }
        }
        if (pivot_col == K) {
            return -1;
        }
        order[row] = (POSITION_T)pivot_col;
        used[pivot_col] = 1;

        const FQ_ELEM inv = fq_inv(work[row][pivot_col]);
        for (uint32_t r = row + 1; r < K; r++) {
            const FQ_ELEM scale = fq_mul(work[r][pivot_col], inv);
            if (scale == 0) {
                continue;
            }
            for (uint32_t col = 0; col < K; col++) {
                work[r][col] = fq_sub(work[r][col], fq_mul(scale, work[row][col]));
            }
        }
    }
    return 0;
}

static int pivot_positions(POSITION_T pivots_out[K],
                           POSITION_T nonpivots_out[K],
                           const uint8_t pivots[N_pad]) {
    uint32_t pivots_count = 0;
    uint32_t nonpivots_count = 0;
    for (uint32_t col = 0; col < N; col++) {
        if (pivots[col]) {
            if (pivots_count >= K) {
                return -1;
            }
            pivots_out[pivots_count++] = (POSITION_T)col;
        } else {
            if (nonpivots_count >= K) {
                return -1;
            }
            nonpivots_out[nonpivots_count++] = (POSITION_T)col;
        }
    }
    return pivots_count == K && nonpivots_count == K ? 0 : -1;
}

static int init_or_update_inverse(FQ_ELEM inv[K][K],
                                  const FQ_ELEM a[K][K],
                                  uint32_t row);
static int root_bit_for_row(uint8_t *root_bit,
                            const FQ_ELEM a[K][K],
                            const FQ_ELEM inv[K][K],
                            uint32_t row);

int celeres_tag_compress_from_rref(uint8_t out[CELERES_TAG_BYTES],
                                   const uint8_t rref[RREF_MAT_PACKEDBYTES]) {
    generator_mat_t full = {0};
    uint8_t pivots[N_pad] = {0};
    POSITION_T pivot_cols[K];
    POSITION_T nonpivot_cols[K];
    POSITION_T order[K];
    FQ_ELEM block[K][K];
    FQ_ELEM ordered[K][K];

    memset(out, 0, CELERES_TAG_BYTES);
    expand_to_rref(&full, rref, pivots);
    if (pivot_positions(pivot_cols, nonpivot_cols, pivots) != 0) {
        return -1;
    }

    for (uint32_t byte = 0; byte < N8; byte++) {
        for (uint32_t bit = 0; bit < 8 && (8u * byte + bit) < N; bit++) {
            out[byte] |= (uint8_t)(pivots[8u * byte + bit] << bit);
        }
    }

    for (uint32_t row = 0; row < K; row++) {
        for (uint32_t col = 0; col < K; col++) {
            block[row][col] = full.values[row][nonpivot_cols[col]];
        }
    }

    if (choose_column_order(order, block) != 0) {
        return -1;
    }

    size_t bit_pos = 0;
    for (uint32_t i = 0; i < K; i++) {
        bit_write(out + N8, bit_pos, order[i], CELERES_TAG_ORDER_BITS_PER_ENTRY);
        bit_pos += CELERES_TAG_ORDER_BITS_PER_ENTRY;
    }

    for (uint32_t col = 0; col < K; col++) {
        for (uint32_t row = 0; row < K; row++) {
            ordered[row][col] = block[row][order[col]];
        }
    }

    FQ_ELEM inv[K][K] = {{0}};
    memset(out + N8 + CELERES_TAG_ORDER_BYTES, 0,
           CELERES_TAG_ROOT_BITS_BYTES);
    for (uint32_t row = 0; row < K; row++) {
        uint8_t root_bit = 0;
        if (root_bit_for_row(&root_bit, ordered, inv, row) != 0) {
            return -1;
        }
        out[N8 + CELERES_TAG_ORDER_BYTES + (row >> 3)] =
            (uint8_t)(out[N8 + CELERES_TAG_ORDER_BYTES + (row >> 3)] |
                      (uint8_t)(root_bit << (row & 7u)));
        if (init_or_update_inverse(inv, ordered, row) != 0) {
            return -1;
        }
    }

    bit_pos = 0;
    for (uint32_t row = 0; row < K; row++) {
        for (uint32_t col = row + 1u; col < K; col++) {
            bit_write(out + N8 + CELERES_TAG_ORDER_BYTES +
                          CELERES_TAG_ROOT_BITS_BYTES,
                      bit_pos, ordered[row][col], BITS_TO_REPRESENT(Q));
            bit_pos += BITS_TO_REPRESENT(Q);
        }
    }
    return 0;
}

static int parse_column_order(POSITION_T order[K], const uint8_t *in) {
    uint8_t seen[K] = {0};
    size_t bit_pos = 0;
    for (uint32_t i = 0; i < K; i++) {
        const uint32_t value = bit_read(in, bit_pos, CELERES_TAG_ORDER_BITS_PER_ENTRY);
        bit_pos += CELERES_TAG_ORDER_BITS_PER_ENTRY;
        if (value >= K || seen[value]) {
            return -1;
        }
        seen[value] = 1;
        order[i] = (POSITION_T)value;
    }
    return 0;
}

static int update_inverse(FQ_ELEM inv[K][K],
                          const FQ_ELEM a[K][K],
                          uint32_t size) {
    FQ_ELEM inv_u[K] = {0};
    FQ_ELEM v_inv[K] = {0};
    FQ_ELEM schur = a[size][size];

    for (uint32_t row = 0; row < size; row++) {
        for (uint32_t col = 0; col < size; col++) {
            inv_u[row] = fq_add(inv_u[row], fq_mul(inv[row][col], a[col][size]));
            v_inv[col] = fq_add(v_inv[col], fq_mul(a[size][row], inv[row][col]));
        }
    }
    for (uint32_t i = 0; i < size; i++) {
        schur = fq_sub(schur, fq_mul(a[size][i], inv_u[i]));
    }
    if (schur == 0) {
        return -1;
    }

    const FQ_ELEM inv_schur = fq_inv(schur);
    for (uint32_t row = 0; row < size; row++) {
        for (uint32_t col = 0; col < size; col++) {
            inv[row][col] = fq_add(inv[row][col],
                                   fq_mul(fq_mul(inv_u[row], v_inv[col]), inv_schur));
        }
    }
    for (uint32_t row = 0; row < size; row++) {
        inv[row][size] = fq_sub(0, fq_mul(inv_u[row], inv_schur));
        inv[size][row] = fq_sub(0, fq_mul(v_inv[row], inv_schur));
    }
    inv[size][size] = inv_schur;
    return 0;
}

static int init_or_update_inverse(FQ_ELEM inv[K][K],
                                  const FQ_ELEM a[K][K],
                                  uint32_t row) {
    if (row == 0) {
        if (a[0][0] == 0) {
            return -1;
        }
        inv[0][0] = fq_inv(a[0][0]);
        return 0;
    }
    if (row + 1u < K) {
        return update_inverse(inv, a, row);
    }
    return 0;
}

static void row_affine_prefix(FQ_ELEM base[K],
                              FQ_ELEM slope[K],
                              FQ_ELEM *tail_square,
                              const FQ_ELEM a[K][K],
                              const FQ_ELEM inv[K][K],
                              uint32_t row) {
    memset(base, 0, K * sizeof(FQ_ELEM));
    memset(slope, 0, K * sizeof(FQ_ELEM));
    *tail_square = 0;

    for (uint32_t col = row + 1u; col < K; col++) {
        *tail_square = fq_add(*tail_square,
                              fq_mul(a[row][col], a[row][col]));
    }

    if (row == 0) {
        return;
    }

    FQ_ELEM rhs[K] = {0};
    FQ_ELEM diagonal_column[K] = {0};
    for (uint32_t prev = 0; prev < row; prev++) {
        diagonal_column[prev] = a[prev][row];
        for (uint32_t col = row + 1u; col < K; col++) {
            rhs[prev] = fq_sub(rhs[prev],
                               fq_mul(a[row][col], a[prev][col]));
        }
    }

    for (uint32_t col = 0; col < row; col++) {
        for (uint32_t prev = 0; prev < row; prev++) {
            base[col] = fq_add(base[col], fq_mul(inv[col][prev], rhs[prev]));
            slope[col] = fq_sub(slope[col],
                                fq_mul(inv[col][prev],
                                       diagonal_column[prev]));
        }
    }
}

static FQ_ELEM row_norm_for_diagonal_candidate(const FQ_ELEM base[K],
                                               const FQ_ELEM slope[K],
                                               FQ_ELEM tail_square,
                                               uint32_t row,
                                               FQ_ELEM t) {
    FQ_ELEM norm = fq_add(tail_square, fq_mul(t, t));
    for (uint32_t col = 0; col < row; col++) {
        const FQ_ELEM value = fq_add(base[col], fq_mul(slope[col], t));
        norm = fq_add(norm, fq_mul(value, value));
    }
    return norm;
}

static int row_diagonal_candidates(FQ_ELEM candidates[2],
                                   const FQ_ELEM a[K][K],
                                   const FQ_ELEM inv[K][K],
                                   uint32_t row) {
    FQ_ELEM base[K] = {0};
    FQ_ELEM slope[K] = {0};
    FQ_ELEM tail_square = 0;
    uint32_t count = 0;

    row_affine_prefix(base, slope, &tail_square, a, inv, row);

    for (uint32_t candidate = 0; candidate < Q; candidate++) {
        const FQ_ELEM t = (FQ_ELEM)candidate;
        if (row_norm_for_diagonal_candidate(base, slope, tail_square,
                                            row, t) == Q - 1u) {
            if (count < 2) {
                candidates[count] = t;
            }
            count++;
        }
    }
    if (count == 0 || count > 2) {
        return -1;
    }
    if (count == 1) {
        candidates[1] = candidates[0];
    }
    if (candidates[1] < candidates[0]) {
        const FQ_ELEM tmp = candidates[0];
        candidates[0] = candidates[1];
        candidates[1] = tmp;
    }
    return 0;
}

static int solve_row_from_self_orthogonality(FQ_ELEM a[K][K],
                                             const FQ_ELEM inv[K][K],
                                             uint32_t row,
                                             uint8_t root_bit) {
    FQ_ELEM base[K] = {0};
    FQ_ELEM slope[K] = {0};
    FQ_ELEM tail_square = 0;
    FQ_ELEM candidates[2];

    if (row_diagonal_candidates(candidates, a, inv, row) != 0) {
        return -1;
    }
    row_affine_prefix(base, slope, &tail_square, a, inv, row);
    (void)tail_square;

    const FQ_ELEM t = candidates[root_bit != 0 ? 1 : 0];
    a[row][row] = t;
    for (uint32_t col = 0; col < row; col++) {
        a[row][col] = fq_add(base[col], fq_mul(slope[col], t));
    }
    return 0;
}

static int root_bit_for_row(uint8_t *root_bit,
                            const FQ_ELEM a[K][K],
                            const FQ_ELEM inv[K][K],
                            uint32_t row) {
    FQ_ELEM candidates[2];
    if (row_diagonal_candidates(candidates, a, inv, row) != 0) {
        return -1;
    }
    if (a[row][row] == candidates[0]) {
        *root_bit = 0;
        return 0;
    }
    if (a[row][row] == candidates[1]) {
        *root_bit = 1;
        return 0;
    }
    return -1;
}

static int reconstruct_ordered_block_diagonal(FQ_ELEM a[K][K],
                                              const uint8_t *root_bits,
                                              const uint8_t *values) {
    FQ_ELEM inv[K][K] = {{0}};
    size_t bit_pos = 0;

    memset(a, 0, K * K * sizeof(FQ_ELEM));
    for (uint32_t row = 0; row < K; row++) {
        for (uint32_t col = row + 1u; col < K; col++) {
            const uint32_t value = bit_read(values, bit_pos,
                                            BITS_TO_REPRESENT(Q));
            bit_pos += BITS_TO_REPRESENT(Q);
            if (value >= Q) {
                return -1;
            }
            a[row][col] = (FQ_ELEM)value;
        }

        const uint8_t root_bit =
            (uint8_t)((root_bits[row >> 3] >> (row & 7u)) & 1u);
        if (solve_row_from_self_orthogonality(a, inv, row, root_bit) != 0) {
            return -1;
        }

        for (uint32_t prev = 0; prev < row; prev++) {
            if (scalar_product(a[row], a[prev], K) != 0) {
                return -1;
            }
        }
        if (inner_square(a[row], K) != Q - 1u) {
            return -1;
        }
        if (init_or_update_inverse(inv, a, row) != 0) {
            return -1;
        }
    }
    return 0;
}

int celeres_tag_expand_to_rref(uint8_t out[RREF_MAT_PACKEDBYTES],
                               const uint8_t tag[CELERES_TAG_BYTES]) {
    uint8_t pivots[N_pad] = {0};
    POSITION_T pivot_cols[K];
    POSITION_T nonpivot_cols[K];
    POSITION_T order[K];
    FQ_ELEM ordered[K][K];
    FQ_ELEM block[K][K];
    generator_mat_t full = {0};

    for (uint32_t byte = 0; byte < N8; byte++) {
        for (uint32_t bit = 0; bit < 8 && (8u * byte + bit) < N; bit++) {
            pivots[8u * byte + bit] = (uint8_t)((tag[byte] >> bit) & 1u);
        }
    }
    if (pivot_positions(pivot_cols, nonpivot_cols, pivots) != 0 ||
        parse_column_order(order, tag + N8) != 0) {
        return -1;
    }
    if (reconstruct_ordered_block_diagonal(
            ordered,
            tag + N8 + CELERES_TAG_ORDER_BYTES,
            tag + N8 + CELERES_TAG_ORDER_BYTES +
                CELERES_TAG_ROOT_BITS_BYTES) != 0) {
        return -1;
    }

    memset(block, 0, sizeof(block));
    for (uint32_t ordered_col = 0; ordered_col < K; ordered_col++) {
        const POSITION_T original_col = order[ordered_col];
        for (uint32_t row = 0; row < K; row++) {
            block[row][original_col] = ordered[row][ordered_col];
        }
    }

    for (uint32_t row = 0; row < K; row++) {
        full.values[row][pivot_cols[row]] = 1;
        for (uint32_t col = 0; col < K; col++) {
            full.values[row][nonpivot_cols[col]] = block[row][col];
        }
    }
    compress_rref(out, &full, pivots);
    return 0;
}

#endif
