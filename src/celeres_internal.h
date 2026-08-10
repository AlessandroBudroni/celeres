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

#pragma once

#include "codes.h"

int celeres_rref(generator_mat_t *g,
                      uint8_t is_pivot_column[N_pad]);

int celeres_rref_with_hint(generator_mat_t *g,
                                uint8_t is_pivot_column[N_pad],
                                const uint8_t pivot_hint[N_pad]);

void celeres_diagonal_scaling_canonical_form(generator_mat_t *g);

void celeres_diagonal_scaling_canonical_form_ct(generator_mat_t *g);

int celeres_cf(generator_mat_t *g);

int celeres_cf_sign_ct(generator_mat_t *g);

int celeres_cf_with_hint(generator_mat_t *g,
                              const uint8_t pivot_hint[N_pad]);

int celeres_sf_from_rref_monomial(
    generator_mat_t *out,
    uint8_t out_pivots[N_pad],
    const uint8_t rref[RREF_MAT_PACKEDBYTES],
    const monomial_t *q);

int celeres_cf_from_rref_monomial_is(
    generator_mat_t *cf,
    const uint8_t rref[RREF_MAT_PACKEDBYTES],
    const monomial_t *q,
    int use_ct_scaling);

void celeres_perm_rep(POSITION_T out[N],
                           const monomial_t *m);
