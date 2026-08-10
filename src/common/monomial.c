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

#include "monomial_mat.h"
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include <string.h>

#define POS_BITS BITS_TO_REPRESENT(N-1)

#define POS_MASK (((POSITION_T) 1 << POS_BITS) - 1)
#define MONOMIAL_SALT_INPUT_BYTES (SEED_LENGTH_BYTES + HASH_DIGEST_LENGTH + sizeof(uint16_t))

void yt_shuffle_state(SHAKE_STATE_STRUCT *shake_monomial_state, POSITION_T permutation[N]) {
   uint64_t rand_u64;
   POSITION_T tmp;
   POSITION_T x;
   int c;

   csprng_randombytes((unsigned char *) &rand_u64,
                             sizeof(rand_u64),
                             shake_monomial_state);
   c = 0;

   for (int i = 0; i < N; i++) {
      do {
         if (c == (64/POS_BITS)-1) {
            csprng_randombytes((unsigned char *) &rand_u64,
                                sizeof(rand_u64),
                                shake_monomial_state);
            c = 0;
         }
         x = rand_u64 & (POS_MASK);
         rand_u64 = rand_u64 >> POS_BITS;
         c = c + 1;
      } while (x >= N);

      tmp = permutation[i];
      permutation[i] = permutation[x];
      permutation[x] = tmp;
   } 
}
/* expands a monomial matrix, given a PRNG seed and a salt (used for ephemeral
 * monomial matrices */
void monomial_sample_salt(monomial_t *res,
                          const unsigned char seed[SEED_LENGTH_BYTES],
                          const unsigned char salt[HASH_DIGEST_LENGTH],
                          const uint16_t round_index) {
    SHAKE_STATE_STRUCT shake_monomial_state = {0};
    uint8_t shake_input_buffer[MONOMIAL_SALT_INPUT_BYTES];
    memcpy(shake_input_buffer, seed, SEED_LENGTH_BYTES);
    memcpy(shake_input_buffer + SEED_LENGTH_BYTES, salt, HASH_DIGEST_LENGTH);
    memcpy(shake_input_buffer + SEED_LENGTH_BYTES + HASH_DIGEST_LENGTH, &round_index, sizeof(uint16_t));

    initialize_csprng(&shake_monomial_state, shake_input_buffer,
                      MONOMIAL_SALT_INPUT_BYTES);
    fq_star_rnd_state_elements(&shake_monomial_state, res->coefficients, N);
    for (uint32_t i = 0; i < N; i++) {
        res->permutation[i] = i;
    }

    /* FY shuffle on the permutation */
    yt_shuffle_state(&shake_monomial_state, res->permutation);
} /* end monomial_mat_seed_expand */

void monomial_sample_prikey(monomial_t *res,
                            const unsigned char seed[PRIVATE_KEY_SEED_LENGTH_BYTES]) {
    SHAKE_STATE_STRUCT shake_monomial_state = {0};
    initialize_csprng(&shake_monomial_state, seed, PRIVATE_KEY_SEED_LENGTH_BYTES);
    fq_star_rnd_state_elements(&shake_monomial_state, res->coefficients, N);
    for (uint32_t i = 0; i < N; i++) {
        res->permutation[i] = i;
    }
    /* FY shuffle on the permutation */
    yt_shuffle_state(&shake_monomial_state, res->permutation);
} /* end monomial_mat_seed_expand */

/* composes a compactly stored action of a monomial on an IS with a regular
 * monomial.
 * NOTE: Only the permutation is computed, as this is the only thing we need
 * since the adaption of canonical forms.
 */
void monomial_compose_action(monomial_action_IS_t *out,
                             const monomial_t *Q_in,
                             const monomial_action_IS_t *in) {
    /* to compose with monomial_action_IS_t, reverse the convention
     * for Q storage: store in permutation[i] the idx of the source column landing
     * as the i-th after the GQ product, and in coefficients[i] the coefficient
     * by which the column is multiplied upon landing */
    monomial_t reverse_Q;
    for (uint32_t i = 0; i < N; i++) {
        reverse_Q.permutation[Q_in->permutation[i]] = i;
    }
    /* compose actions out = Q_in*in */
    for (uint32_t i = 0; i < K; i++) {
        out->permutation[i] = reverse_Q.permutation[in->permutation[i]];
    }
}

void CosetRep(uint8_t *b,
              const monomial_action_IS_t *Q_star) {
    memset(b, 0, N8);
    for (uint32_t i = 0; i < K; i++) {
        const uint32_t limb = (Q_star->permutation[i])/8;
        const uint32_t pos  = (Q_star->permutation[i])%8;
        b[limb] ^= 1u << pos;
    }
}

int CheckCanonicalAction(const uint8_t* const b) {
    uint32_t w = 0;
    for (uint32_t i = 0; i < N8; i++) {
        w += (uint32_t)__builtin_popcount(b[i]);
    }

    return w == K;
}
