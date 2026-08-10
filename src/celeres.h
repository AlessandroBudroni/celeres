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

#include <stddef.h>
#include <stdint.h>

#include "celeres_tag.h"
#include "parameters.h"
#include "permutation_compress.h"

#define CELERES_PUBLIC_MATRIX_SEED_BYTES SEED_LENGTH_BYTES
#define CELERES_SECRET_SEED_BYTES PRIVATE_KEY_SEED_LENGTH_BYTES
#define CELERES_PUBLIC_KEY_BYTES RREF_MAT_PACKEDBYTES
#define CELERES_PERM_BYTES CELERES_PERM_COMPRESSED_BYTES

#ifndef CELERES_RING_SIZE
#error "CELERES_RING_SIZE must be defined by the build system"
#endif

typedef struct {
    uint8_t g_seed[CELERES_PUBLIC_MATRIX_SEED_BYTES];
    uint8_t t_seed[CELERES_PUBLIC_MATRIX_SEED_BYTES];
} celeres_public_parameters_t;

typedef struct {
    uint8_t q_seed[CELERES_SECRET_SEED_BYTES];
    uint8_t q_tilde_seed[CELERES_SECRET_SEED_BYTES];
} celeres_secret_key_t;

typedef struct {
    uint8_t rref[CELERES_PUBLIC_KEY_BYTES];
} celeres_public_key_t;

typedef struct {
    uint8_t tag[CELERES_TAG_BYTES];
} celeres_tag_t;

size_t celeres_signature_bytes(void);

int celeres_generate_parameters(celeres_public_parameters_t *params);

int celeres_keygen(const celeres_public_parameters_t *params,
                        celeres_public_key_t *pk,
                        celeres_tag_t *tag,
                        celeres_secret_key_t *sk);

int celeres_sign(const celeres_public_parameters_t *params,
                      const celeres_public_key_t *ring,
                      const celeres_tag_t *tags,
                      uint32_t signer_index,
                      const celeres_secret_key_t *sk,
                      const uint8_t *message,
                      size_t message_len,
                      uint8_t *signature,
                      size_t *signature_len);

int celeres_verify(const celeres_public_parameters_t *params,
                        const celeres_public_key_t *ring,
                        const uint8_t *message,
                        size_t message_len,
                        const uint8_t *signature,
                        size_t signature_len);

int celeres_link(const uint8_t *signature,
                      size_t signature_len,
                      const uint8_t *other_signature,
                      size_t other_signature_len);
