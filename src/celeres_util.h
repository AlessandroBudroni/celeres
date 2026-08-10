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

#include "parameters.h"

int celeres_randombytes(uint8_t *out, size_t len);

int celeres_log_round_up_u32(uint32_t n);

void celeres_derive_bits(uint8_t out[SEED_LENGTH_BYTES],
                         const uint8_t seed[SEED_LENGTH_BYTES],
                         const uint8_t salt[HASH_DIGEST_LENGTH],
                         uint16_t round,
                         uint32_t ring_index);

void celeres_dummy_commitment(uint8_t out[HASH_DIGEST_LENGTH],
                              const uint8_t salt[HASH_DIGEST_LENGTH],
                              uint16_t round,
                              uint32_t dummy_index);

void celeres_ihmt_build_root_and_path(const uint8_t *commitments_in,
                                      int logn,
                                      int32_t index,
                                      uint8_t root[HASH_DIGEST_LENGTH],
                                      uint8_t *path);

void celeres_ihmt_reconstruct_root(const uint8_t leaf[HASH_DIGEST_LENGTH],
                                   const uint8_t *path,
                                   int logn,
                                   uint8_t root[HASH_DIGEST_LENGTH]);

void celeres_sample_challenge(uint8_t challenge[T],
                              const uint8_t digest[HASH_DIGEST_LENGTH]);
