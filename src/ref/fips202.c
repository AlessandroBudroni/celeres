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

#include <stdint.h>
#include <string.h>

#include "fips202.h"
#include "keccakf1600.h"

#define NROUNDS 24
#define ROL(a, offset) (((a) << (offset)) ^ ((a) >> (64 - (offset))))



static void keccak_absorb(uint64_t *s,
                          uint32_t r,
                          const uint8_t *m, size_t mlen,
                          uint8_t p)
{
   while (mlen >= r) {
      KeccakF1600_StateXORBytes(s, m, 0, r);
      KeccakF1600_StatePermute(s);
      mlen -= r;
      m += r;
   }

   if(mlen > 0) {
      KeccakF1600_StateXORBytes(s, m, 0, mlen);
   }

   if(mlen == r-1) {
      p |= 128;
      KeccakF1600_StateXORBytes(s, &p, mlen, 1);
   } else {
      KeccakF1600_StateXORBytes(s, &p, mlen, 1);
      p = 128;
      KeccakF1600_StateXORBytes(s, &p, r-1, 1);
   }
}



static void keccak_squeezeblocks(uint8_t *h, size_t nblocks,
                                 uint64_t *s, uint32_t r)
{
   while(nblocks > 0) {
      KeccakF1600_StatePermute(s);
      KeccakF1600_StateExtractBytes(s, h, 0, r);
      h += r;
      nblocks--;
   }
}


static void keccak_inc_init(uint64_t *s_inc)
{
   size_t i;

   for (i = 0; i < 25; ++i) {
      s_inc[i] = 0;
   }
   s_inc[25] = 0;
}

static void keccak_inc_absorb(uint64_t *s_inc, uint32_t r, const uint8_t *m,
                              size_t mlen)
{
   /* Recall that s_inc[25] is the non-absorbed bytes xored into the state */
   while (mlen + s_inc[25] >= r) {

      KeccakF1600_StateXORBytes(s_inc, m, s_inc[25], r-s_inc[25]);
      mlen -= (size_t)(r - s_inc[25]);
      m += r - s_inc[25];
      s_inc[25] = 0;

      KeccakF1600_StatePermute(s_inc);
   }

   KeccakF1600_StateXORBytes(s_inc, m, s_inc[25], mlen);
   s_inc[25] += mlen;
}


static void keccak_inc_finalize(uint64_t *s_inc, uint32_t r, uint8_t p)
{
   /* After keccak_inc_absorb, we are guaranteed that s_inc[25] < r,
      so we can always use one more byte for p in the current state. */
   if(s_inc[25] == r-1) {
      p |= 128;
      KeccakF1600_StateXORBytes(s_inc, &p, s_inc[25], 1);
   } else {
      KeccakF1600_StateXORBytes(s_inc, &p, s_inc[25], 1);
      p = 128;
      KeccakF1600_StateXORBytes(s_inc, &p, r-1, 1);
   }
   s_inc[25] = 0;
}


static void keccak_inc_squeeze(uint8_t *h, size_t outlen,
                               uint64_t *s_inc, uint32_t r)
{
   size_t len;
   if(outlen < s_inc[25]) {
      len = outlen;
   } else {
      len = s_inc[25];
   }

   KeccakF1600_StateExtractBytes(s_inc, h, r-s_inc[25], len);
   h += len;
   outlen -= len;
   s_inc[25] -= len;

   /* Then squeeze the remaining necessary blocks */
   while (outlen > 0) {
      KeccakF1600_StatePermute(s_inc);

      if(outlen < r) {
         len = outlen;
      } else {
         len = r;
      }
      KeccakF1600_StateExtractBytes(s_inc, h, 0, len);
      h += len;
      outlen -= len;
      s_inc[25] = r - len;
   }
}

void shake128_inc_init(shake128incctx *state)
{
   keccak_inc_init(state->ctx);
}

void shake128_inc_absorb(shake128incctx *state, const uint8_t *input,
                         size_t inlen)
{
   keccak_inc_absorb(state->ctx, SHAKE128_RATE, input, inlen);
}

void shake128_inc_finalize(shake128incctx *state)
{
   keccak_inc_finalize(state->ctx, SHAKE128_RATE, 0x1F);
}

void shake128_inc_squeeze(uint8_t *output, size_t outlen,
                          shake128incctx *state)
{
   keccak_inc_squeeze(output, outlen, state->ctx, SHAKE128_RATE);
}

void shake256_inc_init(shake256incctx *state)
{
   keccak_inc_init(state->ctx);
}

void shake256_inc_absorb(shake256incctx *state, const uint8_t *input,
                         size_t inlen)
{
   keccak_inc_absorb(state->ctx, SHAKE256_RATE, input, inlen);
}

void shake256_inc_finalize(shake256incctx *state)
{
   keccak_inc_finalize(state->ctx, SHAKE256_RATE, 0x1F);
}

void shake256_inc_squeeze(uint8_t *output, size_t outlen,
                          shake256incctx *state)
{
   keccak_inc_squeeze(output, outlen, state->ctx, SHAKE256_RATE);
}


void shake128_absorb(shake128ctx *state, const uint8_t *input, size_t inlen)
{
   int i;
   for (i = 0; i < 25; i++)
      state->ctx[i] = 0;
   keccak_absorb(state->ctx, SHAKE128_RATE, input, inlen, 0x1F);
}


void shake128_squeezeblocks(uint8_t *output, size_t nblocks,
                            shake128ctx *state)
{
   keccak_squeezeblocks(output, nblocks, state->ctx, SHAKE128_RATE);
}

void shake128(uint8_t *output, size_t outlen, const uint8_t *input,
              size_t inlen)
{
   shake128incctx state;
   keccak_inc_init(state.ctx);

   /* Absorb input */
   keccak_inc_absorb(state.ctx, SHAKE128_RATE, input, inlen);
   keccak_inc_finalize(state.ctx, SHAKE128_RATE, 0x1F);

   /* Squeeze output */
   keccak_inc_squeeze(output, outlen, state.ctx, SHAKE128_RATE);
}

void shake128_ctx_clone(shake128ctx *dest, const shake128ctx *src)
{
   memcpy(dest, src, sizeof(shake128ctx));
}

void shake256_absorb(shake256ctx *state, const uint8_t *input, size_t inlen)
{
   int i;
   for (i = 0; i < 25; i++)
      state->ctx[i] = 0;

   keccak_absorb(state->ctx, SHAKE256_RATE, input, inlen, 0x1F);
}


void shake256_squeezeblocks(uint8_t *output, size_t nblocks,
                            shake256ctx *state)
{
   keccak_squeezeblocks(output, nblocks, state->ctx, SHAKE256_RATE);
}


void shake256(uint8_t *output, size_t outlen,
              const uint8_t *input, size_t inlen)
{
   shake256incctx state;

   keccak_inc_init(state.ctx);

   /* Absorb input */
   keccak_inc_absorb(state.ctx, SHAKE256_RATE, input, inlen);
   keccak_inc_finalize(state.ctx, SHAKE256_RATE, 0x1F);

   /* Squeeze output */
   keccak_inc_squeeze(output, outlen, state.ctx, SHAKE256_RATE);
}

void shake256_ctx_clone(shake256ctx *dest, const shake256ctx *src)
{
   memcpy(dest, src, sizeof(shake256ctx));
}



void sha3_256(uint8_t *output, const uint8_t *input, size_t inlen)
{
   sha3_256incctx state;
   keccak_inc_init(state.ctx);

   /* Absorb input */
   keccak_inc_absorb(state.ctx, SHA3_256_RATE, input, inlen);
   keccak_inc_finalize(state.ctx, SHA3_256_RATE, 0x06);

   /* Squeeze output */
   keccak_inc_squeeze(output, 32, state.ctx, SHA3_256_RATE);
}
void sha3_256_inc_init(sha3_256incctx *state)
{
   keccak_inc_init(state->ctx);
}

void sha3_256_inc_absorb(sha3_256incctx *state, const uint8_t *input,
                         size_t inlen)
{
   keccak_inc_absorb(state->ctx, SHA3_256_RATE, input, inlen);
}

void sha3_256_inc_finalize(uint8_t *output, sha3_256incctx *state)
{
   uint8_t t[SHA3_256_RATE];
   keccak_inc_finalize(state->ctx, SHA3_256_RATE, 0x06);

   keccak_squeezeblocks(t, 1, state->ctx, SHA3_256_RATE);

   for (size_t i = 0; i < 32; i++) {
      output[i] = t[i];
   }
}

void sha3_384_inc_init(sha3_384incctx *state)
{
   keccak_inc_init(state->ctx);
}

void sha3_384_inc_absorb(sha3_384incctx *state, const uint8_t *input,
                         size_t inlen)
{
   keccak_inc_absorb(state->ctx, SHA3_384_RATE, input, inlen);
}

void sha3_384_inc_finalize(uint8_t *output, sha3_384incctx *state)
{
   uint8_t t[SHA3_384_RATE];
   keccak_inc_finalize(state->ctx, SHA3_384_RATE, 0x06);

   keccak_squeezeblocks(t, 1, state->ctx, SHA3_384_RATE);

   for (size_t i = 0; i < 48; i++) {
      output[i] = t[i];
   }
}


void sha3_384(uint8_t *output, const uint8_t *input, size_t inlen)
{
   sha3_384incctx state;
   keccak_inc_init(state.ctx);

   /* Absorb input */
   keccak_inc_absorb(state.ctx, SHA3_384_RATE, input, inlen);
   keccak_inc_finalize(state.ctx, SHA3_384_RATE, 0x06);

   /* Squeeze output */
   keccak_inc_squeeze(output, 48, state.ctx, SHA3_384_RATE);
}


void sha3_512(uint8_t *output, const uint8_t *input, size_t inlen)
{
   sha3_512incctx state;
   keccak_inc_init(state.ctx);

   /* Absorb input */
   keccak_inc_absorb(state.ctx, SHA3_512_RATE, input, inlen);
   keccak_inc_finalize(state.ctx, SHA3_512_RATE, 0x06);

   /* Squeeze output */
   keccak_inc_squeeze(output, 64, state.ctx, SHA3_512_RATE);
}

void sha3_512_inc_init(sha3_512incctx *state)
{
   keccak_inc_init(state->ctx);
}

void sha3_512_inc_absorb(sha3_512incctx *state, const uint8_t *input,
                         size_t inlen)
{
   keccak_inc_absorb(state->ctx, SHA3_512_RATE, input, inlen);
}

void sha3_512_inc_finalize(uint8_t *output, sha3_512incctx *state)
{
   uint8_t t[SHA3_512_RATE];
   keccak_inc_finalize(state->ctx, SHA3_512_RATE, 0x06);

   keccak_squeezeblocks(t, 1, state->ctx, SHA3_512_RATE);

   for (size_t i = 0; i < 64; i++) {
      output[i] = t[i];
   }
}
