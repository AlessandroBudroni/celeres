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

#include "utils.h"
#include <string.h>
#include <stdlib.h>


void cswap(uintptr_t *a,
           uintptr_t *b,
           const uintptr_t mask) {
    *a ^= (mask & *b);
    *b ^= (mask & *a);
    *a ^= (mask & *b);
}

#ifdef USE_AVX2
#include <immintrin.h>

int verify(const uint8_t *a,
           const uint8_t *b,
           size_t len) {
    size_t i;
    uint64_t r;
    __m256i f, g, h;

    h = _mm256_setzero_si256();
    for(i=0;i<len/32;i++) {
        f = _mm256_loadu_si256((__m256i *)&a[32*i]);
        g = _mm256_loadu_si256((__m256i *)&b[32*i]);
        f = _mm256_xor_si256(f,g);
        h = _mm256_or_si256(h,f);
    }
    r = 1u - (uint32_t)_mm256_testz_si256(h,h);

    a += 32*i;
    b += 32*i;
    len -= 32*i;

    for(i=0;i<len;i++) {
        r |= a[i] ^ b[i];
    }

    r = (-r) >> 63;
    return r;
}
#else

int verify(const uint8_t *a,
           const uint8_t *b,
           const size_t len) {
    uint8_t r = 0;

    for(size_t i=0;i<len;i++) {
        r |= a[i] ^ b[i];
    }

    return (-(uint64_t)r) >> 63;
}

#endif
