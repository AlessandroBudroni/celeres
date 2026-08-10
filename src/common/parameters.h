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
#include <stdint.h>

/* Seed tree max size is computed according to Parameter Generation Script in Utilities folder */

/***************************** Common Parameters ******************************/
#define Q (127)
#define Qm1 (Q-1)
#define FQ_ELEM uint8_t
#define FQ_DOUBLEPREC uint16_t
#define POSITION_T uint16_t


/***************************** Toy parameters *********************************/
#if !defined(CATEGORY)
#define N (16)
#define K (8)
#define SEED_LENGTH_BYTES (16)
#define T (106)
#define W (57)

/*************************** Custom Category 1 *******************************/
#elif CATEGORY == 272
#define N (272)
#define K (136)

#define SEED_LENGTH_BYTES (16)

#ifndef TARGET
#define TARGET 192
#elif TARGET != 192
#error CATEGORY 272 only supports TARGET 192 in CELERES
#endif

#define T (192)
#define W (36)
#define TREE_OFFSETS {0, 0, 0, 0, 0, 0, 0, 0, 128}
#define TREE_NODES_PER_LEVEL {1, 2, 4, 8, 16, 32, 64, 128, 128}
#define TREE_LEAVES_PER_LEVEL {0, 0, 0, 0, 0, 0, 0, 64, 128}
#define TREE_SUBROOTS 2
#define TREE_LEAVES_START_INDICES {255, 191}
#define TREE_CONSECUTIVE_LEAVES {128, 64}
#define MAX_PUBLISHED_SEEDS 87

/*************************** Custom Category 3 *******************************/
#elif CATEGORY == 416
#define N (416)
#define K (208)
#define SEED_LENGTH_BYTES (24)

#ifndef TARGET
#define TARGET 220
#elif TARGET != 220
#error CATEGORY 416 only supports TARGET 220 in CELERES
#endif

#define T (220)
#define W (68)

#define TREE_OFFSETS {0, 0, 0, 0, 0, 0, 0, 8, 56}
#define TREE_NODES_PER_LEVEL {1, 2, 4, 8, 16, 32, 64, 120, 192}
#define TREE_LEAVES_PER_LEVEL {0, 0, 0, 0, 0, 0, 4, 24, 192}
#define TREE_SUBROOTS 3
#define TREE_LEAVES_START_INDICES {247, 223, 123}
#define TREE_CONSECUTIVE_LEAVES {192, 24, 4}
#define MAX_PUBLISHED_SEEDS 119

/*************************** Custom Category 5 *******************************/
#elif CATEGORY == 564
#define N (564)
#define K (282)
#define SEED_LENGTH_BYTES (32)

#ifndef TARGET
#define TARGET 345
#elif TARGET != 345
#error CATEGORY 564 only supports TARGET 345 in CELERES
#endif

#define T (345)
#define W (75)

#define TREE_OFFSETS {0, 0, 0, 0, 0, 2, 2, 2, 50, 178}
#define TREE_NODES_PER_LEVEL {1, 2, 4, 8, 16, 30, 60, 120, 192, 256}
#define TREE_LEAVES_PER_LEVEL {0, 0, 0, 0, 1, 0, 0, 24, 64, 256}
#define TREE_SUBROOTS 4
#define TREE_LEAVES_START_INDICES {433, 369, 217, 30}
#define TREE_CONSECUTIVE_LEAVES {256, 64, 24, 1}
#define MAX_PUBLISHED_SEEDS 169

#else
#error define category for parameters
#endif

/* number of bytes needed to store K or N bits */
#define K8 ((K+7u)/8u)
#define N8 ((N+7u)/8u)

#define NEXT_MULTIPLE(x,n) ((((x)+((n)-1u))/(n))*(n))

#if defined(USE_AVX2) || defined(USE_NEON)
#define N_K_pad NEXT_MULTIPLE(N-K, 32)
#define N_pad   NEXT_MULTIPLE(N, 32)
#define K_pad   NEXT_MULTIPLE(K, 32)
#else
#define N_K_pad (N-K)
#define N_pad   N
#define K_pad   K
#endif

#define Q_pad   NEXT_MULTIPLE(Q, 8)

/***************** Derived parameters *****************************************/
/*length of the output of the cryptographic hash, in bytes */
#define HASH_DIGEST_LENGTH (2*SEED_LENGTH_BYTES)
#define SALT_LENGTH_BYTES HASH_DIGEST_LENGTH


/* length of the private key seed doubled to avoid multikey attacks */
#define PRIVATE_KEY_SEED_LENGTH_BYTES (2*SEED_LENGTH_BYTES)

#define MASK_Q ((1 << BITS_TO_REPRESENT(Q)) - 1)
#define MASK_N ((1 << BITS_TO_REPRESENT(N)) - 1)


#define IS_REPRESENTABLE_IN_D_BITS(D, N)                \
  (((unsigned long) N>=(1UL << (D-1)) && (unsigned long) N<(1UL << D)) ? D : -1)

#define BITS_TO_REPRESENT(N)                            \
  (N == 0 ? 1 : (15                                     \
                 + IS_REPRESENTABLE_IN_D_BITS( 1, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 2, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 3, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 4, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 5, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 6, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 7, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 8, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 9, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(10, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(11, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(12, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(13, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(14, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(15, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(16, N)    \
                 )                                      \
   )

#define LOG2(L) ( (BITS_TO_REPRESENT(L) > BITS_TO_REPRESENT(L-1)) ? (BITS_TO_REPRESENT(L-1)) : (BITS_TO_REPRESENT(L)) )

#define NUM_LEAVES_SEED_TREE (T)
#define NUM_NODES_SEED_TREE ((2*NUM_LEAVES_SEED_TREE) - 1)

#define RREF_MAT_PACKEDBYTES ((BITS_TO_REPRESENT(Q)*(N-K)*K + 7)/8 + (N + 7)/8)

#define SEED_TREE_MAX_PUBLISHED_BYTES (MAX_PUBLISHED_SEEDS*SEED_LENGTH_BYTES + 1)

#ifdef USE_AVX2
#if defined(CATEGORY_5)
#define LESS_USE_CUSTOM_HISTOGRAM
#endif
#endif
