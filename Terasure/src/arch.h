//
// Created by syd on 17-4-26.
//

#ifndef LIBMSR_ARCH_H
#define LIBMSR_ARCH_H

#include "gf.h"

#define AVX2

#ifdef AVX2

#include <immintrin.h>

typedef __m256i encode_t;


#define REGION_SIZE 512
#define REGION_BLOCKS (REGION_SIZE/sizeof(encode_t))

static encode_t low_table[256];
static encode_t high_table[256];
static encode_t mask_lo;
static encode_t mask_hi;

static void init_arch() {
    int low_array[32];
    int high_array[32];
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 16; j++)
            low_array[j] = gf_mul(i, j);
        for (int j = 0; j < 16; j++)
            low_array[j + 16] = gf_mul(i, j);
        for (int j = 0; j < 16; j++)
            high_array[j] = gf_mul(i, j << 4);
        for (int j = 0; j < 16; j++)
            high_array[j + 16] = gf_mul(i, j << 4);
        low_table[i] = _mm256_setr_epi8(low_array[0], low_array[1], low_array[2], low_array[3], low_array[4],
                                        low_array[5], low_array[6], low_array[7], low_array[8], low_array[9],
                                        low_array[10], low_array[11], low_array[12], low_array[13], low_array[14],
                                        low_array[15], low_array[16], low_array[17], low_array[18], low_array[19],
                                        low_array[20], low_array[21], low_array[22], low_array[23], low_array[24],
                                        low_array[25], low_array[26], low_array[27], low_array[28], low_array[29],
                                        low_array[30], low_array[31]);
        high_table[i] = _mm256_setr_epi8(high_array[0], high_array[1], high_array[2], high_array[3], high_array[4],
                                         high_array[5], high_array[6], high_array[7], high_array[8], high_array[9],
                                         high_array[10], high_array[11], high_array[12], high_array[13], high_array[14],
                                         high_array[15], high_array[16], high_array[17], high_array[18], high_array[19],
                                         high_array[20], high_array[21], high_array[22], high_array[23], high_array[24],
                                         high_array[25], high_array[26], high_array[27], high_array[28], high_array[29],
                                         high_array[30], high_array[31]);
    }
    mask_lo = _mm256_set1_epi8(0x0f);
    mask_hi = _mm256_set1_epi8(0xf0);
}

__attribute__((always_inline))
static inline encode_t xor_region(encode_t input1, encode_t input2) {
    return _mm256_xor_si256(input1, input2);
}

__attribute__((always_inline))
static inline encode_t multiply_region(encode_t input, uint8_t x) {

    __m256i low = _mm256_and_si256(input, mask_lo);
    __m256i high = _mm256_and_si256(input, mask_hi);

    __m256i low_t = low_table[x], high_t = high_table[x];
    high = _mm256_srli_epi16(high, 4);

    __m256i right = _mm256_shuffle_epi8(low_t, low);
    __m256i left = _mm256_shuffle_epi8(high_t, high);


    return _mm256_xor_si256(left, right);
}

__attribute__((always_inline))
static inline encode_t zero() {
    return _mm256_setzero_si256();
}

__attribute__((always_inline))
static void inline prefetch(const encode_t *src) {
    _mm_prefetch(src, _MM_HINT_NTA);
}

__attribute__((always_inline))
static void inline store(encode_t *dst, encode_t val) {
    _mm256_stream_si256(dst, val);
}

#endif


#endif //LIBMSR_ARCH_H
