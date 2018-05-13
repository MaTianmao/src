#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <sys/time.h>
#include "hh.h"
#include "arch.h"

#define rep(i,n) for(int i=0;i<n;++i)
#define MAX_PARITY 10

__attribute__((always_inline))
static inline void super_fast_encode(int len, const int K, const int R, uint8_t **data, uint8_t **output, uint8_t **Cauchy){
    encode_t **data_ptr = (encode_t **)data;
    encode_t **output_ptr = (encode_t **)output;

    encode_t resA[MAX_PARITY];
    encode_t resB[MAX_PARITY];
    encode_t res[MAX_PARITY];
    for(int i=0;i<MAX_PARITY;i++)
        res[i] = resA[i] = resB[i] = zero();

    int mid = len / sizeof(encode_t)/2;
    for(int i=0; i < mid; ++i){
        int ii = i + mid;
        for(int k = 0; k < K;++k){
            encode_t dataA = data_ptr[k][i], dataB = data_ptr[k][ii];
            for(int r=0;r<R;++r){
               resA[r] = xor_region(resA[r], multiply_region(dataA, Cauchy[r][k]));
               resB[r] = xor_region(resB[r], multiply_region(dataB, Cauchy[r][k]));
            }
        }
        for(int r=1;r<R;r++) {
            encode_t dA = zero();
            for(int j = (K-1)*(r-1)/(R-1); j<(K-1)*r/(R-1);++j)
                dA = xor_region(dA, multiply_region(data_ptr[j][i], Cauchy[1][j]));
            resB[r] = xor_region(resB[r], dA);
        }
        resA[1] = xor_region(resA[1], resB[1]);
        for(int r=0;r<R;r++){
            store(&output_ptr[r][i], resA[r]);
            store(&output_ptr[r][ii], resB[r]);
            resA[r] = resB[r] = zero();
        }
    }
}
__attribute__((always_inline))
static inline void super_fast_regenerate(int len, int K, int R, uint8_t **data, uint8_t *output, int broken, uint8_t **Cauchy) {
    encode_t **data_ptr = (encode_t **)data;
    encode_t *output_ptr = (encode_t *)output;
    int num = len / sizeof(encode_t);
    encode_t A, B;
    A = B = zero();
    encode_t parity[R];
    if (R<=2){
        if (broken<K){
            for(int i=0;i<num;++i){
                A = data_ptr[K][i];
                for(int j=0;j<K;++j)
                    if (j!=broken)
                        A = xor_region(A, multiply_region(data_ptr[j][i], Cauchy[0][j]));
                A = multiply_region(A, gf_div(1, Cauchy[0][broken]));
                store(&output_ptr[i], A);
                //A = zero();
            }
        } else {
            for(int i=0;i<num;++i){
                A = zero();
                for(int j=0;j<K;++j)
                    A = xor_region(A, multiply_region(data_ptr[j][i], Cauchy[broken-K][j]));
                store(&output_ptr[i], A);
            }
        }
    } else {
        int mid = num / 2;
        if (broken < K-1) {
            int block = broken * (R-1) / (K-1);
            while ((block+1)*(K-1)/(R-1)<=broken) ++block;
            int block_st = block*(K-1)/(R-1), block_ed = (block+1)*(K-1)/(R-1);
            for(int i=0;i<mid;++i) {
                int ii = i + mid;
                B = data_ptr[K][i], A=data_ptr[K+block+1][i];
                for(int j=0;j<block_st;++j) {
                    B = xor_region(B, multiply_region(data_ptr[j][i], Cauchy[0][j]));
                    A = xor_region(A, multiply_region(data_ptr[j][i], Cauchy[block+1][j]));
                }
                for(int j=block_st;j<block_ed;++j)
                    if (j!=broken) {
                        B = xor_region(B, multiply_region(data_ptr[j][ii], Cauchy[0][j]));
                        A = xor_region(A, multiply_region(data_ptr[j][ii], Cauchy[block + 1][j]));
                        A = xor_region(A, multiply_region(data_ptr[j][i], Cauchy[1][j]));
                    }
                for(int j=block_ed;j<K;++j) {
                    B = xor_region(B, multiply_region(data_ptr[j][i], Cauchy[0][j]));
                    A = xor_region(A, multiply_region(data_ptr[j][i], Cauchy[block+1][j]));
                }
                B = multiply_region(B, gf_div(1, Cauchy[0][broken]));
                A = xor_region(A, multiply_region(B, Cauchy[block+1][broken]));
                A = multiply_region(A, gf_div(1, Cauchy[1][broken]));
                store(&output_ptr[i], A);
                store(&output_ptr[ii], B);
            }
        } else if (broken == K-1) {
            for(int i=0;i<mid;++i) {
                int ii = i + mid;
                B = data_ptr[K][i];
                for(int j=0;j<K;++j)
                    if (j!=broken)
                        B = xor_region(B, multiply_region(data_ptr[j][i], Cauchy[0][j]));
                B = multiply_region(B, gf_div(1, Cauchy[0][broken]));
                for(int r = 1; r<R; ++r) parity[r] = multiply_region(B, Cauchy[r][K-1]);
                for(int k=0;k<K-1;++k)
                {
                    encode_t tmp = data_ptr[k][i];
                    for(int r = 1; r<R;++r)
                        parity[r] = xor_region(parity[r], multiply_region(tmp, Cauchy[r][k]));//parity[r] ^= Cauchy[r][k] * tmp;
                }
                A = xor_region(data_ptr[K+1][i], parity[1]);
                for(int r = 2;r<R;++r){
                    A = xor_region(A, xor_region(data_ptr[K+r][i], parity[r]));//A ^= data_ptr[K+r][ii] ^ parity[r];
                }
                A = multiply_region(A, gf_div(1, Cauchy[1][broken]));
                store(&output_ptr[i], A);
                store(&output_ptr[ii], B);
            }
        } else {
            for(int i=0; i < mid; ++i){
                int ii = i + mid;
                A = B = zero();
                for(int k = 0; k < K;++k){
                    A = xor_region(A, multiply_region(data_ptr[k][i], Cauchy[broken-K][k]));
                    B = xor_region(B, multiply_region(data_ptr[k][ii], Cauchy[broken-K][k]));
                }
                if (broken > K)
                {
                    for(int j = (K-1)*(broken-K-1)/(R-1); j<(K-1)*(broken-K)/(R-1);++j)
                        B = xor_region(B, multiply_region(data_ptr[j][i], Cauchy[1][j]));
                    if (broken == K+1)
                        A = xor_region(A, B);
                }
                store(&output_ptr[i], A);
                store(&output_ptr[ii], B);
            }
        }
    }
}


void hh_encode(int len, hh_conf *conf, uint8_t **data, uint8_t **output){
    switch (conf->k) {
        case 8:
            switch (conf->r) {
                case 3:
                    super_fast_encode(len, 8, 3, data, output, conf->matrix);
                    break;
                case 4:
                    super_fast_encode(len, 8, 4, data, output, conf->matrix);
                    break;
            }
            break;
        case 10:
            switch (conf->r) {
                case 3:
                    super_fast_encode(len, 10, 3, data, output, conf->matrix);
                    break;
                case 4:
                    super_fast_encode(len, 10, 4, data, output, conf->matrix);
                    break;
            }
            break;
        case 12:
            switch (conf->r) {
                case 3:
                    super_fast_encode(len, 12, 3, data, output, conf->matrix);
                    break;
                case 4:
                    super_fast_encode(len, 12, 4, data, output, conf->matrix);
                    break;
            }
            break;
    }
    rep(i, conf->r) data[i+conf->k] = output[i];
//    if(R==1){
//        encode_batch(1)
//    }else if(R==2){
//        for(int i=0;i<len/sizeof(encode_t);i++){
//            for(int k=0;k<K;k++){
//                encode_t tmp = data_ptr[k][i];
//                for(int r=0;r<R;r++){
//                    res[r] = xor_region(res[r],multiply_region(tmp,Cauchy[r][k]));
//                }
//            }
//            for(int r=0;r<R;r++){
//                store(&output_ptr[r][i],res[r]);
//                res[r] = zero();
//            }
//        }
//    }
}
/*
void hh_regenerate(int len, int N, int K, uint8_t **data, uint8_t *output, int broken) {
    int R = N - K;
    encode_t **data_ptr = (encode_t **)data;
    encode_t *output_ptr = (encode_t *)output;
    int num = len / sizeof(encode_t);
    encode_t A, B;
    A = B = zero();
    encode_t parity[R];
    if (R<=2){
        if (broken<K){
            for(int i=0;i<num;++i){
                A = data_ptr[K][i];
                for(int j=0;j<K;++j)
                    if (j!=broken)
                        A = xor_region(A, multiply_region(data_ptr[j][i], Cauchy[0][j]));
                A = multiply_region(A, gf_div(1, Cauchy[0][broken]));
                store(&output_ptr[i], A);
                //A = zero();
            }
        } else {
            for(int i=0;i<num;++i){
                A = zero();
                for(int j=0;j<K;++j)
                    A = xor_region(A, multiply_region(data_ptr[j][i], Cauchy[broken-K][j]));
                store(&output_ptr[i], A);
            }
        }
    } else {
        int mid = num / 2;
        if (broken < K-1) {
            int block = broken * (R-1) / (K-1);
            while ((block+1)*(K-1)/(R-1)<=broken) ++block;
            for(int i=0;i<mid;++i) {
                int ii = i + mid;
                B = data_ptr[K][ii];
                for(int j=0;j<K;++j)
                    if (j!=broken)
                        B = xor_region(B, multiply_region(data_ptr[j][ii], Cauchy[0][j]));
                B = multiply_region(B, gf_div(1, Cauchy[0][broken])); // B = b[broken]
                //(K-1)x / (R-1) <= broken < (K-1)(x+1) / (R-1)
                A = multiply_region(B, Cauchy[block+1][broken]);
                for(int k=0;k<K;++k)
                    if (k!=broken)
                        A = xor_region(A, multiply_region(data_ptr[k][ii], Cauchy[block+1][k])); // A = f(block+1)b
                A = xor_region(A, data_ptr[K+block+1][ii]); // A = f(block+1)b ^ b[K + block+1 ] = f1(...a...)
                for(int k=block*(K-1)/(R-1);k<(block+1)*(K-1)/(R-1);++k)
                    if(k!=broken)
                        A = xor_region(A, multiply_region(data_ptr[k][i], Cauchy[1][k]));
                A = multiply_region(A, gf_div(1, Cauchy[1][broken]));

                store(&output_ptr[i], A);
                store(&output_ptr[ii], B);
            }
        } else if (broken == K-1) {
            for(int i=0;i<mid;++i) {
                int ii = i + mid;
                B = data_ptr[K][ii];
                for(int j=0;j<K;++j)
                    if (j!=broken)
                        B = xor_region(B, multiply_region(data_ptr[j][ii], Cauchy[0][j]));
                B = multiply_region(B, gf_div(1, Cauchy[0][broken]));
                for(int r = 1; r<R; ++r) parity[r] = multiply_region(B, Cauchy[r][K-1]);
                for(int k=0;k<K-1;++k)
                {
                    encode_t tmp = data_ptr[k][ii];
                    for(int r = 1; r<R;++r)
                        parity[r] = xor_region(parity[r], multiply_region(tmp, Cauchy[r][k]));//parity[r] ^= Cauchy[r][k] * tmp;
                }
                A = xor_region(data_ptr[K+1][i], parity[1]);
                for(int r = 2;r<R;++r){
                    A = xor_region(A, xor_region(data_ptr[K+r][ii], parity[r]));//A ^= data_ptr[K+r][ii] ^ parity[r];
                }
                A = multiply_region(A, gf_div(1, Cauchy[1][broken]));
                store(&output_ptr[i], A);
                store(&output_ptr[ii], B);
            }
        } else {
            for(int i=0; i < mid; ++i){
                int ii = i + mid;
                A = B = zero();
                for(int k = 0; k < K;++k){
                    A = xor_region(A, multiply_region(data_ptr[k][i], Cauchy[broken-K][k]));
                    B = xor_region(B, multiply_region(data_ptr[k][ii], Cauchy[broken-K][k]));
                }
                if (broken > K)
                {
                    for(int j = (K-1)*(broken-K-1)/(R-1); j<(K-1)*(broken-K)/(R-1);++j)
                        B = xor_region(B, multiply_region(data_ptr[j][i], Cauchy[1][j]));
                    if (broken == K+1)
                        A = xor_region(A, B);
                }
                store(&output_ptr[i], A);
                store(&output_ptr[ii], B);
            }
        }
    }
}
*/

void hh_regenerate(int len, hh_regenerate_context *context, hh_conf *conf, uint8_t **data, uint8_t *output) {
    int broken = context->broken;
    switch (conf->k) {
        case 8:
            switch (conf->r) {
                case 3:
                    super_fast_regenerate(len, 8, 3, data, output, broken, conf->matrix);
                    break;
                case 4:
                    super_fast_regenerate(len, 8, 4, data, output, broken, conf->matrix);
                    break;
            }
            break;
        case 10:
            switch (conf->r) {
                case 3:
                    super_fast_regenerate(len, 10, 3, data, output, broken, conf->matrix);
                    break;
                case 4:
                    super_fast_regenerate(len, 10, 4, data, output, broken, conf->matrix);
                    break;
            }
            break;
        case 12:
            switch (conf->r) {
                case 3:
                    super_fast_regenerate(len, 12, 3, data, output, broken, conf->matrix);
                    break;
                case 4:
                    super_fast_regenerate(len, 12, 4, data, output, broken, conf->matrix);
                    break;
            }
            break;
    }
}

void hh_init(hh_conf *conf, int n, int k, void *(*allocator)(size_t), void (*deallocator)(void *)) {
    int r = n - k;

    if (r <= 1 || r > k || n >= GF_SIZE - 1) {
        return;
    }

    gf_init();
    init_arch();

    conf->n = n;
    conf->k = k;
    conf->r = r;

    //Can be improved through padding with a number with is not a power of 2

    conf->allocate = allocator;
    conf->deallocate = deallocator;

    conf->matrix = conf->allocate(conf->r * sizeof(uint8_t*));
    rep(i, conf->r) {
        conf->matrix[i] = conf->allocate(conf->k * sizeof(uint8_t));
        rep(j, conf->k)
            conf->matrix[i][j] = gf_div(1,i ^ (j+conf->r));
    }
}

void hh_get_regenerate_offset(int len, hh_regenerate_context *context, hh_conf *conf, int broken) {
    int k = conf->k, r = conf->r;
    context->broken = broken;
    context->size = conf->allocate(sizeof(int) * conf->n);
    context->offset = conf->allocate(sizeof(int) * conf->n);
    memset(context->size, 0, sizeof(int) * conf->n);
    memset(context->offset, 0, sizeof(int) * conf->n);
    if (r <= 2) {
        rep(i, k)if (i != broken) {
                context->size[i] = len;
                context->offset[i] = 0;
            }
        if (broken < k) {
            context->size[k] = len;
            context->offset[k] = 0;
        }
    } else if (broken < k - 1) {
        int block = broken * (r - 1) / (k - 1);
        while ((block + 1) * (k - 1) / (r - 1) <= broken) block++;
        rep(i, k + 1)if (i != broken) {
                context->size[i] = len / 2;
                context->offset[i] = len / 2;
            }
        context->size[k + block + 1] = len / 2;
        context->offset[k + block + 1] = len / 2;
        for (int i = block * (k - 1) / (r - 1); i < (block + 1) * (k - 1) / (r - 1); i++)
            if (i != broken) {
                context->size[i] += len / 2;
                context->offset[i] = 0;
            }
    } else if (broken == k - 1) {
        rep(i, k + 1)
            if (i != broken) {
                context->offset[i] = len / 2;
                context->size[i] = len / 2;
            }
        for (int i = k + 2; i < conf->n; i++) {
            context->offset[i] = len / 2;
            context->size[i] = len / 2;
        }
        context->offset[k + 1] = 0;
        context->size[k + 1] = len / 2;
    } else {
        rep(i, k) {
            context->offset[i] = 0;
            context->size[i] = len;
        }
    }
}

void hh_free_regenerate_context(hh_conf *conf, hh_regenerate_context *context) {
    conf->deallocate(context->offset);
    conf->deallocate(context->size);
}

void hh_free_conf(hh_conf *conf) {
    rep(i, conf->r) {
        conf->deallocate(conf->matrix[i]);
    }
    conf->deallocate(conf->matrix);
}
