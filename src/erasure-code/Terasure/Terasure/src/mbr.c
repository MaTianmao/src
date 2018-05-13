//
// Created by msn on 17-10-31.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <sys/time.h>
#include "mbr.h"
#include "arch.h"

#define MAX_PARITY 16
void mbr_init(mbr_conf *conf, int n, int k, int d, void *(allocator)(size_t), void (*deallocator)(void *)){


	gf_init();
	init_arch();

	conf->n = n;
	conf->k = k;
	conf->d = d;

	conf->allocator = allocator;
	conf->deallocator = deallocator;

    conf->matrix = conf->allocator((conf->n - conf->k) * sizeof(uint8_t *));
    for(int i = 0; i < conf->n - conf->k; i++){
        conf->matrix[i] = conf->allocator(conf->d * sizeof(uint8_t));
        for(int j = 0; j < conf->d; j++)
            conf->matrix[i][j] = gf_div((uint8_t)1, i ^ (j + conf->n - conf->k));
    }
}

__attribute__((always_inline))
static inline void super_fast_cal(int len, const int K, const int R, 
    uint8_t **data, uint8_t **output, uint8_t **Matrix){

	encode_t **data_ptr = (encode_t **)data;
    encode_t **output_ptr = (encode_t **)output;

    encode_t res[R*K][REGION_BLOCKS];
    for(int i = 0; i < R; i++)
        for(int j = 0; j < K; j++)
            for(int k = 0; k < REGION_BLOCKS; k++)
                res[i*K+j][k] = zero();

    for(int i = 0; i < len / sizeof(encode_t); i += REGION_BLOCKS){ 
    	for(int j = 0; j < K; j++){
            for(int k = 0; k < K; k++){ 
                for(int x = 0; x < REGION_BLOCKS; x++){ 
                    encode_t tmp = data_ptr[j*K+k][i + x]; 
                    for(int r = 0; r < R; r++) 
                        res[r*K+k][x] = xor_region(res[r*K+k][x], multiply_region(tmp, Matrix[r][j]));
                } 
            } 
        }
        for(int r = 0; r < R; r++){ 
            for(int k = 0; k < K; k++){
                for(int x = 0; x < REGION_BLOCKS; x++){ 
                    store(&output_ptr[r*K+k][i + x], res[r*K+k][x]);
                    res[r*K+k][x] = zero(); 
                } 
            }   
        }
    }    
}

void mbr_encode(int len, mbr_conf *conf, uint8_t **data, 
    uint8_t **output, uint8_t **memory_pre_allocated){

    switch ((conf->n - conf->k)) {
        case 2:
            switch (conf->d) {
                case 7:
                    super_fast_cal(len, 7, 2, data, memory_pre_allocated, conf->matrix);
                    break;
                case 9:
                    super_fast_cal(len, 9, 2, data, memory_pre_allocated, conf->matrix);
                    break;
                case 11:
                    super_fast_cal(len, 11, 2, data, memory_pre_allocated, conf->matrix);
                    break;
                case 13:
                    super_fast_cal(len, 13, 2, data, memory_pre_allocated, conf->matrix);
                    break;
            }
            break;
        case 3:
            switch (conf->d) {
                case 7:
                    super_fast_cal(len, 7, 3, data, memory_pre_allocated, conf->matrix);
                    break;
                case 9:
                    super_fast_cal(len, 9, 3, data, memory_pre_allocated, conf->matrix);
                    break;
                case 11:
                    super_fast_cal(len, 11, 3, data, memory_pre_allocated, conf->matrix);
                    break;
                case 13:
                    super_fast_cal(len, 13, 3, data, memory_pre_allocated, conf->matrix);
                    break;
            }
            break;
        case 4:
            switch (conf->d) {
                case 7:
                    super_fast_cal(len, 7, 4, data, memory_pre_allocated, conf->matrix);
                    break;
                case 9:
                    super_fast_cal(len, 9, 4, data, memory_pre_allocated, conf->matrix);
                    break;
                case 11:
                    super_fast_cal(len, 11, 4, data, memory_pre_allocated, conf->matrix);
                    break;
                case 13:
                    super_fast_cal(len, 13, 4, data, memory_pre_allocated, conf->matrix);
                    break;
            }
            break;
        case 5:
            switch (conf->d) {
                case 7:
                    super_fast_cal(len, 7, 5, data, memory_pre_allocated, conf->matrix);
                    break;
                case 9:
                    super_fast_cal(len, 9, 5, data, memory_pre_allocated, conf->matrix);
                    break;
                case 11:
                    super_fast_cal(len, 11, 5, data, memory_pre_allocated, conf->matrix);
                    break;
                case 13:
                    super_fast_cal(len, 13, 5, data, memory_pre_allocated, conf->matrix);
                    break;
            }
            break;
        case 6:
            switch (conf->d) {
                case 7:
                    super_fast_cal(len, 7, 6, data, memory_pre_allocated, conf->matrix);
                    break;
                case 9:
                    super_fast_cal(len, 9, 6, data, memory_pre_allocated, conf->matrix);
                    break;
                case 11:
                    super_fast_cal(len, 11, 6, data, memory_pre_allocated, conf->matrix);
                    break;
                case 13:
                    super_fast_cal(len, 13, 6, data, memory_pre_allocated, conf->matrix);
                    break;
            }
            break;
    }

    for(int i = 0; i < conf->n; i++){
        for(int j = 0; j < conf->d; j++){
            if(i < conf->k) 
                output[i*conf->d+j] = data[i*conf->d+j];
            else
                output[i*conf->d+j] = memory_pre_allocated[(i - conf->k)*conf->d+j];
        }
    }
}

bool mbr_check_broken(int x, int *broken_id, int broken_num){
	for(int i = 0; i < broken_num; i++){
		if(broken_id[i] == x) 
            return true;
	}
	return false;
}

void mbr_invert_matrix(uint8_t **in_mat, uint8_t **out_mat, int n){
	int i, j, k;
	for(i = 0; i < n; i++){
		for(j = i; j < n; j++){
			if(in_mat[j][i]) 
                break;
		}

		if(j != i){
			for(k = 0; k < n; k++){
				uint8_t temp = in_mat[i][k];
				in_mat[i][k] = in_mat[j][k];
				in_mat[j][k] = temp;
				temp = out_mat[i][k];
				out_mat[i][k] = out_mat[j][k];
				out_mat[j][k] = temp;
			}
		}

		uint8_t temp = in_mat[i][i];
		for(j = 0; j < n; j++){
			in_mat[i][j] = gf_div(in_mat[i][j], temp);
			out_mat[i][j] = gf_div(out_mat[i][j], temp);
		}

		for(j = 0; j < n; j++){
			if(j != i){
				temp = in_mat[j][i];
				for(k = 0; k < n; k++){
					in_mat[j][k] ^= gf_mul(temp, in_mat[i][k]);
					out_mat[j][k] ^= gf_mul(temp, out_mat[i][k]);
				}
			}
		}
	}
}

void mbr_recover_context_init(mbr_conf *conf, 
    mbr_recover_context *context, int *broken_id, int broken_num){

	context->broken_id = broken_id;
	context->broken_num = broken_num;
	context->recover_id = conf->allocator(conf->k * sizeof(int));
	uint8_t **matrix;
	matrix = conf->allocator(conf->k * sizeof(uint8_t *));
    context->matrix_half = conf->allocator(conf->k * sizeof(uint8_t *));
	int re_num = 0;

	for(int i = 0; i < conf->n; i++){
		if(re_num == conf->k) break;
		if(!mbr_check_broken(i, broken_id, broken_num)){
			matrix[re_num] = conf->allocator(conf->k * sizeof(uint8_t));
            context->matrix_half[re_num] = conf->allocator((conf->d - conf->k) * sizeof(uint8_t));
            if(i >= conf->k){
                for(int j = 0; j < conf->d; j++){
                    if(j < conf->k){
                        matrix[re_num][j] = conf->matrix[i - conf->k][j];
                    }
                    else 
                        context->matrix_half[re_num][j - conf->k] = conf->matrix[i-conf->k][j];
                }
            }
            else{
                for(int j = 0; j < conf->d; j++){
                    if(j < conf->k) 
                        matrix[re_num][j] = (i == j) ? 1 : 0;
                    else 
                        context->matrix_half[re_num][j-conf->k] = 0;
                }
            }
            context->recover_id[re_num] = i;
            re_num++;
		}
	}

	context->inv_matrix = conf->allocator(conf->k * sizeof(uint8_t *));
	for(int i = 0; i < conf->k; i++){
		context->inv_matrix[i] = conf->allocator(conf->k * sizeof(uint8_t));
		for(int j = 0; j < conf->k; j++){
			context->inv_matrix[i][j] = (i==j) ? 1 : 0;
		}
	}

	mbr_invert_matrix(matrix, context->inv_matrix, conf->k);
	for(int i = 0; i < conf->k; i++) 
        conf->deallocator(matrix[i]);
	conf->deallocator(matrix);
}

void mbr_recover(int len, mbr_conf *conf, mbr_recover_context *context, 
     uint8_t **data, uint8_t **output){
	
    encode_t **data_ptr = (encode_t **)data;
    encode_t **output_ptr = (encode_t **)output;
    const int K = conf->k;
    const int D = conf->d - conf->k;
    encode_t res1[K*D][REGION_BLOCKS];
    encode_t res2[K*K][REGION_BLOCKS];
    encode_t res3[K*K][REGION_BLOCKS];

    for(int i = 0; i < K; i++)
        for(int j = 0; j < D; j++)
            for(int x = 0; x < REGION_BLOCKS; x++)
                res1[i*D+j][x] = zero();
    for(int i = 0; i < len/sizeof(encode_t); i+=REGION_BLOCKS){
        for(int j = 0; j < K; j++){
            for(int k = 0; k < D; k++){
                for(int x = 0; x < REGION_BLOCKS; x++){
                    encode_t tmp = data_ptr[context->recover_id[j]*conf->d+k + K][i + x];
                    for(int r = 0; r < K; r++){
                        res1[r*D+k][x] = xor_region(res1[r*D+k][x], 
                            multiply_region(tmp, context->inv_matrix[r][j]));
                    }
                }
            }
        }
        for(int r = 0; r < K; r++){
            for(int j = 0; j < D; j++){
                for(int x = 0; x < REGION_BLOCKS; x++){
                    store(&output_ptr[r*conf->d+j + K][i + x], res1[r*D+j][x]);
                    store(&output_ptr[(j + K)*conf->d+r][i + x], res1[r*D+j][x]);
                    res1[r*D+j][x] = zero();
                }
            }
        }
    }

    for(int i = 0; i < K; i++)
        for(int j = 0; j < K; j++)
            for(int x = 0; x < REGION_BLOCKS; x++){
                res2[i*K+j][x] = zero();
                res3[i*K+j][x] = zero();
            }
    for(int i = 0; i < len/sizeof(encode_t); i += REGION_BLOCKS){
        for(int j = 0; j < D; j++){
            for(int k = 0; k < K; k++){
                for(int x = 0; x < REGION_BLOCKS; x++){
                    encode_t tmp = output_ptr[(j + K)*conf->d+k][i + x];
                    for(int r = 0; r < K; r++){
                        res2[r*K+k][x] = xor_region(res2[r*K+k][x], 
                            multiply_region(tmp, context->matrix_half[r][j]));
                    }
                }
            }
        }
        for(int j = 0; j < K; j++){
            for(int k = 0; k < K; k++){
                for(int x = 0; x < REGION_BLOCKS; x++){
                    encode_t tmp1 = data_ptr[context->recover_id[j]*conf->d+k][i + x];
                    encode_t tmp2 = res2[j*K+k][x];
                    encode_t tmp = xor_region(tmp1, tmp2);
                    for(int r = 0; r < K; r++){
                        res3[r*K+k][x] = xor_region(res3[r*K+k][x], 
                            multiply_region(tmp, context->inv_matrix[r][j]));
                    }
                }
            }
        }
        for(int j = 0; j < K; j++){
            for(int k = 0; k < K; k++){
                for(int x = 0; x < REGION_BLOCKS; x++){
                    store(&output_ptr[j*conf->d+k][i + x], res3[j*K+k][x]);
                    res2[j*K+k][x] = zero();
                    res3[j*K+k][x] = zero();
                }
            }
        }
    }
}

void mbr_free_recover_context(mbr_conf *conf, mbr_recover_context *context){
    for(int i = 0; i < conf->k; i++){
        conf->deallocator(context->inv_matrix[i]);
        conf->deallocator(context->matrix_half[i]);
    }
    conf->deallocator(context->inv_matrix);
    conf->deallocator(context->matrix_half);
    conf->deallocator(context->recover_id);
}

void mbr_regenerate_context_init(mbr_conf *conf, mbr_regenerate_context *context,int broken_id){
    context->broken_id = broken_id;
    uint8_t **matrix = conf->allocator(conf->d * sizeof(uint8_t *));
    int xcount = 0;
    context->broken_vector = conf->allocator(conf->d * sizeof(uint8_t));
    if(broken_id < conf->k){
        for(int i = 0; i < conf->d; i++)
            context->broken_vector[i] = (i == broken_id) ? 1 :0;
    }
    else{
        memcpy(context->broken_vector, conf->matrix[broken_id - conf->k], conf->d * sizeof(uint8_t));
    }

    for(int i = 0; i < conf->n; i++){
        if(i == broken_id) continue;
        matrix[xcount] = conf->allocator(conf->d * sizeof(uint8_t));
        if(i < conf->k){
            for(int j = 0; j < conf->d; j++) 
                matrix[xcount][j] = (i == j) ? 1 : 0;
        }
        else{
            memcpy(matrix[xcount], conf->matrix[i - conf->k], conf->d * sizeof(uint8_t));
        }
        xcount++;
    }

    context->inv_matrix = conf->allocator(conf->d * sizeof(uint8_t *));
    for(int i = 0; i < conf->d; i++){
        context->inv_matrix[i] = conf->allocator(conf->d * sizeof(uint8_t));
        for(int j = 0; j < conf->d; j++)
            context->inv_matrix[i][j] = (i == j) ? 1 : 0;
    }
    mbr_invert_matrix(matrix, context->inv_matrix, conf->d);
    for(int i = 0; i < conf->d; i++) 
        conf->deallocator(matrix[i]);
    conf->deallocator(matrix);
}

void mbr_regenerate_one_node(int len, mbr_conf *conf, mbr_regenerate_context *context, uint8_t **data, int id, uint8_t **output){
    encode_t **output_ptr = (encode_t **)output;
    encode_t **data_ptr = (encode_t **)data;
    const int D = conf->d;
    encode_t res[REGION_BLOCKS];

    for(int x = 0; x < REGION_BLOCKS; x++)
        res[x] = zero();
    for(int i = 0; i < len/sizeof(encode_t); i+=REGION_BLOCKS){
        for(int j = 0; j < D; j++){
            for(int x = 0; x<REGION_BLOCKS; x++){
                res[x] = xor_region(res[x], multiply_region(data_ptr[id*D+j][i+x], context->broken_vector[j]));
            }
        }
        int bk_id = (id > context->broken_id) ? id - 1 : id;
        for(int x = 0; x < REGION_BLOCKS; x++){
            store(&output_ptr[bk_id][i + x], res[x]);
            res[x] = zero();
        }
    }
}
void mbr_regenerate(int len, mbr_conf *conf, mbr_regenerate_context *context, uint8_t **data, uint8_t **output){
    for(int i = 0; i < conf->n; i++){
        if(i == context->broken_id) continue;
        mbr_regenerate_one_node(len, conf, context, data, i, output);
    }
    encode_t **output_ptr = (encode_t **)output;
    const int D = conf->d;
    encode_t res[D][REGION_BLOCKS];
    for(int i = 0; i < D; i++){
        for(int x = 0; x < REGION_BLOCKS; x++)
             res[i][x] = zero();
    }
    for(int i = 0; i < len / sizeof(encode_t); i+=REGION_BLOCKS){
        for(int j = 0; j < D; j++){
            for(int x = 0; x < REGION_BLOCKS; x++){
                encode_t tmp = output_ptr[j][i + x];
                for(int r = 0; r < D; r++){
                    res[r][x] = xor_region(res[r][x], multiply_region(tmp, context->inv_matrix[r][j]));
                }
            }
        }
        for(int j = 0; j < D; j++){
            for(int x = 0; x < REGION_BLOCKS; x++){
                store(&output_ptr[j][i + x], res[j][x]);
                res[j][x] = zero();
            }
        }
    }
}
void mbr_free_regenerate_context(mbr_conf *conf, mbr_regenerate_context *context){
    for(int i = 0; i < conf->d; i++){
        conf->deallocator(context->inv_matrix[i]);
    }
	conf->deallocator(context->inv_matrix);
	conf->deallocator(context->broken_vector);
}

void mbr_free_conf(mbr_conf *conf){
    for(int i = 0; i < conf->n - conf->k; i++){
        conf->deallocator(conf->matrix[i]);
    }
	conf->deallocator(conf->matrix);
}
