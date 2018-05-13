//
// Created by msn on 17-10-14.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <sys/time.h>
#include "reed_solomon.h"
#include "arch.h"
#define MAX_PARITY 16
void rs_init(rs_conf *conf, int n, int k, void *(allocator)(size_t), void (*deallocator)(void *)){
	int r = n - k;

	gf_init();

	init_arch();

	conf->n = n;
	conf->k = k;
	conf->r = r;

	conf->allocator = allocator;
	conf->deallocator = deallocator;
}

void vandermonde_matrix(rs_conf *conf){
	conf->matrix = conf->allocator(conf->n * sizeof(uint8_t*));
    conf->matrix_type = 1;

	for(uint8_t i = 0; i < conf->n; i++){
		conf->matrix[i] = conf->allocator(conf->k * sizeof(uint8_t));
		for(uint8_t j = 0; j < conf->k; j++){
            if(i == 0) 
                conf->matrix[i][j] = 1;
			else 
                conf->matrix[i][j] = gf_mul(conf->matrix[i-1][j], j + 1);
		}
	}
}

void cauchy_matrix(rs_conf *conf){
	conf->matrix = conf->allocator(conf->r * sizeof(uint8_t*));
    conf->matrix_type = 0;

	for(uint8_t i = 0; i < conf->r; i++){
		conf->matrix[i] = conf->allocator(conf->k * sizeof(uint8_t));
		for(uint8_t j = 0; j < conf->k; j++){
			conf->matrix[i][j] = gf_div((uint8_t)1, i ^ (j + conf->r));
		}
	}
}

__attribute__((always_inline))
static inline void super_fast_cal(int len, const int K, const int R, uint8_t **data, uint8_t **output, uint8_t **Matrix){
	encode_t **data_ptr = (encode_t **)data;
    encode_t **output_ptr = (encode_t **)output;
   
    encode_t res[R][REGION_BLOCKS];
    for(int i = 0; i < R; i++)
        for(int j = 0; j < REGION_BLOCKS; j++)
            res[i][j] = zero();

    for(int i = 0; i < len / sizeof(encode_t); i += REGION_BLOCKS){ 
    	
        for(int k = 0; k < K; k++){ 
            for(int x = 0; x < REGION_BLOCKS; x++){ 
                encode_t tmp = data_ptr[k][i + x]; 
                for(int r = 0; r < R; r++) 
                    res[r][x] = xor_region(res[r][x], multiply_region(tmp, Matrix[r][k]));
            } 
        } 
     		
        for(int r = 0; r < R; r++){ 
            for(int x = 0; x < REGION_BLOCKS; x++){ 
                store(&output_ptr[r][i + x], res[r][x]); 
                res[r][x] = zero(); 
            } 
        } 

    }
    
}

__attribute__((always_inline))
static inline void super_fast_reg(int len, const int K, const int R,
 uint8_t **data, uint8_t **output, uint8_t **Matrix, int *recover_id){
    encode_t **data_ptr = (encode_t **)data;
    encode_t **output_ptr = (encode_t **)output;
   
    encode_t res[R][REGION_BLOCKS];
    for(int i = 0; i < R; i++)
        for(int j = 0; j < REGION_BLOCKS; j++)
            res[i][j] = zero();

    for(int i = 0; i < len / sizeof(encode_t); i += REGION_BLOCKS){ 
        
        for(int k = 0; k < K; k++){ 
            for(int x = 0; x < REGION_BLOCKS; x++){ 
                encode_t tmp = data_ptr[recover_id[k]][i + x]; 
                for(int r = 0; r < R; r++) 
                    res[r][x] = xor_region(res[r][x], multiply_region(tmp, Matrix[r][k]));
            } 
        } 
            
        for(int r = 0; r < R; r++){ 
            for(int x = 0; x < REGION_BLOCKS; x++){ 
                store(&output_ptr[r][i + x], res[r][x]); 
                res[r][x] = zero(); 
            } 
        } 

    }
    
}

void rs_encode(int len, rs_conf *conf, uint8_t **data, 
    uint8_t **output, uint8_t **memory_pre_allocated){

    if(conf->matrix_type == 0){
        switch (conf->k) {
        case 8:
            switch (conf->r) {
                case 1:
                    super_fast_cal(len, 8, 1, data, memory_pre_allocated, conf->matrix);
                    break;
                case 2:
                    super_fast_cal(len, 8, 2, data, memory_pre_allocated, conf->matrix);
                    break;
                case 3:
                    super_fast_cal(len, 8, 3, data, memory_pre_allocated, conf->matrix);
                    break;
                case 4:
                    super_fast_cal(len, 8, 4, data, memory_pre_allocated, conf->matrix);
                    break;
            }
            break;
        case 10:
            switch (conf->r) {
                case 1:
                    super_fast_cal(len, 10, 1, data, memory_pre_allocated, conf->matrix);
                    break;
                case 2:
                    super_fast_cal(len, 10, 2, data, memory_pre_allocated, conf->matrix);
                    break;
                case 3:
                    super_fast_cal(len, 10, 3, data, memory_pre_allocated, conf->matrix);
                    break;
                case 4:
                    super_fast_cal(len, 10, 4, data, memory_pre_allocated, conf->matrix);
                    break;
            }
            break;
        case 12:
            switch (conf->r) {
                case 1:
                    super_fast_cal(len, 12, 1, data, memory_pre_allocated, conf->matrix);
                    break;
                case 2:
                    super_fast_cal(len, 12, 2, data, memory_pre_allocated, conf->matrix);
                    break;
                case 3:
                    super_fast_cal(len, 12, 3, data, memory_pre_allocated, conf->matrix);
                    break;
                case 4:
                    super_fast_cal(len, 12, 4, data, memory_pre_allocated, conf->matrix);
                    break;
            }
            break;
        case 32:
            switch (conf->r) {
                case 16:
                    super_fast_cal(len, 32, 16, data, memory_pre_allocated, conf->matrix);
                    break;
            }
            break;
        default:
            super_fast_cal(len, conf->k, conf->r, data, memory_pre_allocated, conf->matrix);
            break;
        }
        if(output != NULL){
            for(int i = 0; i < conf->n; i++){
                if(i < conf->k) 
                    output[i] = data[i];
                else
                    output[i] = memory_pre_allocated[i - conf->k];
            }
        }
    }
    else if(conf->matrix_type == 1){
        switch(conf->k){
            case 8:
                switch(conf->r){
                    case 1:
                        super_fast_cal(len, 8, 9, data, output, conf->matrix);
                        break;
                    case 2:
                        super_fast_cal(len, 8, 10, data, output, conf->matrix);
                        break;
                    case 3:
                        super_fast_cal(len, 8, 11, data, output, conf->matrix);
                        break;
                    case 4:
                        super_fast_cal(len, 8, 12, data, output, conf->matrix);
                        break;
                }
                break;
            case 10:
                switch(conf->r){
                    case 1:
                        super_fast_cal(len, 10, 11, data, output, conf->matrix);
                        break;
                    case 2:
                        super_fast_cal(len, 10, 12, data, output, conf->matrix);
                        break;
                    case 3:
                        super_fast_cal(len, 10, 13, data, output, conf->matrix);
                        break;
                    case 4:
                        super_fast_cal(len, 10, 14, data, output, conf->matrix);
                        break;
                }
                break;
            case 12:
                switch(conf->r){
                    case 1:
                        super_fast_cal(len, 12, 13, data, output, conf->matrix);
                        break;
                    case 2:
                        super_fast_cal(len, 12, 14, data, output, conf->matrix);
                        break;
                    case 3:
                        super_fast_cal(len, 12, 15, data, output, conf->matrix);
                        break;
                    case 4:
                        super_fast_cal(len, 12, 16, data, output, conf->matrix);
                        break;
                }
                break;
            default:
                super_fast_cal(len, conf->k, conf->k + conf->r, data, output, conf->matrix);
                break;
        }
    }
}

bool rs_check_broken(int x, int *broken_id, int broken_num){
	for(int i = 0; i < broken_num; i++){
		if(broken_id[i] == x) 
            return true;
	}
	return false;
}

void rs_invert_matrix(uint8_t **in_mat, uint8_t **out_mat, int n){
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

void rs_decode(int len, rs_conf *conf, rs_regenerate_context *context,uint8_t **data, uint8_t **output){
	
	switch (conf->k) {
        case 8:
        	super_fast_reg(len, 8, 8, data, output, context->inv_matrix, context->recover_id);
            break;
        case 10:
            super_fast_reg(len, 10, 10, data, output, context->inv_matrix, context->recover_id);
            break;
        case 12:
            super_fast_reg(len, 12, 12, data, output, context->inv_matrix, context->recover_id);
            break;
        default:
            super_fast_reg(len, conf->k, conf->k, data, output, context->inv_matrix, context->recover_id);
            break;
    }
}

void rs_regenerate_context_init(rs_conf *conf, rs_regenerate_context *context, int *broken_id, int broken_num){
	context->broken_id = broken_id;
	context->broken_num = broken_num;
	context->recover_id = conf->allocator(conf->k * sizeof(int));
	uint8_t **matrix;
	matrix = conf->allocator(conf->k * sizeof(uint8_t*));
	int re_num = 0;

	for(int i = 0; i < conf->n; i++){
		if(re_num == conf->k) break;
		if(!rs_check_broken(i, broken_id, broken_num)){
			matrix[re_num] = conf->allocator(conf->k * sizeof(uint8_t));
            if(conf->matrix_type == 0){
                if(i >= conf->k) 
                    memcpy(matrix[re_num], conf->matrix[i - conf->k], conf->k);
                else
                    for(int j = 0; j < conf->k; j++)
                        matrix[re_num][j] = (i == j) ? 1 : 0;
                context->recover_id[re_num] = i;
                re_num++;
            }
            else if(conf->matrix_type == 1){
                memcpy(matrix[re_num], conf->matrix[i], conf->k);
                context->recover_id[re_num] = i;
                re_num++;
            }
		}
	}

	context->inv_matrix = conf->allocator(conf->k * sizeof(uint8_t*));
	for(int i = 0; i < conf->k; i++){
		context->inv_matrix[i] = conf->allocator(conf->k * sizeof(uint8_t));
		for(int j = 0; j < conf->k; j++){
			context->inv_matrix[i][j] = (i==j) ? 1 : 0;
		}
	}

	rs_invert_matrix(matrix, context->inv_matrix, conf->k);
	for(int i = 0; i < conf->k; i++) 
        conf->deallocator(matrix[i]);
	conf->deallocator(matrix);
}

void rs_regenerate(int len, rs_regenerate_context *context, rs_conf *conf, uint8_t **data, uint8_t **output){
	uint8_t **matrix = conf->allocator(context->broken_num * sizeof(uint8_t *));
	uint8_t **mat = conf->allocator(context->broken_num * sizeof(uint8_t *));

	for(int i = 0; i < context->broken_num; i++){
		matrix[i] = conf->allocator(conf->k * sizeof(uint8_t));
		mat[i] = conf->allocator(conf->k * sizeof(uint8_t));

        if(conf->matrix_type == 0){
            if(context->broken_id[i] >= conf->k) 
                memcpy(matrix[i], conf->matrix[context->broken_id[i] - conf->k], conf->k * sizeof(uint8_t));
            else
                for(int j = 0; j < conf->k; j++) 
                    matrix[i][j] = (j == context->broken_id[i]) ? 1 : 0;
        }
        else if(conf->matrix_type == 1){
            memcpy(matrix[i], conf->matrix[context->broken_id[i]], conf->k * sizeof(uint8_t));
        }
		memset(mat[i], 0, conf->k * sizeof(uint8_t));
	}

	for(int i = 0; i < context->broken_num; i++){
		for(int j = 0; j < conf->k; j++){
			for(int k = 0; k < conf->k; k++){
				mat[i][j] ^= gf_mul(matrix[i][k], context->inv_matrix[k][j]);
			}
		}
	}

	switch (conf->k) {
        case 8:
            switch (context->broken_num) {
            	case 1:
                    super_fast_reg(len, 8, 1, data, output, mat, context->recover_id);
                    break;
                case 2:
                    super_fast_reg(len, 8, 2, data, output, mat, context->recover_id);
                    break;
                case 3:
                    super_fast_reg(len, 8, 3, data, output, mat, context->recover_id);
                    break;
                case 4:
                    super_fast_reg(len, 8, 4, data, output, mat, context->recover_id);
                    break;
            }
            break;
        case 10:
            switch (context->broken_num) {
            	case 1:
                    super_fast_reg(len, 10, 1, data, output, mat, context->recover_id);
                    break;
                case 2:
                    super_fast_reg(len, 10, 2, data, output, mat, context->recover_id);
                    break;
                case 3:
                    super_fast_reg(len, 10, 3, data, output, mat, context->recover_id);
                    break;
                case 4:
                    super_fast_reg(len, 10, 4, data, output, mat, context->recover_id);
                    break;
            }
            break;
        case 12:
            switch (context->broken_num) {
            	case 1:
                    super_fast_reg(len, 12, 1, data, output, mat, context->recover_id);
                    break;
                case 2:
                    super_fast_reg(len, 12, 2, data, output, mat, context->recover_id);
                    break;
                case 3:
                    super_fast_reg(len, 12, 3, data, output, mat, context->recover_id);
                    break;
                case 4:
                    super_fast_reg(len, 12, 4, data, output, mat, context->recover_id);
                    break;
            }
            break;
        default:
            super_fast_reg(len, conf->k, context->broken_num, data, output, mat, context->recover_id);
            break;
    }
    for(int i = 0; i < context->broken_num; i++){
    	conf->deallocator(matrix[i]);
    	conf->deallocator(mat[i]);
    }
    conf->deallocator(matrix);
    conf->deallocator(mat);
}

void rs_free_regenerate_context(rs_conf *conf, rs_regenerate_context *context){
    if(context->inv_matrix) {
        for(int i = 0; i < conf->k; i++){
            conf->deallocator(context->inv_matrix[i]);
        }
        conf->deallocator(context->inv_matrix);
    }
    if(context->recover_id)
	    conf->deallocator(context->recover_id);
}

void rs_free_conf(rs_conf *conf){
    if(conf->matrix) {
        if(conf->matrix_type == 0){
            for(int i = 0; i < conf->r; i++){
                conf->deallocator(conf->matrix[i]);
            }
        }
        else if(conf->matrix_type == 1){
            for(int i = 0; i < conf->n; i++){
                conf->deallocator(conf->matrix[i]);
            }
        }
        conf->deallocator(conf->matrix);
    }
}
