//
// Created by msn on 17-10-31.
//

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mbr.h"
#include <stdio.h>
#include <mm_malloc.h>
#include <math.h>

const int n = 14;
const int d = 13;
const int k = 10;

#define DATA_SIZE (1 << 17)
#define REGION_SIZE 512

int broken_num = 4;
int broken_id[4];
int test_round = 1;

void test_recover_correctness(mbr_conf conf, uint8_t **data, uint8_t **encode_data,
	uint8_t **memory_pre_allocated, uint8_t **decode_data){

	printf("/**********test recover correctness**********/\n");
	mbr_encode(DATA_SIZE, &conf, data, encode_data, memory_pre_allocated);

	// for(int i = 0; i < d; i++){
	// 	for(int j = 0; j < d; j++) printf("%d ", data[i*d+j][0]);
	// 	printf("\n");
	// }
	// printf("********\n");
	// for(int i = 0; i < n; i++){
	// 	for(int j = 0; j < d; j++) printf("%d ", encode_data[i*d+j][0]);
	// 	printf("\n");
	// }
	// printf("********\n");
	mbr_recover_context context;
	int flag = 0;

	for(int test_id = 0; test_id < test_round; test_id++){
		for(int i = 0; i < broken_num; i++){
			broken_id[i] = rand() % 4 + i * ((n + broken_num - 1) / broken_num);
			broken_id[i] = (broken_id[i] < n - 1) ? broken_id[i] : n - 1;
			//printf("%d\n", broken_id[i]);
		}

		mbr_recover_context_init(&conf, &context, broken_id, broken_num);

		// for(int i = 0; i < conf.n-conf.k; i++){
		// 	for(int j = 0; j < conf.d; j++) printf("%d ",conf.matrix[i][j]);
		// 	printf("\n");
		// }
		// printf("/*************\n");
		// for(int i = 0; i < conf.k; i++){
		// 	for(int j = 0; j < conf.k; j++) printf("%d ",context.inv_matrix[i][j]);
		// 	printf("\n");
		// }
		// printf("**********\n");
		// for(int i = 0; i < conf.k; i++){
		// 	for(int j = 0;  j < conf.d-conf.k; j++) printf("%d ",context.matrix_half[i][j]);
		// 	printf("\n");
		// }
		// printf("***********\n");
		// for(int i = 0; i < conf.k; i++) printf("%d ",context.recover_id[i]);
		// printf("\n");
		// for(int i = 0; i < k; i++){
		// 	for(int j = 0; j < k; j++) printf("%d ", data[i][j]);
		// 	printf("\n");
		// }
		// printf("********\n");
		// for(int i = 0; i < k; i++){
		// 	for(int j = 0; j < k; j++) printf("%d ", decode_data[i][j]);
		// 	printf("\n");
		// }
		// printf("********\n");
		//printf("%d %d\n",conf.d, conf.k);
		mbr_recover(DATA_SIZE, &conf, &context, encode_data, decode_data);
		//printf("recover over\n");
		// for(int i = 0; i < k; i++){
		// 	for(int j = 0; j < k; j++) printf("%d ", decode_data[i][j]);
		// 	printf("\n");
		// }
		// printf("********\n");
		
		for(int i = 0; i < d; i++){
			for(int j = 0; j < d; j++){
				if(memcmp(data[i*d+j], decode_data[i*d+j], DATA_SIZE) != 0){
					printf("Fail to recover!\n");
					flag = 1;
					break;
				}
			}
		}
		mbr_free_recover_context(&conf, &context);
	}
	if(!flag) 
		printf("Successful recover!\n");
}


void test_regenerate_correctness(mbr_conf conf, uint8_t **data, uint8_t **encode_data,
	uint8_t **memory_pre_allocated, uint8_t **regenerate_data){

	printf("/**********test regenerate correctness**********/\n");
	mbr_encode(DATA_SIZE, &conf, data, encode_data, memory_pre_allocated);

	mbr_recover_context context;
	int flag = 0;

	for(int test_id = 0; test_id < test_round; test_id++){
		int bk_id = rand() % n;
		mbr_regenerate_context_init(&conf, &context, bk_id);

		mbr_regenerate(DATA_SIZE, &conf, &context, encode_data, regenerate_data);
		// printf("%d\n",bk_id);
		// for(int i = 0; i < d; i++){
		// 	printf("%d ",encode_data[bk_id*d+i][0]);
		// }
		// printf("********\n");
		// for(int i = 0; i < d; i++){
		// 	printf("%d ",regenerate_data[i][0]);
		// }
		// printf("************\n");	
		for(int i = 0; i < d; i++){
			if(memcmp(encode_data[bk_id*d+i], regenerate_data[i], DATA_SIZE) != 0){
				printf("Fail to regenerate!\n");
				flag = 1;
				break;
			}
		}
		mbr_free_regenerate_context(&conf, &context);
	}
	if(!flag) 
		printf("Successful regenerate!\n");
}

int main(int argc, char **argv) {
	srand(time(0));

	mbr_conf conf;
	mbr_init(&conf, n, k, d, malloc, free);

	uint8_t *data[d*d];
	uint8_t *memory_pre_allocated[(n-k)*d];
	uint8_t *encode_data[n*d];
	uint8_t *decode_data[d*d];
	uint8_t *regenerate_data[d];

	for(int i = 0; i < d; i++){
		for(int j = 0; j < d; j++){
			posix_memalign((void **)(&data[i*d+j]), 64, sizeof(uint8_t) * DATA_SIZE);
			posix_memalign((void **)(&decode_data[i*d+j]), 64, sizeof(uint8_t) * DATA_SIZE);
			for(int x = 0; x < DATA_SIZE; x++){
				if(i <= j) data[i*d+j][x] = (uint8_t)(rand() & 0xff);
				else{
					data[i*d+j][x] = data[j*d+i][x];
				}
				if(i >= k && j >= k) data[i*d+j][x] = 0;
				decode_data[i*d+j][x] = 0;
			}
		}
	}

	for(int i = 0; i < n-k; i++){
		for(int j = 0; j < d; j++){
			posix_memalign((void **)(&memory_pre_allocated[i*d+j]), 64, sizeof(uint8_t) * DATA_SIZE);
		}
	}

	for(int i = 0; i < n; i++){
		for(int j = 0; j < d; j++){
			posix_memalign((void **)(&encode_data[i*d+j]), 64, sizeof(uint8_t) * DATA_SIZE);
		}
	}

	for(int i = 0; i < d; i++){
		posix_memalign((void **)(&regenerate_data[i]), 64, sizeof(uint8_t) * DATA_SIZE);
	}

	test_recover_correctness(conf, data, encode_data, memory_pre_allocated, decode_data);

	test_regenerate_correctness(conf, data, encode_data, memory_pre_allocated, regenerate_data);

	mbr_free_conf(&conf);
	for(int i = 0; i < d * d; i++) free(decode_data[i]);
	for(int i = 0; i < n * d; i++) free(encode_data[i]);
	for(int i = 0; i < d; i++) free(regenerate_data[i]);
	return 0;
}