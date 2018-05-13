//
// Created by msn on 17-10-15.
//

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "reed_solomon.h"
#include <stdio.h>
#include <mm_malloc.h>
#include <math.h>

const int n = 14;
const int k = 10;

#define DATA_SIZE (1 << 20)
#define REGION_SIZE 512

int broken_num = 4;
int broken_id[4];
int test_round = 1;

void test_decode_correctness(rs_conf conf, uint8_t **data, 
	uint8_t **encode_data, uint8_t **decode_data){

	int r = n - k;
	printf("/**********test correctness**********/\n");
	rs_encode(DATA_SIZE, &conf, data, encode_data, NULL);

	// for(int i = 0; i < k; i++){
	// 	for(int j = 0; j < k; j++) printf("%d ", data[i][j]);
	// 	printf("\n");
	// }
	// printf("********\n");
	rs_regenerate_context context;
	int flag = 0;

	for(int test_id = 0; test_id < test_round; test_id++){
		for(int i = 0; i < broken_num; i++){
			broken_id[i] = rand() % 4 + i * ((n + broken_num - 1) / broken_num);
			broken_id[i] = (broken_id[i] < n - 1) ? broken_id[i] : n - 1;
			//printf("%d\n", broken_id[i]);
		}

		rs_regenerate_context_init(&conf, &context, broken_id, broken_num);

		// printf("%d %d %d\n", conf.n, conf.k, conf.r);
		// for(int i = 0; i < conf.r; i++){
		// 	for(int j = 0; j < conf.k; j++) printf("%d ",conf.matrix[i][j]);
		// 	printf("\n");
		// }
		// for(int i = 0; i < conf.k; i++){
		// 	for(int j = 0; j < conf.k; j++) printf("%d ",context.inv_matrix[i][j]);
		// 	printf("\n");
		// }
		// printf("**********\n");
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

		rs_decode(DATA_SIZE, &conf, &context, encode_data, decode_data);

		// for(int i = 0; i < k; i++){
		// 	for(int j = 0; j < k; j++) printf("%d ", decode_data[i][j]);
		// 	printf("\n");
		// }
		// printf("********\n");
		
		for(int i = 0; i < k; i++){
			if(memcmp(data[i], decode_data[i], DATA_SIZE) != 0){
				printf("Fail to decode!\n");
				flag = 1;
				break;
			}
		}
		rs_free_regenerate_context(&conf, &context);
	}
	if(!flag) 
		printf("Successful decode!\n");
}


void test_regenerate_correctness(rs_conf conf, uint8_t **data, 
		uint8_t **encode_data, uint8_t **decode_data){
	int r = n - k;
	rs_encode(DATA_SIZE, &conf, data, encode_data, NULL);

	rs_regenerate_context context;
	int flag = 0;

	for(int test_id = 0; test_id < test_round; test_id++){
		for(int i = 0; i < broken_num; i++){
			broken_id[i] = rand() % 4 + i * ((n + broken_num - 1) / broken_num);
			broken_id[i] = (broken_id[i] < n - 1) ? broken_id[i] : n - 1;
		}

		rs_regenerate_context_init(&conf, &context, broken_id, broken_num);
		rs_regenerate(DATA_SIZE, &context, &conf, encode_data, decode_data);
		
		for(int i = 0; i < broken_num; i++){
			if(memcmp(encode_data[broken_id[i]], decode_data[i], DATA_SIZE) != 0){
				printf("Fail to decode!\n");
				flag = 1;
				break;
			}
		}
	}
	if(!flag) 
		printf("Successful decode!\n");
	rs_free_regenerate_context(&conf, &context);
}
int main(int argc, char **argv) {
	srand(time(0));

	int r = n - k;
	rs_conf conf;
	rs_init(&conf, n, k, malloc, free);

	uint8_t *data[k];
	uint8_t *encode_data[n];
	uint8_t *decode_data[k];

	for(int i = 0; i < k; i++) data[i] = NULL;

	for(int i = 0; i < k; i++){
		posix_memalign((void **)&(data[i]), 64, sizeof(uint8_t) * DATA_SIZE);
		posix_memalign((void **)&(decode_data[i]), 64, sizeof(uint8_t) * DATA_SIZE);
		for(int j = 0; j < DATA_SIZE; j++){
			data[i][j] = (uint8_t)(rand() & 0xff);
			decode_data[i][j] = 0;
		}
	}

	for(int i = 0; i < n; i++)
		posix_memalign((void **)&(encode_data[i]), 64, sizeof(uint8_t) * DATA_SIZE);
	printf("/*********Vandermonde Matrix***********/\n");
	vandermonde_matrix(&conf);
	test_decode_correctness(conf, data, encode_data, decode_data);
	rs_free_conf(&conf);

	printf("/*******Regenerate Correctness********/\n");
	vandermonde_matrix(&conf);
	test_regenerate_correctness(conf, data, encode_data, decode_data);
	rs_free_conf(&conf);
	for(int i = 0; i < k; i++) free(data[i]);
	for(int i = 0; i < k; i++) free(decode_data[i]);
	for(int i = 0; i < n; i++) free(encode_data[i]);
	return 0;
}