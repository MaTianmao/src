//
// Created by msn on 17-10-15.
//

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "reed_solomon.h"
#include <stdio.h>
#include <mm_malloc.h>
#include <math.h>

const int n = 14;
const int k = 10;

#define DATA_SIZE (1 << 26)
#define REGION_SIZE 512

int test_round = 100;
int broken_num = 4;
int broken_id[4];

void test_encode_speed(rs_conf conf, uint8_t **data, 
	uint8_t **encode_data,uint8_t **memory_pre_allocated){
	int r = n - k;	

	clock_t start = clock();
	for(int test_id = 0; test_id < test_round; test_id++){
		rs_encode(DATA_SIZE, &conf, data, encode_data, memory_pre_allocated);
	}
	clock_t end = clock();

	printf("Total Clock Time: %.2fs\n", (end - start) / (double) CLOCKS_PER_SEC);
    printf("Encode Throughput: %.2fGB/s\n",
           test_round * (double) DATA_SIZE* n / ((end - start) / (double) CLOCKS_PER_SEC) / 1024 / 1024 / 1024);
}

void test_decode_speed(rs_conf conf, uint8_t **encode_data, uint8_t **decode_data){
	int r = n - k;

	clock_t start = clock();
	for(int i = 0; i < broken_num; i++){
		broken_id[i] = rand() % 4 + i * ((n + broken_num - 1) / broken_num);
		broken_id[i] = (broken_id[i] < n - 1) ? broken_id[i] : n - 1;
	}

	rs_regenerate_context context;
	rs_regenerate_context_init(&conf, &context, broken_id, broken_num);
	for(int test_id = 0; test_id < test_round; test_id++){
			rs_regenerate(DATA_SIZE, &context, &conf, encode_data, decode_data);
	}
	clock_t end = clock();
	rs_free_regenerate_context(&conf, &context);
	printf("Total Clock Time: %.2fs\n", (end - start) / (double) CLOCKS_PER_SEC);
    printf("Regenerate Throughput: %.2fGB/s\n",
           test_round * (double) DATA_SIZE * n / ((end - start) / (double) CLOCKS_PER_SEC) / 1024 / 1024 / 1024);
}

int main(int argc, char **argv) {
	int r = n - k;
	rs_conf conf;
	rs_init(&conf, n, k, malloc, free);

	uint8_t *data[k];
	uint8_t *memory_pre_allocated[r];
	uint8_t *encode_data[n];
	uint8_t *decode_data[k];

	for(int i = 0; i < k; i++) data[i] = NULL;

	for(int i = 0; i < k; i++){
		posix_memalign((void **)&(data[i]), 64, sizeof(uint8_t) * DATA_SIZE);
		posix_memalign((void **)&(decode_data[i]), 64, sizeof(uint8_t) * DATA_SIZE);
		for(int j = 0; j < DATA_SIZE; j++)
			data[i][j] = (uint8_t)(rand() & 0xff);
	}

	for(int i = 0; i < r; i++)
		posix_memalign((void **)&(memory_pre_allocated[i]), 64, sizeof(uint8_t) * DATA_SIZE);

	for(int i = 0; i < n; i++)
		posix_memalign((void **)&(encode_data[i]), 64, sizeof(uint8_t) * DATA_SIZE);

	printf("/**********Cauchy Performance************/\n");
	cauchy_matrix(&conf);
	test_encode_speed(conf, data, encode_data, memory_pre_allocated);

	test_decode_speed(conf, encode_data, decode_data);
	rs_free_conf(&conf);
	for(int i = 0; i < k; i++) free(decode_data[i]);
	for(int i = 0; i < n; i++) free(encode_data[i]);
	return 0;
}