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

#define DATA_SIZE (1 << 22)
#define REGION_SIZE 512

int broken_num = 4;
int broken_id[4];
int test_round = 100;

void recover_speed(mbr_conf conf, uint8_t **data, uint8_t **encode_data,
	uint8_t **memory_pre_allocated, uint8_t **decode_data){

	mbr_recover_context context;
	for(int i = 0; i < broken_num; i++){
		broken_id[i] = rand() % 4 + i * ((n + broken_num - 1) / broken_num);
		broken_id[i] = (broken_id[i] < n - 1) ? broken_id[i] : n - 1;
	}

	mbr_recover_context_init(&conf, &context, broken_id, broken_num);

	clock_t start = clock();
	for(int test_id = 0; test_id < test_round; test_id++){
		mbr_recover(DATA_SIZE, &conf, &context, encode_data, decode_data);
	}
	clock_t end = clock();
	mbr_free_recover_context(&conf, &context);

	printf("Total Clock Time: %.2fs\n", (end - start) / (double) CLOCKS_PER_SEC);
    printf("Recover Throughput: %.2fGB/s\n",
           (double) DATA_SIZE / ((end - start) / (double) CLOCKS_PER_SEC) / 1024 / 1024 / 1024 * test_round * d * d);

}

void encode_speed(mbr_conf conf, uint8_t **data, uint8_t **encode_data,
	uint8_t **memory_pre_allocated, uint8_t **decode_data){

	clock_t start = clock();
	for(int test_id = 0; test_id < test_round; test_id++){
		mbr_encode(DATA_SIZE, &conf, data, encode_data, memory_pre_allocated);
	}
	clock_t end = clock();

	printf("Total Clock Time: %.2fs\n", (end - start) / (double) CLOCKS_PER_SEC);
    printf("Encode Throughput: %.2fGB/s\n",
           (double) DATA_SIZE / ((end - start) / (double) CLOCKS_PER_SEC) / 1024 / 1024 / 1024 * test_round * n * d);

}

void regenerate_speed(mbr_conf conf, uint8_t **data, uint8_t **encode_data,
	uint8_t **memory_pre_allocated, uint8_t **regenerate_data){

	mbr_regenerate_context context;
	int bk_id = rand() % n;

	mbr_regenerate_context_init(&conf, &context, bk_id);

	clock_t start = clock();
	for(int test_id = 0; test_id < test_round; test_id++){
		mbr_regenerate(DATA_SIZE, &conf, &context, encode_data, regenerate_data);
	}
	clock_t end = clock();
	mbr_free_regenerate_context(&conf, &context);

	printf("Total Clock Time: %.2fs\n", (end - start) / (double) CLOCKS_PER_SEC);
    printf("Regenerate Throughput: %.2fGB/s\n",
            (double) DATA_SIZE / ((end - start) / (double) CLOCKS_PER_SEC) / 1024 / 1024 / 1024 *test_round *d);

}


// void test_regenerate_correctness(rs_conf conf, uint8_t **data, uint8_t **encode_data,
// 	uint8_t **memory_pre_allocated, uint8_t **decode_data){
// 	int r = n - k;
// 	rs_encode(DATA_SIZE, &conf, data, encode_data, memory_pre_allocated);

// 	rs_regenerate_context context;
// 	int flag = 0;

// 	for(int test_id = 0; test_id < test_round; test_id++){
// 		for(int i = 0; i < broken_num; i++){
// 			broken_id[i] = rand() % 4 + i * ((n + broken_num - 1) / broken_num);
// 			broken_id[i] = (broken_id[i] < n - 1) ? broken_id[i] : n - 1;
// 		}

// 		rs_regenerate_context_init(&conf, &context, broken_id, broken_num);
// 		rs_regenerate(DATA_SIZE, &context, &conf, encode_data, decode_data);
		
// 		for(int i = 0; i < broken_num; i++){
// 			if(memcmp(encode_data[broken_id[i]], decode_data[i], DATA_SIZE) != 0){
// 				printf("Fail to decode!\n");
// 				flag = 1;
// 				break;
// 			}
// 		}
// 	}
// 	if(!flag) 
// 		printf("Successful decode!\n");
// 	rs_free_regenerate_context(&conf, &context);
// }
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

	encode_speed(conf, data, encode_data, memory_pre_allocated, decode_data);
	recover_speed(conf, data, encode_data, memory_pre_allocated, decode_data);
	regenerate_speed(conf, data, encode_data, memory_pre_allocated, regenerate_data);

	mbr_free_conf(&conf);
	for(int i = 0; i < d * d; i++) free(decode_data[i]);
	for(int i = 0; i < n * d; i++) free(encode_data[i]);
	for(int i = 0; i < d; i++) free(regenerate_data[i]);
	return 0;
}