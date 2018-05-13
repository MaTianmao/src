//
// Created by syd on 16-11-11.
//
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "msr.h"
#include <stdio.h>
#include <mm_malloc.h>
#include <malloc.h>

const int n = 14;
const int k = 10;


#define REGION_SIZE 512

#define TEST_LOOP (1)




int main() {

    const int r = n - k;

    msr_conf conf;
    msr_init(&conf, n, k, malloc, free);

    int stripe_size = k * conf.alpha * REGION_SIZE;

    int DATA_SIZE = k * (1<<27);
    DATA_SIZE-= DATA_SIZE%stripe_size;

    uint8_t *data[n];
    uint8_t *memory_pre_allocated[r];

    printf("n:%d r:%d\n", n, r);

    for (int i = 0; i < n; i++)
        data[i] = NULL;

    for (int i = 0; i < k; i++) {
        posix_memalign((void **) &(data[i]), 64, sizeof(uint8_t) * DATA_SIZE / k);
        //Warm the memory up.
        memset(data[i], 0xaa, sizeof(uint8_t) * DATA_SIZE / k);
    }

    for (int i = 0; i < r; i++) {
        posix_memalign((void **) &(memory_pre_allocated[i]), 64, sizeof(uint8_t) * DATA_SIZE / k);
        //Warm the memory up.
        memset(memory_pre_allocated[i], 0x00, sizeof(uint8_t) * DATA_SIZE / k);
    }


    clock_t start = clock();

    for (int loop = 0; loop < TEST_LOOP; loop++) {
        for (int i = 0; i < r; i++)
            data[i + k] = NULL;

        msr_encode_context context;
        msr_fill_encode_context(&context, &conf, data);

        uint8_t *buf;
        posix_memalign((void **) &buf, 64, context.encoding_buf_size * sizeof(uint8_t));

        msr_encode(DATA_SIZE / k, &context, &conf, buf, data, memory_pre_allocated);

        free(buf);
        msr_free_encode_context(&conf,&context);
    }

    printf("Total Clock Time: %.2fs\n", (clock() - start) / (double) CLOCKS_PER_SEC);

    printf("Encode Throughput: %.2fMB/s\n",
           TEST_LOOP * (double) DATA_SIZE / k * n / ((clock() - start) / (double) CLOCKS_PER_SEC) * 1e-6);


    start = clock();

    for (int loop = 0; loop < TEST_LOOP; loop++) {

        uint8_t *input[n];

        for (int j = 0; j < n; j++)
            input[j] = NULL;

        int survive_cnt =  k;

        int ok_cnt = 0;
        while (ok_cnt < survive_cnt) {
            int ok_id = rand() % n;
            if (!input[ok_id]) {
                input[ok_id] = data[ok_id];
                ok_cnt++;
            }
        }
        msr_encode_context context;
        msr_fill_encode_context(&context, &conf, input);

        for(int i=0;i<r;i++) {
            posix_memalign((void *)&(memory_pre_allocated[i]),64,sizeof(uint8_t) * DATA_SIZE / k);
        }

        uint8_t *buf;
        posix_memalign((void **) &buf, 64, context.encoding_buf_size * sizeof(uint8_t));

        msr_encode(DATA_SIZE / k, &context, &conf, buf, input, memory_pre_allocated);

        for(int i=0;i<r;i++) {
           free(memory_pre_allocated[i]);
        }

        free(buf);
        msr_free_encode_context(&conf,&context);
    }

    printf("Total Clock Time: %.2fs\n", (clock() - start) / (double) CLOCKS_PER_SEC);

    printf("Decode Throughput: %.2fMB/s\n",
           TEST_LOOP * (double) DATA_SIZE / k * n / ((clock() - start) / (double) CLOCKS_PER_SEC) * 1e-6);


    int error = 1;

    msr_regenerate_context context;
    msr_fill_regenerate_context(&context, &conf, error);

    uint8_t *buf;

    posix_memalign((void **) &buf, 64, sizeof(uint8_t) * context.regenerate_buf_size);

    uint8_t *input[n];
    int offsets[conf.beta];

    msr_get_regenerate_offset(DATA_SIZE / k, &context, &conf, offsets);


    for (int j = 0; j < n; j++) {
        if (j != error) {
            posix_memalign((void **) (&input[j]), 64, sizeof(uint8_t) * DATA_SIZE / k / r);
            for (int z = 0; z < conf.beta; z++) {
                memcpy(input[j] + DATA_SIZE / k / conf.alpha * z, data[j] + offsets[z],
                       DATA_SIZE / k / conf.alpha * sizeof(uint8_t));
            }

        } else
            input[j] = NULL;
    }

    uint8_t *memory;

    posix_memalign((void **) (&memory), 64, sizeof(uint8_t) * DATA_SIZE / k);

    start = clock();

    for (int i = 0; i < TEST_LOOP; i++)
        msr_regenerate(DATA_SIZE / k / r, &context, &conf, buf, input, memory);

    printf("Total Clock Time: %.2fs\n", (clock() - start) / (double) CLOCKS_PER_SEC);

    printf("Regenerate Throughput: %.2fMB/s\n",
           TEST_LOOP * (double) DATA_SIZE / k / ((clock() - start) / (double) CLOCKS_PER_SEC) * 1e-6);

    free(memory);
    free(buf);


    msr_free_regenerate_context(&conf,&context);
    for (int j = 0; j < n; j++)
        if (input[j] != NULL)
            free(input[j]);


    for (int j = 0; j < n; j++)
        free(data[j]);

    msr_free_conf(&conf);

}