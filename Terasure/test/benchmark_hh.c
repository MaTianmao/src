//
// Created by gty on 17-8-29.
//

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "hh.h"
#include <stdio.h>
#include <mm_malloc.h>
#include <malloc.h>

const int n = 14;
const int k = 10;

#define REGION_SIZE 512

#define TEST_LOOP (1)

int main() {

    const int r = n - k;

    hh_conf conf;
    hh_init(&conf, n, k, malloc, free);

    int stripe_size = k;

    int DATA_SIZE = k * (1 << 27);
    DATA_SIZE -= DATA_SIZE % stripe_size;

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
        hh_encode(DATA_SIZE / k, &conf, data, memory_pre_allocated);
    }
    printf("Total Clock Time: %.2fs\n", (clock() - start) / (double) CLOCKS_PER_SEC);
    printf("Encode Throughput: %.2fMB/s\n",
           TEST_LOOP * (double) DATA_SIZE / k * n / ((clock() - start) / (double) CLOCKS_PER_SEC) * 1e-6);

    int error = k;
    hh_regenerate_context context;
    uint8_t *input[n];
    hh_get_regenerate_offset(DATA_SIZE / k, &context, &conf, error);
    for (int j = 0; j < n; j++) {
        if (j != error) {
            posix_memalign((void **) (&input[j]), 64, sizeof(uint8_t) * context.size[j]);
            memcpy(input[j], data[j] + context.offset[j], context.size[j] * sizeof(uint8_t));
        } else
            input[j] = NULL;
    }
    uint8_t *memory;
    posix_memalign((void **) (&memory), 64, sizeof(uint8_t) * DATA_SIZE / k);
    start = clock();
    for (int i = 0; i < TEST_LOOP; i++)
        hh_regenerate(DATA_SIZE / k, &context, &conf, input, memory);
    printf("Total Clock Time: %.2fs\n", (clock() - start) / (double) CLOCKS_PER_SEC);
    printf("Regenerate Throughput: %.2fMB/s\n",
           TEST_LOOP * (double) DATA_SIZE / k / ((clock() - start) / (double) CLOCKS_PER_SEC) * 1e-6);
    free(memory);
    hh_free_regenerate_context(&conf, &context);
    for (int j = 0; j < n; j++)
        if (input[j] != NULL)
            free(input[j]);


    for (int j = 0; j < n; j++)
        free(data[j]);

    hh_free_conf(&conf);
    return 0;
}

