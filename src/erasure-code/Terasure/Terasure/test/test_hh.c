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

const int n = 14;
const int k = 10;

#define DATA_SIZE (k * (1 << 22))
#define REGION_SIZE 512

int main(int argc, char **argv) {
    srand(time(0));

    int r = n - k;
    hh_conf conf;
    hh_init(&conf, n, k, malloc, free);


    uint8_t *data[n];
    uint8_t *memory_pre_allocated[r];

    for (int i = 0; i < n; i++)
        data[i] = NULL;

    for (int i = 0; i < k; i++) {
        posix_memalign((void **) &data[i], 64, sizeof(uint8_t) * DATA_SIZE / k);
        for (int j = 0; j < DATA_SIZE / k; j++)
            data[i][j] = (uint8_t) (rand() & 0xff);
    }
    for (int i = 0; i < r; i++) {
        posix_memalign((void **) &memory_pre_allocated[i], 64, sizeof(uint8_t) * DATA_SIZE / k);
    }
    hh_encode(DATA_SIZE / k, &conf, data, memory_pre_allocated);
    for (int i = 0; i < n; i++) {
        printf("Original %d: ", i);
        for (int s = 0; s < 8; s++)
            printf("%x ", data[i][s]);
        printf("\n");
    }

    printf("-----------Begin to test regenerate-----------\n");
    int test_turn = n;

    for (int i = 0; i < test_turn; i++) {
        printf("Turn %d:\n", i);

        hh_regenerate_context regenerate_context;


        uint8_t *input[n];

        hh_get_regenerate_offset(DATA_SIZE / k, &regenerate_context, &conf, i);
        for (int j = 0; j < n; j++) {
            if (j != i) {
                posix_memalign((void **) (&input[j]), 64, sizeof(uint8_t) * regenerate_context.size[j]);
                memcpy(input[j], data[j] + regenerate_context.offset[j], regenerate_context.size[j] * sizeof(uint8_t));
            } else
                input[j] = NULL;
        }
        uint8_t *memory;
        posix_memalign((void **) (&memory), 64, sizeof(uint8_t) * DATA_SIZE / k);
        hh_regenerate(DATA_SIZE / k, &regenerate_context, &conf, input, memory);
        assert(!memcmp(memory, data[i], DATA_SIZE / k));
        printf("Check OK!\n");
        free(memory);
        for (int j = 0; j < n; j++)
            if (input[j] != NULL)
                free(input[j]);
        hh_free_regenerate_context(&conf, &regenerate_context);
    }
    for (int i = 0; i < n; i++)
        free(data[i]);
    hh_free_conf(&conf);
    return 0;
}