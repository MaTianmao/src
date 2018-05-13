//
// Created by syd on 16-10-23.
//
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "msr.h"
#include <stdio.h>
#include <mm_malloc.h>

const int n = 14;
const int k = 10;

#define DATA_SIZE (k * (1 << 22))
#define REGION_SIZE 512

int main(int argc, char **argv) {
    srand(time(0));

    int r = n-k;
    msr_conf conf;
    msr_init(&conf,n,k,malloc,free);


    uint8_t *data[n];
    uint8_t *memory_pre_allocated[r];

    for(int i=0;i<n;i++)
        data[i] = NULL;

    for(int i=0;i<k;i++) {
        posix_memalign((void **)&data[i],64,sizeof(uint8_t) * DATA_SIZE / k);
        for(int j=0;j<DATA_SIZE/k;j++)
            data[i][j] = (uint8_t)(rand() & 0xff);
    }

    for(int i=0;i<r;i++) {
        posix_memalign((void **)&memory_pre_allocated[i],64,sizeof(uint8_t) * DATA_SIZE / k);
    }

    msr_encode_context context;
    msr_fill_encode_context(&context, &conf, data);


    uint8_t *buf ;
    posix_memalign((void **)&buf,64,context.encoding_buf_size * sizeof(uint8_t));

    msr_encode(DATA_SIZE/k,&context,&conf,buf,data,memory_pre_allocated);

    msr_free_encode_context(&conf,&context);


    for (int i = 0; i < n; i++) {
        printf("Original %d: ", i);
        for (int s = 0; s < 8; s++)
            printf("%x ", data[i][s]);
        printf("\n");
    }


    printf("-----------Begin to test decode-----------\n");
    int test_turn = 100;

    for (int t = 0; t < test_turn; t++) {
        printf("Turn %d:\n", t);

        uint8_t *input[n];
        for (int j = 0; j < n; j++)
            input[j] = NULL;

        int survive_cnt =  k + rand() % r;

        int ok_cnt = 0;
        while (ok_cnt < survive_cnt) {
            int ok_id = rand() % n;
            if (!input[ok_id]) {
                input[ok_id] = data[ok_id];
                ok_cnt++;
            }
        }

        for(int i=0;i<r;i++) {
            posix_memalign((void **)&(memory_pre_allocated[i]),64,sizeof(uint8_t) * DATA_SIZE / k);
        }
        printf("mtm\n");


        msr_encode_context encode_context;
        msr_fill_encode_context(&encode_context, &conf, input);


        msr_encode(DATA_SIZE/k,&encode_context,&conf,buf,input,memory_pre_allocated);


        for(int j=0;j<n;j++)
            assert(!memcmp(input[j],data[j],DATA_SIZE / k));

        printf("Check OK!\n");

        for(int i=0;i<r;i++) {
            free(memory_pre_allocated[i]);
        }

        msr_free_encode_context(&conf,&encode_context);

    }

    free(buf);



    printf("-----------Begin to test regenerate-----------\n");
    test_turn = n;


    for (int i = 0; i < test_turn; i++) {
        printf("Turn %d:\n", i);

        msr_regenerate_context regenerate_context;
        msr_fill_regenerate_context(&regenerate_context,&conf,i);

        uint8_t *regenerate_buf;

        posix_memalign((void **)&regenerate_buf,64,sizeof(uint8_t)* regenerate_context.regenerate_buf_size);

        uint8_t *input[n];
        int offsets[conf.beta];

        msr_get_regenerate_offset(DATA_SIZE/k, &regenerate_context, &conf, offsets);
            if(j!=i){
                posix_memalign((void **)(&input[j]),64,sizeof(uint8_t) * DATA_SIZE / k / r);
                for(int z=0;z<conf.beta;z++){
                    memcpy(input[j] + DATA_SIZE/k/conf.alpha * z,data[j] + offsets[z], DATA_SIZE/k/conf.alpha * sizeof(uint8_t));
                }

            }else
                input[j] = NULL;

        uint8_t * memory;

        posix_memalign((void **)(&memory),64,sizeof(uint8_t) * DATA_SIZE / k);

        msr_regenerate(DATA_SIZE/k/r,&regenerate_context,&conf,regenerate_buf,input,memory);



        assert(!memcmp(memory,data[i],DATA_SIZE/k));

        printf("Check OK!\n");

        free(memory);


        for(int j=0;j<n;j++)
            if(input[j] != NULL)
                free(input[j]);

        free(regenerate_buf);
        msr_free_regenerate_context(&conf,&regenerate_context);

    }

    for(int i=0;i<n;i++)
        free(data[i]);

    msr_free_conf(&conf);

    return 0;
}