//
// Created by msn on 17-10-14.
//
#ifndef REED_SOLOMON_H
#define REED_SOLOMON_H

#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
struct rs_conf_t {
	int n;
	int k;
	int r;

	void *(*allocator)(size_t);

	void (*deallocator)(void *);

	uint8_t **matrix; //matrix of coeficients for encoding
	int matrix_type; //0 cauchy, 1 vandermonde, 2 cauchy-vandermonde
};
typedef struct rs_conf_t rs_conf;

struct rs_regenerate_context_t {
	uint8_t **inv_matrix;
	int *broken_id;
	int broken_num;
	int *recover_id;
};
typedef struct rs_regenerate_context_t rs_regenerate_context;


void rs_init(rs_conf *conf, int n, int k, void *(allocator)(size_t), void (*deallocator)(void *));
void vandermonde_matrix(rs_conf *conf);
void cauchy_matrix(rs_conf *conf);
void rs_encode(int len, rs_conf *conf, uint8_t **data, uint8_t **output, uint8_t **memory_pre_allocated);
void rs_decode(int len, rs_conf *conf, rs_regenerate_context *context,uint8_t **data, uint8_t **output);
void rs_regenerate(int len, rs_regenerate_context *context, rs_conf *conf, uint8_t **data, uint8_t **output);
void rs_regenerate_context_init(rs_conf *conf, rs_regenerate_context *context,int *broken_id, int broken_num);
void rs_free_regenerate_context(rs_conf *conf, rs_regenerate_context *context);
void rs_free_conf(rs_conf *conf);

#endif //REED_SOLOMON_H
