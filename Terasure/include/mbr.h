//
// Created by msn on 17-10-31.
//
#ifndef MBR_H
#define MBR_H

#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
struct mbr_conf_t {
	int n;
	int k;
	int d;

	void *(*allocator)(size_t);

	void (*deallocator)(void *);

	uint8_t **matrix; //matrix of coeficients for encoding
};
typedef struct mbr_conf_t mbr_conf;
struct mbr_recover_context_t{
	uint8_t **inv_matrix;
	uint8_t **matrix_half;
	int *broken_id;
	int broken_num;
	int *recover_id;
};
typedef struct mbr_recover_context_t mbr_recover_context;

struct mbr_regenerate_context_t {//regenerate single node
	uint8_t **inv_matrix;
	int broken_id;
	uint8_t *broken_vector;
};
typedef struct mbr_regenerate_context_t mbr_regenerate_context;


void mbr_init(mbr_conf *conf, int n, int k, int d, void *(allocator)(size_t), void (*deallocator)(void *));
void mbr_encode(int len, mbr_conf *conf, uint8_t **data, uint8_t **output, uint8_t **memory_pre_allocated);
void mbr_recover(int len, mbr_conf *conf, mbr_recover_context *context, uint8_t **data, uint8_t **output);
void mbr_recover_context_init(mbr_conf *conf, mbr_recover_context *context, int *broken_id, int broken_num);
void mbr_regenerate(int len, mbr_conf *conf, mbr_regenerate_context *context, uint8_t **data, uint8_t **output);
void mbr_regenerate_context_init(mbr_conf *conf, mbr_regenerate_context *context,int broken_id);
void mbr_free_regenerate_context(mbr_conf *conf, mbr_regenerate_context *context);
void mbr_free_conf(mbr_conf *conf);
void mbr_free_recover_context(mbr_conf *conf, mbr_recover_context *context);
#endif //MBR_H
