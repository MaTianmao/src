//
// Created by syd on 17-6-12.
//
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <stdio.h>
#include "msr.h"
#include "arch.h"
//Magic number for CL-MSR.
const static uint8_t u = 3;
const static uint8_t a = 71;
const static uint8_t b = 201;

//To avoid false cache conflict
#define PADDING 3

static int ipow(int base, int exp) {
    int result = 1;
    while (exp) {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

static int get_bit(int z, int y, int q, int t) {
    return z / ipow(q, t - y - 1) % q;
}

static int permute(int z, int remove_bit_id, int added_bit, int q, int t) {
    int result = 0;
    int power = ipow(q, t - 1);
    int i;
    for (i = 0; i < t; i++) {
        if (i != remove_bit_id)
            result = result * q + (z / power) % q;
        else
            result = result * q + added_bit;
        power /= q;
    }
    return result;
}

static void init_companion(msr_conf *conf) {

    int i, z;

    for (i = 0; i < conf->nodes_round_up; i++)
        for (z = 0; z < conf->alpha; z++) {
            int x = i % conf->r;
            int y = i / conf->r;
            conf->node_companion[i * conf->alpha + z] = (uint8_t) (get_bit(z, y, conf->r, conf->groups) + y * conf->r);
            conf->z_companion[i * conf->alpha + z] = (uint8_t) (permute(z, y, x, conf->r, conf->groups));
        }
}

static void init_theta(msr_conf *conf) {
    int i, j;

    int nodes_round_up = conf->groups * conf->r;
    memset(conf->theta, 0, sizeof(uint8_t) * conf->r * nodes_round_up);

    for (i = 0; i < conf->r; i++)
        for (j = 0; j < conf->k; j++) {
            conf->theta[i * nodes_round_up + j] = gf_div(1, (uint8_t) (j ^ (i + conf->n)));
        }

    for (i = 0; i < conf->r; i++) {
        conf->theta[i * nodes_round_up + i + conf->k] = 1;
    }

    for (i = 0; i < conf->r; i++)
        for (j = conf->n; j < nodes_round_up; j++) {
            conf->theta[i * nodes_round_up + j] = gf_div(1, (uint8_t) (j - conf->r ^ (i + conf->n)));
        }

}

static void inline inverse_matrix(uint8_t *matrix, int n, uint8_t *inv) {
    int i = 0, j = 0;

    memset(inv, 0, n * n * sizeof(uint8_t));

    for (i = 0; i < n; i++)
        inv[i * n + i] = 1;

    for (i = 0; i < n; i++) {
        if (!matrix[i * n + i]) {
            for (j = i + 1; j < n; j++)
                if (matrix[j * n + i])
                    break;
            assert(j != n);

            for (int t = 0; t < n; t++) {
                uint8_t tmp = matrix[i * n + t];
                matrix[i * n + t] = matrix[j * n + t];
                matrix[j * n + t] = tmp;

                tmp = inv[i * n + t];
                inv[i * n + t] = inv[j * n + t];
                inv[j * n + t] = tmp;

            }
        }

        uint8_t f = matrix[i * n + i];

        for (j = 0; j < n; j++) {
            matrix[i * n + j] = gf_div(matrix[i * n + j], f);
            inv[i * n + j] = gf_div(inv[i * n + j], f);
        }

        for (j = 0; j < n; j++)
            if (i != j) {
                f = matrix[j * n + i];
                for (int t = 0; t < n; t++) {
                    matrix[j * n + t] ^= gf_mul(matrix[i * n + t], f);
                    inv[j * n + t] ^= gf_mul(inv[i * n + t], f);
                }
            }
    }
}

static int compute_sigma(int z, const int *errors, int error_cnt, int r, int groups) {
    int sigma = 0;
    for (int i = 0; i < error_cnt; i++) {
        if (errors[i] % r == get_bit(z, errors[i] / r, r, groups)) {
            sigma++;
        }
    }
    return sigma;
}

int msr_init(msr_conf *conf, int n, int k, void *(*allocator)(size_t), void (*deallocator)(void *)) {

    int r = n - k;

    if (r <= 1 || r > k || n >= GF_SIZE - 1) {
        return -1;
    }

    gf_init();
    init_arch();

    conf->n = n;
    conf->k = k;
    conf->r = r;

    conf->groups = (n + r - 1) / r;
    conf->alpha = ipow(r, conf->groups);
    conf->beta = conf->alpha / r;

    //Can be improved through padding with a number with is not a power of 2

    conf->allocate = allocator;
    conf->deallocate = deallocator;

    conf->nodes_round_up = conf->groups * conf->r;

    conf->node_companion = allocator(sizeof(uint8_t) * conf->nodes_round_up * conf->alpha);
    conf->z_companion = allocator(sizeof(uint8_t) * conf->nodes_round_up * conf->alpha);
    conf->theta = allocator(sizeof(uint8_t) * conf->r * conf->nodes_round_up);

    init_companion(conf);
    init_theta(conf);

    return 0;
}

void msr_fill_encode_context(msr_encode_context *context, const msr_conf *conf, uint8_t **survived) {
    int survived_cnt = 0;
    int erased_cnt = 0;
    int i, j, z;

    context->is_erased = conf->allocate(conf->nodes_round_up * sizeof(bool));
    context->erase_id = conf->allocate(conf->n * sizeof(int));

    memset(context->is_erased, 0, conf->nodes_round_up * sizeof(bool));

    for (i = 0; i < conf->n; i++)
        if (survived[i] == NULL)
            context->is_erased[i] = true, erased_cnt++;
        else
            context->is_erased[i] = false, survived_cnt++;

    survived_cnt += conf->nodes_round_up - conf->n;
    context->survive_cnt = survived_cnt;
    context->erase_cnt = erased_cnt;

    context->survived = conf->allocate((survived_cnt + 1) * sizeof(int));
    context->erased = conf->allocate((erased_cnt + 1) * sizeof(int));

    survived_cnt = erased_cnt = 0;

    for (i = 0; i < conf->n; i++)
        if (survived[i] == NULL) {
            context->erase_id[i] = erased_cnt;
            context->erased[erased_cnt++] = i;
        } else
            context->survived[survived_cnt++] = i;

    for (i = conf->n; i < conf->nodes_round_up; i++)
        context->survived[survived_cnt++] = i;

    //Trick to prevent from overflow.
    context->survived[survived_cnt] = context->erased[erased_cnt] = 0;

    context->sigmas = conf->allocate(conf->alpha * sizeof(int));
    context->sigma_max = 0;
    for (z = 0; z < conf->alpha; z++) {
        context->sigmas[z] = compute_sigma(z, context->erased, context->erase_cnt, conf->r, conf->groups);
        if (context->sigmas[z] > context->sigma_max)
            context->sigma_max = context->sigmas[z];
    }


    uint8_t *encode_matrix = conf->allocate(conf->r * conf->r * sizeof(uint8_t));
    uint8_t *inv_matrix = conf->allocate(conf->r * conf->r * sizeof(uint8_t));

    for (i = 0; i < conf->r; i++)
        for (j = 0; j < conf->r; j++) {
            if (j < erased_cnt)
                encode_matrix[i * conf->r + j] = conf->theta[i * conf->nodes_round_up + context->erased[j]];
            else
                encode_matrix[i * conf->r + j] = conf->theta[i * conf->nodes_round_up +
                                                             context->survived[j - erased_cnt]];
        }

    inverse_matrix(encode_matrix, conf->r, inv_matrix);

    context->matrix = conf->allocate(conf->r * conf->nodes_round_up * sizeof(uint8_t));

    memset(context->matrix, 0, conf->r * conf->nodes_round_up * sizeof(uint8_t));


    for (i = 0; i < conf->r; i++)
        for (j = 0; j < conf->nodes_round_up; j++) {
            for (int t = 0; t < conf->r; t++) {
                context->matrix[i * conf->nodes_round_up + j] ^= gf_mul(inv_matrix[i * conf->r + t],
                                                                        conf->theta[t * conf->nodes_round_up + j]);
            }
        }

    conf->deallocate(encode_matrix);
    conf->deallocate(inv_matrix);

    context->encoding_buf_size =
            context->erase_cnt * (conf->alpha * REGION_SIZE * sizeof(uint8_t) + PADDING * sizeof(encode_t));
}

void msr_fill_regenerate_context(msr_regenerate_context *context, const msr_conf *conf, int broken) {
    int i, j;

    int nodes_round_up = conf->groups * conf->r;

    context->matrix = conf->allocate(sizeof(uint8_t) * conf->r * nodes_round_up);
    context->u_matrix = conf->allocate(sizeof(uint8_t) * conf->r * nodes_round_up);

    context->broken = broken;

    uint8_t *reg_matrix = conf->allocate(sizeof(uint8_t) * conf->r * conf->r);
    uint8_t *inv_matrix = conf->allocate(sizeof(uint8_t) * conf->r * conf->r);

    int y0 = broken / conf->r;
    int x0 = broken % conf->r;


    for (i = 0; i < conf->r; i++)
        for (j = 0; j < conf->r; j++) {
            int id = y0 * conf->r + j;
            if (j == x0)
                reg_matrix[i * conf->r + j] = conf->theta[i * conf->nodes_round_up + id];
            else
                reg_matrix[i * conf->r + j] = gf_mul(u, conf->theta[i * conf->nodes_round_up + id]);
        }


    inverse_matrix(reg_matrix, conf->r, inv_matrix);

    memset(context->matrix, 0, conf->r * nodes_round_up * sizeof(uint8_t));

    for (i = 0; i < conf->r; i++)
        for (j = 0; j < nodes_round_up; j++) {
            for (int t = 0; t < conf->r; t++) {
                context->matrix[i * nodes_round_up + j] ^= gf_mul(inv_matrix[i * conf->r + t],
                                                                  conf->theta[t * conf->nodes_round_up + j]);

                context->u_matrix[i * nodes_round_up + j] = gf_mul(context->matrix[i * nodes_round_up + j], u);
            }
        }


    conf->deallocate(reg_matrix);
    conf->deallocate(inv_matrix);

    context->z_num = conf->allocate(conf->beta * sizeof(int));
    context->z_pos = conf->allocate(conf->alpha * sizeof(int));

    memset(context->z_pos, -1, conf->alpha * sizeof(int));
    for (int z_id = 0; z_id < conf->beta; z_id++) {
        context->z_num[z_id] =
                (z_id / ipow(conf->r, conf->groups - y0 - 1) * conf->r + x0) * ipow(conf->r, conf->groups - y0 - 1) +
                z_id % ipow(conf->r, conf->groups - y0 - 1);
        context->z_pos[context->z_num[z_id]] = z_id;
    }

    context->z_comp_pos = conf->allocate(conf->n * conf->alpha * sizeof(int));

    for (i = 0; i < conf->n; i++)
        for (int z = 0; z < conf->alpha; z++)
            context->z_comp_pos[i * conf->alpha + z] = context->z_pos[conf->z_companion[i * conf->alpha + z]];

    context->regenerate_buf_size = conf->r * REGION_SIZE * sizeof(uint8_t);

}


__attribute__((always_inline))
static inline void
decode_plane(const msr_encode_context *context, const msr_conf *conf, encode_t **data, encode_t *buf, size_t buf_len,
             int index, int block_size, int z, int erase_cnt) {
    int j = 0, e;
    uint8_t *matrix_ptr = context->matrix;


    for (j = 0; j < context->survive_cnt; j++) {
        int node_id = context->survived[j];

        int z_index = (block_size * z + index) * REGION_BLOCKS;
        int companion = conf->node_companion[node_id * conf->alpha + z];

        int new_node = context->survived[j + 1];
        int new_comp = conf->node_companion[context->survived[j + 1] * conf->alpha + z];
        int new_comp_z = conf->z_companion[context->survived[j + 1] * conf->alpha + z];
        int new_comp_z_index = (block_size * new_comp_z + index) * REGION_BLOCKS;

        if (new_comp >= conf->n)
            new_comp = 0;

        if (new_node >= conf->n)
            new_node = 0;

        uint8_t matrix_tmp[erase_cnt];
        for (e = 0; e < erase_cnt; e++)
            matrix_tmp[e] = matrix_ptr[e * conf->nodes_round_up + node_id];


        if (node_id < conf->n) {
            if (companion < conf->n && companion != node_id && data[companion]) {
                int new_z = conf->z_companion[node_id * conf->alpha + z];
                int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

                for (int w = 0; w < REGION_BLOCKS; w++) {
                    encode_t a_cur = xor_region(data[node_id][z_index + w],
                                                multiply_region(data[companion][new_z_index + w], u));
                    for (e = 0; e < erase_cnt; e++) {
                        buf[e * buf_len + z * REGION_BLOCKS + w] = xor_region(buf[e * buf_len + z * REGION_BLOCKS + w],
                                                                              multiply_region(a_cur, matrix_tmp[e]));
                    }
                    prefetch(&data[new_node][z_index + w]);
                    prefetch(&data[new_comp][new_comp_z_index + w]);
                }
            } else {
                for (int w = 0; w < REGION_BLOCKS; w++) {
                    for (e = 0; e < erase_cnt; e++) {
                        buf[e * buf_len + z * REGION_BLOCKS + w] = xor_region(buf[e * buf_len + z * REGION_BLOCKS + w],
                                                                              multiply_region(
                                                                                      data[node_id][z_index + w],
                                                                                      matrix_tmp[e]));
                    }
                    prefetch(&data[new_node][z_index + w]);
                    prefetch(&data[new_comp][new_comp_z_index + w]);
                }
            }
        } else {
            if (companion < conf->n && data[companion] && companion != node_id) {
                int new_z = conf->z_companion[node_id * conf->alpha + z];
                int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

                for (int w = 0; w < REGION_BLOCKS; w++) {
                    encode_t a_cur = multiply_region(data[companion][new_z_index + w], u);
                    for (e = 0; e < erase_cnt; e++) {
                        buf[e * buf_len + z * REGION_BLOCKS + w] = xor_region(buf[e * buf_len + z * REGION_BLOCKS + w],
                                                                              multiply_region(a_cur, matrix_tmp[e]));
                    }
                }
            }
        }
    }

    for (j = 0; j < erase_cnt; j++) {
        int erased = context->erased[j];
        int companion = conf->node_companion[erased * conf->alpha + z];

        if (erased != companion && !context->is_erased[companion] && companion < conf->n) {
            int new_z = conf->z_companion[erased * conf->alpha + z];
            int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

            for (int w = 0; w < REGION_BLOCKS; w++) {
                buf[j * buf_len + z * REGION_BLOCKS + w] = xor_region(buf[j * buf_len + z * REGION_BLOCKS + w],
                                                                      multiply_region(data[companion][new_z_index + w],
                                                                                      u));
            }
        }
    }

}

__attribute__((always_inline))
static inline void
super_fast_decode_plane(const msr_encode_context *context, const msr_conf *conf, encode_t **data, encode_t *buf,
                        int index, int block_size, int z) {
    switch (conf->alpha) {
        case 4:
            switch (context->erase_cnt) {
                case 1:
                    decode_plane(context, conf, data, buf, 4 * REGION_BLOCKS + PADDING, index, block_size, z, 1);
                    break;
                case 2:
                    decode_plane(context, conf, data, buf, 4 * REGION_BLOCKS + PADDING, index, block_size, z, 2);
                    break;
                case 3:
                    decode_plane(context, conf, data, buf, 4 * REGION_BLOCKS + PADDING, index, block_size, z, 3);
                    break;
                case 4:
                    decode_plane(context, conf, data, buf, 4 * REGION_BLOCKS + PADDING, index, block_size, z, 4);
                    break;
                default:
                    decode_plane(context, conf, data, buf, 4 * REGION_BLOCKS + PADDING, index, block_size, z,
                                 context->erase_cnt);
            }
            break;
        case 8:
            switch (context->erase_cnt) {
                case 1:
                    decode_plane(context, conf, data, buf, 8 * REGION_BLOCKS + PADDING, index, block_size, z, 1);
                    break;
                case 2:
                    decode_plane(context, conf, data, buf, 8 * REGION_BLOCKS + PADDING, index, block_size, z, 2);
                    break;
                case 3:
                    decode_plane(context, conf, data, buf, 8 * REGION_BLOCKS + PADDING, index, block_size, z, 3);
                    break;
                case 4:
                    decode_plane(context, conf, data, buf, 8 * REGION_BLOCKS + PADDING, index, block_size, z, 4);
                    break;
                default:
                    decode_plane(context, conf, data, buf, 8 * REGION_BLOCKS + PADDING, index, block_size, z,
                                 context->erase_cnt);
            }
            break;
        case 9:
            switch (context->erase_cnt) {
                case 1:
                    decode_plane(context, conf, data, buf, 9 * REGION_BLOCKS + PADDING, index, block_size, z, 1);
                    break;
                case 2:
                    decode_plane(context, conf, data, buf, 9 * REGION_BLOCKS + PADDING, index, block_size, z, 2);
                    break;
                case 3:
                    decode_plane(context, conf, data, buf, 9 * REGION_BLOCKS + PADDING, index, block_size, z, 3);
                    break;
                case 4:
                    decode_plane(context, conf, data, buf, 9 * REGION_BLOCKS + PADDING, index, block_size, z, 4);
                    break;
                default:
                    decode_plane(context, conf, data, buf, 9 * REGION_BLOCKS + PADDING, index, block_size, z,
                                 context->erase_cnt);
            }
            break;
        case 16:
            switch (context->erase_cnt) {
                case 1:
                    decode_plane(context, conf, data, buf, 16 * REGION_BLOCKS + PADDING, index, block_size, z, 1);
                    break;
                case 2:
                    decode_plane(context, conf, data, buf, 16 * REGION_BLOCKS + PADDING, index, block_size, z, 2);
                    break;
                case 3:
                    decode_plane(context, conf, data, buf, 16 * REGION_BLOCKS + PADDING, index, block_size, z, 3);
                    break;
                case 4:
                    decode_plane(context, conf, data, buf, 16 * REGION_BLOCKS + PADDING, index, block_size, z, 4);
                    break;
                default:
                    decode_plane(context, conf, data, buf, 16 * REGION_BLOCKS + PADDING, index, block_size, z,
                                 context->erase_cnt);
            }
            break;
        case 27:
            switch (context->erase_cnt) {
                case 1:
                    decode_plane(context, conf, data, buf, 27 * REGION_BLOCKS + PADDING, index, block_size, z, 1);
                    break;
                case 2:
                    decode_plane(context, conf, data, buf, 27 * REGION_BLOCKS + PADDING, index, block_size, z, 2);
                    break;
                case 3:
                    decode_plane(context, conf, data, buf, 27 * REGION_BLOCKS + PADDING, index, block_size, z, 3);
                    break;
                case 4:
                    decode_plane(context, conf, data, buf, 27 * REGION_BLOCKS + PADDING, index, block_size, z, 4);
                    break;
                default:
                    decode_plane(context, conf, data, buf, 27 * REGION_BLOCKS + PADDING, index, block_size, z,
                                 context->erase_cnt);
            }
            break;
        case 64:
            switch (context->erase_cnt) {
                case 1:
                    decode_plane(context, conf, data, buf, 64 * REGION_BLOCKS + PADDING, index, block_size, z, 1);
                    break;
                case 2:
                    decode_plane(context, conf, data, buf, 64 * REGION_BLOCKS + PADDING, index, block_size, z, 2);
                    break;
                case 3:
                    decode_plane(context, conf, data, buf, 64 * REGION_BLOCKS + PADDING, index, block_size, z, 3);
                    break;
                case 4:
                    decode_plane(context, conf, data, buf, 64 * REGION_BLOCKS + PADDING, index, block_size, z, 4);
                    break;
                default:
                    decode_plane(context, conf, data, buf, 64 * REGION_BLOCKS + PADDING, index, block_size, z,
                                 context->erase_cnt);
            }
            break;
        case 81:
            switch (context->erase_cnt) {
                case 1:
                    decode_plane(context, conf, data, buf, 81 * REGION_BLOCKS + PADDING, index, block_size, z, 1);
                    break;
                case 2:
                    decode_plane(context, conf, data, buf, 81 * REGION_BLOCKS + PADDING, index, block_size, z, 2);
                    break;
                case 3:
                    decode_plane(context, conf, data, buf, 81 * REGION_BLOCKS + PADDING, index, block_size, z, 3);
                    break;
                case 4:
                    decode_plane(context, conf, data, buf, 81 * REGION_BLOCKS + PADDING, index, block_size, z, 4);
                    break;
                default:
                    decode_plane(context, conf, data, buf, 81 * REGION_BLOCKS + PADDING, index, block_size, z,
                                 context->erase_cnt);
            }
            break;
        case 256:
            switch (context->erase_cnt) {
                case 1:
                    decode_plane(context, conf, data, buf, 256 * REGION_BLOCKS + PADDING, index, block_size, z, 1);
                    break;
                case 2:
                    decode_plane(context, conf, data, buf, 256 * REGION_BLOCKS + PADDING, index, block_size, z, 2);
                    break;
                case 3:
                    decode_plane(context, conf, data, buf, 256 * REGION_BLOCKS + PADDING, index, block_size, z, 3);
                    break;
                case 4:
                    decode_plane(context, conf, data, buf, 256 * REGION_BLOCKS + PADDING, index, block_size, z, 4);
                    break;
                default:
                    decode_plane(context, conf, data, buf, 256 * REGION_BLOCKS + PADDING, index, block_size, z,
                                 context->erase_cnt);
            }
            break;
        default:
            switch (context->erase_cnt) {
                case 1:
                    decode_plane(context, conf, data, buf, conf->alpha * REGION_BLOCKS + PADDING, index, block_size, z,
                                 1);
                    break;
                case 2:
                    decode_plane(context, conf, data, buf, conf->alpha * REGION_BLOCKS + PADDING, index, block_size, z,
                                 2);
                    break;
                case 3:
                    decode_plane(context, conf, data, buf, conf->alpha * REGION_BLOCKS + PADDING, index, block_size, z,
                                 3);
                    break;
                case 4:
                    decode_plane(context, conf, data, buf, conf->alpha * REGION_BLOCKS + PADDING, index, block_size, z,
                                 4);
                    break;
                default:
                    decode_plane(context, conf, data, buf, conf->alpha * REGION_BLOCKS + PADDING, index, block_size, z,
                                 context->erase_cnt);
            }
            break;
    }
}

static inline void
write_plane(const msr_encode_context *context, const msr_conf *conf, encode_t **data, encode_t *buf, size_t buf_len,
            int index, int block_size, int z) {
    int j;
    for (j = 0; j < context->erase_cnt; j++) {
        int erased = context->erased[j];

        int companion = conf->node_companion[erased * conf->alpha + z];
        int new_z = conf->z_companion[erased * conf->alpha + z];

        int comp_id = context->erase_id[companion];

        if (companion < erased && context->is_erased[companion]) {
            for (int w = 0; w < REGION_BLOCKS; w++) {
                int z_index = z * REGION_BLOCKS + w;
                int new_z_index = new_z * REGION_BLOCKS + w;

                encode_t a_cur = xor_region(multiply_region(buf[j * buf_len + z_index], a),
                                            multiply_region(buf[comp_id * buf_len + new_z_index], b));
                encode_t a_companion = xor_region(multiply_region(buf[j * buf_len + z_index], b),
                                                  multiply_region(buf[comp_id * buf_len + new_z_index], a));

                buf[j * buf_len + z_index] = buf[comp_id * buf_len + new_z_index] = zero();

                store(&data[erased][(block_size * z + index) * REGION_BLOCKS + w], a_cur);
                store(&data[companion][(block_size * new_z + index) * REGION_BLOCKS + w], a_companion);

            }
        } else if (companion == erased || (companion >= conf->n || !context->is_erased[companion])) {
            for (int w = 0; w < REGION_BLOCKS; w++) {
                store(&data[erased][(block_size * z + index) * REGION_BLOCKS + w],
                      buf[j * buf_len + z * REGION_BLOCKS + w]);
                buf[j * buf_len + z * REGION_BLOCKS + w] = zero();
            }
        }

    }
}

__attribute__((always_inline))
static inline void
regenerate_plane(const msr_regenerate_context *context, const msr_conf *conf, encode_t **data, encode_t *buf, int index,
                 int block_size, int z_id, int r) {
    int j;
    int nodes_round_up = conf->groups * conf->r;

    for (j = 0; j < nodes_round_up; j++) {
        int z = context->z_num[z_id];
        int z_index = (block_size * z_id + index) * REGION_BLOCKS;
        int companion = conf->node_companion[j * conf->alpha + z];
        if (j < conf->n && data[j]) {
            if (companion != j && companion < conf->n && data[companion]) {
                int new_z = context->z_pos[conf->z_companion[j * conf->alpha + z]];
                int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;
                for (int w = 0; w < REGION_BLOCKS; w++) {
                    encode_t a_cur = xor_region(data[j][z_index + w],
                                                multiply_region(data[companion][new_z_index + w], u));
                    for (int e = 0; e < r; e++)
                        buf[e * REGION_BLOCKS + w] = xor_region(buf[e * REGION_BLOCKS + w], multiply_region(a_cur,
                                                                                                            context->matrix[
                                                                                                                    e *
                                                                                                                    nodes_round_up +
                                                                                                                    j]));
                }
            } else {
                for (int w = 0; w < REGION_BLOCKS; w++) {
                    for (int e = 0; e < r; e++)
                        buf[e * REGION_BLOCKS + w] = xor_region(buf[e * REGION_BLOCKS + w],
                                                                multiply_region(data[j][z_index + w],
                                                                                context->matrix[e * nodes_round_up +
                                                                                                j]));
                }
            }
        } else if (j != companion && companion < conf->n && data[companion]) {
            int new_z = context->z_pos[conf->z_companion[j * conf->alpha + z]];
            int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;
            for (int w = 0; w < REGION_BLOCKS; w++) {
                for (int e = 0; e < r; e++)
                    buf[e * REGION_BLOCKS + w] = xor_region(buf[e * REGION_BLOCKS + w],
                                                            multiply_region(data[companion][new_z_index + w],
                                                                            context->u_matrix[e * nodes_round_up + j]));
            }
        }
    }

}

__attribute__((always_inline))
static inline void
super_fast_regenerate_plane(const msr_regenerate_context *context, const msr_conf *conf, encode_t **data, encode_t *buf,
                            int index, int block_size, int z) {
    switch (conf->r) {
        case 2:
            regenerate_plane(context, conf, data, buf, index, block_size, z, 2);
            break;
        case 3:
            regenerate_plane(context, conf, data, buf, index, block_size, z, 3);
            break;
        case 4:
            regenerate_plane(context, conf, data, buf, index, block_size, z, 4);
            break;
        default:
            regenerate_plane(context, conf, data, buf, index, block_size, z, conf->r);
    }

}


__attribute__((always_inline))
static inline void
write_regenerate_plane(const msr_regenerate_context *context, const msr_conf *conf, encode_t *output, encode_t *buf,
                       int index, int block_size, int z_id, int r) {
    int j;
    int y0 = context->broken / conf->r;
    for (j = 0; j < r; j++) {
        int z = (z_id / ipow(r, conf->groups - y0 - 1) * conf->r + j) * ipow(conf->r, conf->groups - y0 - 1) +
                z_id % ipow(conf->r, conf->groups - y0 - 1);
        int z_index = (z * block_size + index) * REGION_BLOCKS;


        for (int w = 0; w < REGION_BLOCKS; w++) {
            store(output + z_index + w, buf[j * REGION_BLOCKS + w]);
            buf[j * REGION_BLOCKS + w] = zero();
        }
    }
}

void msr_encode(int len, const msr_encode_context *context, const msr_conf *conf, uint8_t *buf, uint8_t **data,
                uint8_t **output) {
    int z = 0, s = 0;

    int block_size = len / conf->alpha / REGION_SIZE;

    assert(len % (conf->alpha * REGION_SIZE) == 0);

    memset(buf, 0, context->encoding_buf_size * sizeof(uint8_t));

    encode_t *input_ptr[conf->n];
    encode_t *buf_ptr = (encode_t *) buf;

    int erase_cnt = 0;
    for (int i = 0; i < conf->n; i++) {
        if (context->is_erased[i]) {
            data[i] = output[erase_cnt++];
        }
        input_ptr[i] = (encode_t *) data[i];
    }

    for (int index = 0; index < block_size; index++) {
        s = 0;
        while (s <= context->sigma_max) {
            for (z = 0; z < conf->alpha; z++) {
                if (context->sigmas[z] == s) {
                    super_fast_decode_plane(context, conf, input_ptr, buf_ptr, index, block_size, z);
                }
            }

            for (z = 0; z < conf->alpha; z++) {
                if (context->sigmas[z] == s) {
                    write_plane(context, conf, input_ptr, buf_ptr, conf->alpha * REGION_BLOCKS + PADDING, index,
                                block_size, z);
                }
            }
            s++;
        }
    }
}

void msr_regenerate(int len, const msr_regenerate_context *context, const msr_conf *conf, uint8_t *buf, uint8_t **data,
                    uint8_t *output) {
    int block_size = len / conf->beta / REGION_SIZE;

    assert(len % (conf->beta * REGION_SIZE) == 0);

    encode_t *input_ptr[conf->n];
    encode_t *buf_ptr = (encode_t *) buf;
    encode_t *output_ptr = (encode_t *) output;

    for (int i = 0; i < conf->n; i++)
        input_ptr[i] = (encode_t *) data[i];

    memset(buf, 0, sizeof(uint8_t) * context->regenerate_buf_size);


    for (int index = 0; index < block_size; index++) {
        for (int z_id = 0; z_id < conf->beta; z_id++) {
            super_fast_regenerate_plane(context, conf, input_ptr, buf_ptr, index, block_size, z_id);
            write_regenerate_plane(context, conf, output_ptr, buf_ptr, index, block_size, z_id, conf->r);
        }
    }


}

void msr_get_regenerate_offset(int len, const msr_regenerate_context *context, const msr_conf *conf, int *offsets) {
    int y0 = context->broken / conf->r;
    int x0 = context->broken % conf->r;

    for (int i = 0; i < conf->beta; i++) {
        int z = (i / ipow(conf->r, conf->groups - y0 - 1) * conf->r + x0) * ipow(conf->r, conf->groups - y0 - 1) +
                i % ipow(conf->r, conf->groups - y0 - 1);
        offsets[i] = len / conf->alpha * z;
    }
}

void msr_free_conf(msr_conf *conf){ 
    conf->deallocate(conf->theta);
    conf->deallocate(conf->node_companion);
    conf->deallocate(conf->z_companion);
}
void msr_free_encode_context(const msr_conf *conf,msr_encode_context *context){
    conf->deallocate(context->survived);
    conf->deallocate(context->is_erased);
    conf->deallocate(context->matrix);
    conf->deallocate(context->erase_id);
    conf->deallocate(context->erased);
    conf->deallocate(context->sigmas);
}
void msr_free_regenerate_context(const msr_conf *conf,msr_regenerate_context *context){
    conf->deallocate(context->matrix);
    conf->deallocate(context->z_pos);
    conf->deallocate(context->z_num);
    conf->deallocate(context->z_comp_pos);
    conf->deallocate(context->u_matrix);
}