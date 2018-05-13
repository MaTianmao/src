// -*- mode:C++; tab-width:8; c-basic-offset:2; indent-tabs-mode:t -*- 
// vim: ts=8 sw=2 smarttab
/*
 * Ceph distributed storage system
 *
 * Copyright (C) 2013,2014 Cloudwatt <libre.licensing@cloudwatt.com>
 * Copyright (C) 2014 Red Hat <contact@redhat.com>
 *
 * Author: Loic Dachary <loic@dachary.org>
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 * 
 */

#include "common/debug.h"
#include "ErasureCodeTerasure.h"

using namespace std;

extern "C"{
  #include <stdlib.h>
  #include <mm_malloc.h>
}

#define LARGEST_VECTOR_WORDSIZE 16

#define dout_context g_ceph_context
#define dout_subsys ceph_subsys_osd
#undef dout_prefix
#define dout_prefix _prefix(_dout)

static ostream& _prefix(std::ostream* _dout)
{
  return *_dout << "ErasureCodeTerasure: ";
}


int ErasureCodeTerasure::init(ErasureCodeProfile& profile, ostream *ss)
{
  int err = 0;
  dout(10) << "technique=" << technique << dendl;
  profile["technique"] = technique;
  err |= parse(profile, ss);
  if (err)
    return err;
  prepare();
  return ErasureCode::init(profile, ss);
}

int ErasureCodeTerasure::parse(ErasureCodeProfile &profile,
			       ostream *ss)
{
  int err = ErasureCode::parse(profile, ss);
  err |= to_int("k", profile, &k, DEFAULT_K, ss);
  err |= to_int("m", profile, &m, DEFAULT_M, ss);
  if (chunk_mapping.size() > 0 && (int)chunk_mapping.size() != k + m) {
    *ss << "mapping " << profile.find("mapping")->second
	<< " maps " << chunk_mapping.size() << " chunks instead of"
	<< " the expected " << k + m << " and will be ignored" << std::endl;
    chunk_mapping.clear();
    err = -EINVAL;
  }
  err |= sanity_check_k(k, ss);
  return err;
}

unsigned int ErasureCodeTerasure::get_chunk_size(unsigned int object_size) const
{
    return object_size / k;
}

int ErasureCodeTerasure::encode_chunks(const set<int> &want_to_encode,
				       map<int, bufferlist> *encoded)
{
  char *chunks[k + m];
  for (int i = 0; i < k + m; i++)
    chunks[i] = (*encoded)[i].c_str();
  Terasure_encode(&chunks[0], &chunks[k], (*encoded)[0].length());
  return 0;
}

int ErasureCodeTerasure::decode_chunks(const set<int> &want_to_read,
				       const map<int, bufferlist> &chunks,
				       map<int, bufferlist> *decoded)
{
  unsigned blocksize = (*chunks.begin()).second.length();
  int erasures[k + m + 1];
  int erasures_count = 0;
  char *data[k];
  char *coding[m];
  for (int i =  0; i < k + m; i++) {
    if (chunks.find(i) == chunks.end()) {
      erasures[erasures_count] = i;
      erasures_count++;
    }
    if (i < k)
      data[i] = (*decoded)[i].c_str();
    else
      coding[i - k] = (*decoded)[i].c_str();
  }
  erasures[erasures_count] = -1;

  assert(erasures_count > 0);
  return Terasure_decode(erasures, data, coding, blocksize);
}

bool ErasureCodeTerasure::is_prime(int value)
{
  int prime55[] = {
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
    73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,
    151,157,163,167,173,179,
    181,191,193,197,199,211,223,227,229,233,239,241,251,257
  };
  int i;
  for (i = 0; i < 55; i++)
    if (value == prime55[i])
      return true;
  return false;
}

// // 
// // ErasureCodeTerasureReedSolomonVandermonde
// //
// void ErasureCodeTerasureReedSolomonVandermonde::Terasure_encode(char **data,
//                                                                 char **coding,
//                                                                 int blocksize)
// {
//     char *output[k+m];
//     for(int i = 0; i < k + m; i++)
//         posix_memalign((void **)&output[i], 64, sizeof(char) * blocksize);
//     rs_encode(blocksize, &conf, (uint8_t **)data, (uint8_t **)output, NULL);
//     for(int i = 0; i < k + m; i++){
//         if(i < k) memcpy(data[i], output[i], blocksize);
//         else memcpy(coding[i], output[i - k], blocksize);
//         free(output[i]);
//     }

// }

// int ErasureCodeTerasureReedSolomonVandermonde::Terasure_decode(int *erasures,
//                                                                 char **data,
//                                                                 char **coding,
//                                                                 int blocksize)
// {
//     rs_regenerate_context context;
//     int erasure_num = 0;
//     for(int i = 0; erasures[i] != -1; i++) erasure_num++;
//     if(erasure_num > m) return -1;
//     rs_regenerate_context_init(&conf, &context, erasures, erasure_num);
//     char *input[k+m];
//     char *output[m];
//     for(int i = 0; i < k + m; i++){
//         input[i] = (i < k) ? data[i]:coding[i-k];
//     }
//     for(int i = 0; i < erasure_num; i++){
//         if(erasures[i] < k) output[i] = data[erasures[i]];
//         else output[i] = coding[erasures[i] - k];
//     }
//     rs_regenerate(blocksize, &context, &conf, (uint8_t **)input, (uint8_t **)output);

//     rs_free_regenerate_context(&conf, &context);
//     return 0;
// }

// unsigned ErasureCodeTerasureReedSolomonVandermonde::get_alignment() const
// {
//   if (per_chunk_alignment) {
//     return w * LARGEST_VECTOR_WORDSIZE;
//   } else {
//     unsigned alignment = k*w*sizeof(int);
//     if ( ((w*sizeof(int))%LARGEST_VECTOR_WORDSIZE) )
//       alignment = k*w*LARGEST_VECTOR_WORDSIZE;
//     return alignment;
//   }
// }

// int ErasureCodeTerasureReedSolomonVandermonde::parse(ErasureCodeProfile &profile,
// 						     ostream *ss)
// {
//   int err = 0;
//   err |= ErasureCodeTerasure::parse(profile, ss);
//     rs_init(&conf, k+m, k, malloc, free);
//     vandermonde_matrix(&conf);
//   if (w != 8 && w != 16 && w != 32) {
//     *ss << "ReedSolomonVandermonde: w=" << w
// 	<< " must be one of {8, 16, 32} : revert to " << DEFAULT_W << std::endl;
//     profile["w"] = "8";
//     err |= to_int("w", profile, &w, DEFAULT_W, ss);
//     err = -EINVAL;
//   }
//   err |= to_bool("Terasure-per-chunk-alignment", profile,
// 		 &per_chunk_alignment, "false", ss);
//   return err;
// }

// void ErasureCodeTerasureReedSolomonVandermonde::prepare()
// {

// }


// 
// ErasureCodeTerasureReedSolomonCauchy
//
void ErasureCodeTerasureReedSolomonCauchy::Terasure_encode(char **data,
						char **coding,
						int blocksize)
{
    rs_encode(blocksize, &conf, (uint8_t **)data, NULL, (uint8_t**)coding);
}

int ErasureCodeTerasureReedSolomonCauchy::Terasure_decode(int *erasures,
					       char **data,
					       char **coding,
					       int blocksize)
{
    int erasure_num = 0;
    for(int i = 0; erasures[i] != -1; i++) erasure_num++;
    if(erasure_num > m) return -1;
    rs_regenerate_context context;
    rs_regenerate_context_init(&conf, &context, erasures, erasure_num);
    char *input[k+m];
    char *output[m];
    for(int i = 0; i < k + m; i++){
        input[i] = (i < k) ? data[i]:coding[i-k];
    }
    for(int i = 0; i < erasure_num; i++){
        if(erasures[i] < k) output[i] = data[erasures[i]];
        else output[i] = coding[erasures[i] - k];
    }
    rs_regenerate(blocksize, &context, &conf, (uint8_t **)input, (uint8_t **)output);

    rs_free_regenerate_context(&conf, &context);
    return 0;
}

int ErasureCodeTerasureReedSolomonCauchy::parse(ErasureCodeProfile &profile,
				     ostream *ss)
{
  int err = 0;
  err |= ErasureCodeTerasure::parse(profile, ss);
  return err;
}

void ErasureCodeTerasureReedSolomonCauchy::prepare()
{
  rs_init(&conf, k+m, k, malloc, free);
  cauchy_matrix(&conf);
}


//ErasureCodeTerasureCLMSR

void ErasureCodeTerasureCLMSR::Terasure_encode(char **data,
            char **coding,
            int blocksize){
  uint8_t *memery_pre_allocated[m];
  uint8_t *input[k+m];
  for(int i = 0; i < k + m; i++){
    if(i < k) input[i] = (uint8_t *)data[i];
    else{
      input[i] = NULL;
      posix_memalign((void **)&(memery_pre_allocated[i - k]), SIMD_ALIGN, sizeof(uint8_t) * blocksize);
    }
  }

  msr_encode_context context;
  msr_fill_encode_context(&context, &conf, input);

  uint8_t *buf;
  posix_memalign((void **)&buf, SIMD_ALIGN, context.encoding_buf_size * sizeof(uint8_t));

  msr_encode(blocksize, &context, &conf, buf, input, memery_pre_allocated);

  for(int i = 0; i < m; i++){
    memcpy(coding[i], memery_pre_allocated[i], blocksize);
    free(memery_pre_allocated[i]);
  }
  free(buf);
  msr_free_encode_context(&conf, &context);

  return;
}

int ErasureCodeTerasureCLMSR::regenerating(int erasure_node,
                    map<int, bufferlist> &chunks,
                    bufferlist *out,
                    int blocksize){
  dout(10)<<"Regenerating Code: regenerating"<<dendl;
    if((int)chunks.size() != k+m-1 || chunks.find(erasure_node) != chunks.end()) 
      return -EIO;
    msr_regenerate_context regenerate_context;
    msr_fill_regenerate_context(&regenerate_context,&conf,erasure_node);

    uint8_t *regenerate_buf;
    posix_memalign((void **)&regenerate_buf,SIMD_ALIGN,sizeof(uint8_t)* regenerate_context.regenerate_buf_size);

    uint8_t *input[k+m];
    for(int i = 0; i < k + m; i++){
      if(i != erasure_node){
        posix_memalign((void **)&(input[i]), SIMD_ALIGN, blocksize * sizeof(uint8_t));
        memcpy(input[i], chunks[i].c_str(), blocksize * sizeof(uint8_t));
      }
      else input[i] = NULL;
    }
    uint8_t * memory;
    posix_memalign((void **)(&memory),SIMD_ALIGN,sizeof(uint8_t) * blocksize * conf.r);
    msr_regenerate(blocksize, &regenerate_context, &conf, regenerate_buf, input, memory);

    out->append((char *)memory, blocksize * conf.r);

    for(int i = 0; i < k+m; i++){
      if(i != erasure_node) free(input[i]);
    }
    free(regenerate_buf);
    free(memory);
    msr_free_regenerate_context(&conf, &regenerate_context);
    return 0;
}

int ErasureCodeTerasureCLMSR::Terasure_decode(int *erasures,
                 char **data,
                 char **coding,
                 int blocksize){
  int erasure_num = 0;
  int erased[k+m];
  dout(10)<<"Decode not regenerating"<<dendl;
  memset(erased, 0 ,sizeof(int) * (k + m));
  for(int i = 0; erasures[i] != -1; i++){
    erasure_num++;
    erased[erasures[i]] = 1;
  }
  if(erasure_num > m) return -1;

  uint8_t *input[k+m];
  uint8_t *memery_pre_allocated[m];
  for(int i = 0; i < k + m; i++){
    if(erased[i]){
      input[i] = NULL;
      continue;
    }
    if(i < k) input[i] = (uint8_t *)data[i];
    else input[i] = (uint8_t *)coding[i - k];
  }

  for(int i = 0; i < m; i++){
    posix_memalign((void **)&(memery_pre_allocated[i]), SIMD_ALIGN, sizeof(uint8_t) * blocksize);
  }
  msr_encode_context context;
  msr_fill_encode_context(&context, &conf, input);

  uint8_t *buf;
  posix_memalign((void **)&buf, SIMD_ALIGN, context.encoding_buf_size * sizeof(uint8_t));

  msr_encode(blocksize, &context, &conf, buf, input, memery_pre_allocated);

  int num = 0;
  for(int i = 0; i < k + m; i++){
    if(erased[i]){
      if(i < k)
        memcpy(data[i], memery_pre_allocated[num], blocksize);
      else 
        memcpy(coding[i - k], memery_pre_allocated[num], blocksize);
      free(memery_pre_allocated[num]);
      num++;
    }
  }
  free(buf);
  msr_free_encode_context(&conf, &context);

  return 0;
}

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

int ErasureCodeTerasureCLMSR::minimum_to_decode(const set<int> &want_to_read,
                                   const set<int> &available_chunks,
                                   map<int, vector<pair<int, int> > > *minimum){
  dout(10)<<"available_chunks number is "<<(int)available_chunks.size()<<dendl;
  if((int)available_chunks.size() == conf.n-1){
    dout(10)<<"Regenerating Code: minimum_to_decode"<<dendl;
    vector<pair<int, int> > default_subchunks;
    int x0,y0;
    for(int i = 0; i < k + m; i++){
      if(available_chunks.find(i) == available_chunks.end()){
        x0 = i / conf.r;
        y0 = i % conf.r;
        break;
      }
    }
    for(int i = 0; i < conf.beta; i++){
      int z = (i / ipow(conf.r, conf.groups - y0 - 1) * conf.r + x0) * ipow(conf.r, conf.groups - y0 - 1) +
                i % ipow(conf.r, conf.groups - y0 - 1);
      default_subchunks.push_back(make_pair(z, 1));
    }
    for(auto &&i:available_chunks){
      minimum->insert(make_pair(i, default_subchunks));
    }
    return 0;
  }
  else return ErasureCode::minimum_to_decode(want_to_read, available_chunks, minimum);
}

int ErasureCodeTerasureCLMSR::get_subchunks_count(){
  return conf.alpha;
}

int ErasureCodeTerasureCLMSR::parse(ErasureCodeProfile &profile,
             ostream *ss)
{
  int err = 0;
  err |= ErasureCodeTerasure::parse(profile, ss);
  return err;
}

void ErasureCodeTerasureCLMSR::prepare()
{
  msr_init(&conf, k+m, k, malloc, free);
}

