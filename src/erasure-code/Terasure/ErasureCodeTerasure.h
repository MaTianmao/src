// -*- mode:C++; tab-width:8; c-basic-offset:2; indent-tabs-mode:t -*- 
// vim: ts=8 sw=2 smarttab
/*
 * Ceph distributed storage system
 *
 * Copyright (C) 2013, 2014 Cloudwatt <libre.licensing@cloudwatt.com>
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

#ifndef CEPH_ERASURE_CODE_TERASURE_H
#define CEPH_ERASURE_CODE_TERASURE_H

#include "erasure-code/ErasureCode.h"

extern "C"{
#include "reed_solomon.h"
#include "msr.h"
}
class ErasureCodeTerasure : public ErasureCode {
public:
  int k;
  std::string DEFAULT_K;
  int m;
  std::string DEFAULT_M;
  const char *technique;
  std::string rule_root;
  std::string rule_failure_domain;

  explicit ErasureCodeTerasure(const char *_technique) :
    k(0),
    DEFAULT_K("2"),
    m(0),
    DEFAULT_M("1"),
    technique(_technique)
  {}

  ~ErasureCodeTerasure() override {}
  
  unsigned int get_chunk_count() const override {
    return k + m;
  }

  unsigned int get_data_chunk_count() const override {
    return k;
  }

  unsigned int get_chunk_size(unsigned int object_size) const override;

  int encode_chunks(const std::set<int> &want_to_encode,
			    std::map<int, bufferlist> *encoded) override;

  int decode_chunks(const std::set<int> &want_to_read,
			    const std::map<int, bufferlist> &chunks,
			    std::map<int, bufferlist> *decoded) override;

  int init(ErasureCodeProfile &profile, std::ostream *ss) override;

  virtual void Terasure_encode(char **data,
                               char **coding,
                               int blocksize) = 0;
  virtual int Terasure_decode(int *erasures,
                               char **data,
                               char **coding,
                               int blocksize) = 0;

  virtual int get_subchunks_count(){
    return 1;
  }

  virtual void prepare() = 0;
  static bool is_prime(int value);
protected:
  virtual int parse(ErasureCodeProfile &profile, std::ostream *ss);
};

class ErasureCodeTerasureReedSolomonCauchy : public ErasureCodeTerasure {
public:
  rs_conf conf;

  ErasureCodeTerasureReedSolomonCauchy() :
    ErasureCodeTerasure("reed_sol_cauchy")
  {
    DEFAULT_K = "8";
    DEFAULT_M = "4";
  }
  ~ErasureCodeTerasureReedSolomonCauchy() override {
    rs_free_conf(&conf);
  }

  void Terasure_encode(char **data,
                               char **coding,
                               int blocksize) override;
  int Terasure_decode(int *erasures,
                               char **data,
                               char **coding,
                               int blocksize) override;
  void prepare() override;
private:
  int parse(ErasureCodeProfile &profile, std::ostream *ss) override;
};

class ErasureCodeTerasureCLMSR : public ErasureCodeTerasure {
public : 
  msr_conf conf;

  ErasureCodeTerasureCLMSR() :
    ErasureCodeTerasure("cl_msr")
  {
    DEFAULT_K = "10";
    DEFAULT_M = "4";
  }

  ~ErasureCodeTerasureCLMSR() override {
    msr_free_conf(&conf);
  }

  void Terasure_encode(char **data,
                               char **coding,
                               int blocksize) override;
  int Terasure_decode(int *erasures,
                               char **data,
                               char **coding,
                               int blocksize) override;

  int minimum_to_decode(const std::set<int> &want_to_read,
                                   const std::set<int> &available_chunks,
                                   std::map<int, std::vector<std::pair<int, int> > > *minimum) override;
  virtual int get_subchunks_count();


  int regenerating(int erasure_node,
                    map<int, bufferlist> &chunks,
                    bufferlist *out,
                    int blocksize) override;

  void prepare() override;
private:
  int parse(ErasureCodeProfile &profile, std::ostream *ss) override;

};


#endif
