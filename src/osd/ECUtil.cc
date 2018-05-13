// -*- mode:C++; tab-width:8; c-basic-offset:2; indent-tabs-mode:t -*-
#include "common/debug.h"
#include <errno.h>
#include "include/encoding.h"
#include "ECUtil.h"
using namespace std;

#define dout_context g_ceph_context
#define dout_subsys ceph_subsys_osd
#undef dout_prefix
#define dout_prefix _prefix(_dout)


static ostream& _prefix(std::ostream* _dout)
{
  return *_dout << "ECUtil: ";
}

int ECUtil::decode(
  const stripe_info_t &sinfo,
  ErasureCodeInterfaceRef &ec_impl,
  map<int, bufferlist> &to_decode,
  bufferlist *out) {
  assert(to_decode.size());

  uint64_t total_data_size = to_decode.begin()->second.length();
  assert(total_data_size % sinfo.get_chunk_size() == 0);

  assert(out);
  assert(out->length() == 0);

  for (map<int, bufferlist>::iterator i = to_decode.begin();
       i != to_decode.end();
       ++i) {
    assert(i->second.length() == total_data_size);
  }

  if (total_data_size == 0)
    return 0;

  for (uint64_t i = 0; i < total_data_size; i += sinfo.get_chunk_size()) {
    map<int, bufferlist> chunks;
    for (map<int, bufferlist>::iterator j = to_decode.begin();
	 j != to_decode.end();
	 ++j) {
      chunks[j->first].substr_of(j->second, i, sinfo.get_chunk_size());
    }
    bufferlist bl;
    int r = ec_impl->decode_concat(chunks, &bl);
    assert(bl.length() == sinfo.get_stripe_width());
    assert(r == 0);
    out->claim_append(bl);
  }
  return 0;
}

int ECUtil::decode(
  const stripe_info_t &sinfo,
  ErasureCodeInterfaceRef &ec_impl,
  map<int, bufferlist> &to_decode,
  map<int, bufferlist*> &out) {
  assert(to_decode.size());
  dout(10)<<"recovery decode"<<dendl;
  set<int> have;
  for(auto &&i:to_decode){
    assert(i.second == 0);
    have.insert(i.first);
  }

  set<int> need;
  for (map<int, bufferlist*>::iterator i = out.begin();
       i != out.end();
       ++i) {
    assert(i->second);
    assert(i->second->length() == 0);
    need.insert(i->first);
  }

  map<int, vector<pair<int, int> > > min;
  int r = ec_impl->minimum_to_decode(need, have, &min);
  assert(r == 0);

  if(min.find(*(have.begin()))->second.size() == 1){
    uint64_t total_data_size = to_decode.begin()->second.length();
    assert(total_data_size % sinfo.get_chunk_size() == 0);
    assert(total_data_size != 0);
    for(auto &&i : to_decode){
      assert(total_data_size == i.second.length());
    }

    for (uint64_t i = 0; i < total_data_size; i += sinfo.get_chunk_size()) {
      map<int, bufferlist> chunks;
      for (map<int, bufferlist>::iterator j = to_decode.begin();
     j != to_decode.end();
     ++j) {
        chunks[j->first].substr_of(j->second, i, sinfo.get_chunk_size());
      }
      map<int, bufferlist> out_bls;
      int r = ec_impl->decode(need, chunks, &out_bls);
      assert(r == 0);
      for (map<int, bufferlist*>::iterator j = out.begin();
     j != out.end();
     ++j) {
        assert(out_bls.count(j->first));
        assert(out_bls[j->first].length() == sinfo.get_chunk_size());
        j->second->claim_append(out_bls[j->first]);
      }
    }
    for (map<int, bufferlist*>::iterator i = out.begin();
         i != out.end();  
         ++i) {
      assert(i->second->length() == total_data_size);
    }
  }
  else{
    assert(need.size() == 1);
    dout(10)<<"recovery regenerating"<<dendl;
    set<int>::iterator erasure_node = need.begin();
    int subchunk_size = sinfo.get_chunk_size() / ec_impl->get_subchunks_count();
    int repair_count_size = subchunk_size * (int)min.begin()->second.size();
    uint64_t total_data_size = to_decode.begin()->second.length();
    for(uint64_t i = 0; i < total_data_size; i += repair_count_size){
      map<int, bufferlist> chunks;
      for(auto &&j : to_decode){
        chunks[j.first].substr_of(j.second, i, repair_count_size);
      }
      bufferlist out_bl;
      int r = ec_impl->regenerating(*erasure_node, chunks, &out_bl, repair_count_size);
      assert(r == 0);
      assert(out_bl.length() == sinfo.get_chunk_size());
      out[*erasure_node]->claim_append(out_bl);
    }
    assert(out[*erasure_node]->length() == total_data_size * (sinfo.get_chunk_size() / repair_count_size));
  }

  return 0;
}

int ECUtil::encode(
  const stripe_info_t &sinfo,
  ErasureCodeInterfaceRef &ec_impl,
  bufferlist &in,
  const set<int> &want,
  map<int, bufferlist> *out) {

  uint64_t logical_size = in.length();

  assert(logical_size % sinfo.get_stripe_width() == 0);
  assert(out);
  assert(out->empty());

  if (logical_size == 0)
    return 0;

  for (uint64_t i = 0; i < logical_size; i += sinfo.get_stripe_width()) {
    map<int, bufferlist> encoded;
    bufferlist buf;
    buf.substr_of(in, i, sinfo.get_stripe_width());
    int r = ec_impl->encode(want, buf, &encoded);
    assert(r == 0);
    for (map<int, bufferlist>::iterator i = encoded.begin();
	 i != encoded.end();
	 ++i) {
      assert(i->second.length() == sinfo.get_chunk_size());
      (*out)[i->first].claim_append(i->second);
    }
  }

  for (map<int, bufferlist>::iterator i = out->begin();
       i != out->end();
       ++i) {
    assert(i->second.length() % sinfo.get_chunk_size() == 0);
    assert(
      sinfo.aligned_chunk_offset_to_logical_offset(i->second.length()) ==
      logical_size);
  }
  return 0;
}

void ECUtil::HashInfo::append(uint64_t old_size,
			      map<int, bufferlist> &to_append) {
  assert(old_size == total_chunk_size);
  uint64_t size_to_append = to_append.begin()->second.length();
  if (has_chunk_hash()) {
    assert(to_append.size() == cumulative_shard_hashes.size());
    for (map<int, bufferlist>::iterator i = to_append.begin();
	 i != to_append.end();
	 ++i) {
      assert(size_to_append == i->second.length());
      assert((unsigned)i->first < cumulative_shard_hashes.size());
      uint32_t new_hash = i->second.crc32c(cumulative_shard_hashes[i->first]);
      cumulative_shard_hashes[i->first] = new_hash;
    }
  }
  total_chunk_size += size_to_append;
}

void ECUtil::HashInfo::encode(bufferlist &bl) const
{
  ENCODE_START(1, 1, bl);
  ::encode(total_chunk_size, bl);
  ::encode(cumulative_shard_hashes, bl);
  ENCODE_FINISH(bl);
}

void ECUtil::HashInfo::decode(bufferlist::iterator &bl)
{
  DECODE_START(1, bl);
  ::decode(total_chunk_size, bl);
  ::decode(cumulative_shard_hashes, bl);
  projected_total_chunk_size = total_chunk_size;
  DECODE_FINISH(bl);
}

void ECUtil::HashInfo::dump(Formatter *f) const
{
  f->dump_unsigned("total_chunk_size", total_chunk_size);
  f->open_object_section("cumulative_shard_hashes");
  for (unsigned i = 0; i != cumulative_shard_hashes.size(); ++i) {
    f->open_object_section("hash");
    f->dump_unsigned("shard", i);
    f->dump_unsigned("hash", cumulative_shard_hashes[i]);
    f->close_section();
  }
  f->close_section();
}

void ECUtil::HashInfo::generate_test_instances(list<HashInfo*>& o)
{
  o.push_back(new HashInfo(3));
  {
    bufferlist bl;
    bl.append_zero(20);
    map<int, bufferlist> buffers;
    buffers[0] = bl;
    buffers[1] = bl;
    buffers[2] = bl;
    o.back()->append(0, buffers);
    o.back()->append(20, buffers);
  }
  o.push_back(new HashInfo(4));
}

const string HINFO_KEY = "hinfo_key";

bool ECUtil::is_hinfo_key_string(const string &key)
{
  return key == HINFO_KEY;
}

const string &ECUtil::get_hinfo_key()
{
  return HINFO_KEY;
}
