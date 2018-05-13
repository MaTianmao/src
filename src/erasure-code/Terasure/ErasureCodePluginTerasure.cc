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

#include "ceph_ver.h"
#include "common/debug.h"
#include "ErasureCodeTerasure.h"
#include "ErasureCodePluginTerasure.h"

#define dout_context g_ceph_context
#define dout_subsys ceph_subsys_osd
#undef dout_prefix
#define dout_prefix _prefix(_dout)

static ostream& _prefix(std::ostream* _dout)
{
  return *_dout << "ErasureCodePluginTerasure: ";
}

int ErasureCodePluginTerasure::factory(const std::string& directory,
		      ErasureCodeProfile &profile,
		      ErasureCodeInterfaceRef *erasure_code,
		      std::ostream *ss) {
    ErasureCodeTerasure *interface;
    std::string t;
    if (profile.find("technique") != profile.end())
      t = profile.find("technique")->second;
    if (t == "reed_sol_cauchy") {
      interface = new ErasureCodeTerasureReedSolomonCauchy();
    }else if(t == "cl_msr"){
      interface = new ErasureCodeTerasureCLMSR();
    }else {
      derr << "technique=" << t << " is not a valid coding technique. "
	   << " Choose one of the following: "
	   << "reed_sol_cauchy"
	   << dendl;
      return -ENOENT;
    }
    dout(20) << __func__ << ": " << profile << dendl;
    int r = interface->init(profile, ss);
    if (r) {
      delete interface;
      return r;
    }
    *erasure_code = ErasureCodeInterfaceRef(interface);
    return 0;
}

#ifndef BUILDING_FOR_EMBEDDED

const char *__erasure_code_version() { return CEPH_GIT_NICE_VER; }

int __erasure_code_init(char *plugin_name, char *directory)
{
  ErasureCodePluginRegistry &instance = ErasureCodePluginRegistry::instance();

  return instance.add(plugin_name, new ErasureCodePluginTerasure());
}


#endif