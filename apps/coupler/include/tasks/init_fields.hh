/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "state.hh"
#include "mpasoflecsi/specialization/cartmesh.hh"

namespace coupler {
namespace task {

void init_fields(cartmesh::accessor<ro> m,
  field<double>::accessor<wo, wo> ua,
  field<double>::accessor<wo, wo> va); 

} // namespace task
} // namespace coupler
