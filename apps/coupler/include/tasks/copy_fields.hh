/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "specialization/cartmesh.hh"

namespace coupler {
namespace task {

void copy_fields(cartmesh::accessor<ro> m,
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<wo, ro> va); 

} // namespace task
} // namespace coupler
