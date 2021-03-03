/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "state.hh"
#include "mpasoflecsi/specialization/cartmesh.hh"

namespace coupler {
namespace task {

void io(cartmesh::accessor<ro> m, 
        field<double>::accessor<ro, ro> ua,
        field<double>::accessor<ro, ro> va);

} // namespace task
} // namespace coupler
