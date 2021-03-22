/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "state.hh"
#include "mpasoflecsi/specialization/cartmesh.hh"

namespace coupler {
namespace task {

void copy_fields(cartmesh::accessor<ro> mesh_src,
	         cartmesh::accessor<ro> mesh_trg,
  		 field<double>::accessor<ro, ro> ua,
  		 field<double>::accessor<wo, ro> va); 

} // namespace task
} // namespace coupler