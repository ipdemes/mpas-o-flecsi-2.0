/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "tasks/copy_fields.hh"

#include <flecsi/flog.hh>
using namespace flecsi; 

void coupler::task::copy_fields(cartmesh::accessor<ro> m,
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<wo, ro> va) {

  auto u = m.mdspan<cartmesh::cells>(ua);
  auto v = m.mdspan<cartmesh::cells>(va);

  for(auto j : m.cells<cartmesh::y_axis>()) {
    for(auto i : m.cells<cartmesh::x_axis>()) {
          v[j][i] = u[j][i]; 
    } // for
  } // for
  
} //copy fields

