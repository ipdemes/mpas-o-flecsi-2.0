/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "tasks/init_fields.hh"

using namespace flecsi; 

void coupler::task::init_fields(cartmesh::accessor<ro> m,
  field<double>::accessor<wo, wo> ua,
  field<double>::accessor<wo, wo> va) {

  auto u = m.mdspan<cartmesh::cells>(ua);
  auto v = m.mdspan<cartmesh::cells>(va);
  const auto dsqr = pow(m.delta(), 2);

  //init u
  for(auto j : m.cells<cartmesh::y_axis>()) {
    for(auto i : m.cells<cartmesh::x_axis>()) {
          u[j][i] = dsqr; 
    } // for
  } // for

  //init v
  for(auto j : m.cells<cartmesh::y_axis>()) {
    for(auto i : m.cells<cartmesh::x_axis>()) {
          v[j][i] = 0.0; 
    } // for
  } // for
  
} //init fields

