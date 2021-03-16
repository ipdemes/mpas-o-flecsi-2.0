/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "tasks/init_fields.hh"

using namespace flecsi; 

void coupler::task::init_fields(
  cartmesh::accessor<ro> mesh_src,
  cartmesh::accessor<ro> mesh_trg, 
  field<double>::accessor<wo, wo> ua,
  field<double>::accessor<wo, wo> va) {

  auto u = mesh_src.mdspan<cartmesh::cells>(ua);
  auto v = mesh_trg.mdspan<cartmesh::cells>(va);

  const auto dsqr = pow(mesh_src.delta(), 2);

  //init u on src
  for(auto j : mesh_src.cells<cartmesh::y_axis>()) {
    for(auto i : mesh_src.cells<cartmesh::x_axis>()) {
          u[i][j] = dsqr; 
    } // for
  } // for

  //init v on trg
  for(auto j : mesh_trg.cells<cartmesh::y_axis>()) {
    for(auto i : mesh_trg.cells<cartmesh::x_axis>()) {
          v[i][j] = 0.0; 
    } // for
  } // for
  
} //init fields

