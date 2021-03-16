/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "tasks/copy_fields.hh"

#include <flecsi/flog.hh>
using namespace flecsi; 

void coupler::task::copy_fields(
  cartmesh::accessor<ro> mesh_src,
  cartmesh::accessor<ro> mesh_trg, 
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<wo, ro> va) {

  auto u = mesh_src.mdspan<cartmesh::cells>(ua);
  auto v = mesh_trg.mdspan<cartmesh::cells>(va);

  //copy field value from source to target 
  for(auto j : mesh_trg.cells<cartmesh::y_axis>()) {
    for(auto i : mesh_trg.cells<cartmesh::x_axis>()) {
          v[i][j] = u[i][j]; 
    } // for
  } // for
  
} //copy fields

