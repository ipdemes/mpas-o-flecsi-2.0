/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "tasks/io.hh"

#include <flecsi/flog.hh>
#include <sstream>

using namespace flecsi;

void
coupler::task::io(
  cartmesh::accessor<ro> mesh_src, 
  cartmesh::accessor<ro> mesh_trg, 
  field<double>::accessor<ro, ro> ua, 
  field<double>::accessor<ro, ro> va) {

  auto u = mesh_src.mdspan<cartmesh::cells>(ua);
  auto v = mesh_trg.mdspan<cartmesh::cells>(va);

  std::stringstream ss;
  if(processes() == 1) {
    ss << "remap.dat";
  }
  else {
    ss << "remap-" << process() << ".dat";
  } // if

  std::ofstream solution(ss.str(), std::ofstream::out);
  solution << " SOURCE MESH AND FIELD "<<std::endl;
  solution << "    x        y        U"<<std::endl;
 
  for(auto j : mesh_src.cells<cartmesh::y_axis>()) {
    const double y = mesh_src.value<cartmesh::y_axis>(j);
    for(auto i : mesh_src.cells<cartmesh::x_axis>()) {
      const double x = mesh_src.value<cartmesh::x_axis>(i);
          solution << x <<" "<< y <<" "<< u[i][j] <<std::endl; 
    } // for
  } // for

  solution << "\n\n TARGET MESH AND FIELD "<<std::endl;
  solution << "    x        y        V"<<std::endl;
 
  for(auto j : mesh_trg.cells<cartmesh::y_axis>()) {
    const double y = mesh_trg.value<cartmesh::y_axis>(j);
    for(auto i : mesh_trg.cells<cartmesh::x_axis>()) {
      const double x = mesh_trg.value<cartmesh::x_axis>(i);
          solution << x <<" "<< y <<" "<< v[i][j] <<std::endl; 
    } // for
  } // for

} // io
