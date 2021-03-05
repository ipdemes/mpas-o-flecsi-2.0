/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "tasks/io.hh"

#include <sstream>

using namespace flecsi;

void
coupler::task::io(cartmesh::accessor<ro> m, 
  field<double>::accessor<ro, ro> ua, 
  field<double>::accessor<ro, ro> va) {

  auto u = m.mdspan<cartmesh::cells>(ua);
  auto v = m.mdspan<cartmesh::cells>(va);

  std::stringstream ss;
  if(processes() == 1) {
    ss << "solution.dat";
  }
  else {
    ss << "solution-" << process() << ".dat";
  } // if

  std::ofstream solution(ss.str(), std::ofstream::out);
  solution << " U        V"<<std::endl;
 
  for(auto j : m.cells<cartmesh::y_axis>()) {
    for(auto i : m.cells<cartmesh::x_axis>()) {
          solution << u[j][i] << " "<<v[j][i] <<std::endl; 
          //solution << u[i][j] << " "<<v[i][j] <<std::endl; 
    } // for
  } // for

  /*for(auto j : m.vertices<mesh::y_axis, mesh::logical>()) {
    const double y = m.value<mesh::y_axis>(j);
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      const double x = m.value<mesh::x_axis>(i);
      solution << x << " " << y << " " << u[j][i] << std::endl;
    } // for
  } // for
  */
} // io
