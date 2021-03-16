/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "initialize.hh"
#include "options.hh"
#include "state.hh"

#include <flecsi/flog.hh>

using namespace flecsi;

int coupler::action::init_mesh() {
  flog(info) << "Initializing " << x_extents.value() << "x" << y_extents.value()
             << " mesh" << std::endl;
  flecsi::log::flush();

  std::vector<std::size_t> axis_extents{x_extents.value(), y_extents.value()};

  // Distribute the number of processes over the axis colors.
  auto axis_colors = cartmesh::distribute(flecsi::processes(), axis_extents);

  //source mesh
  coloring_src.allocate(axis_colors, axis_extents);
  cartmesh::grect src_bnds;
  src_bnds[0][0] = 0.0;
  src_bnds[0][1] = 1.0;
  src_bnds[1] = src_bnds[0];
  mesh_src.allocate(coloring_src.get(), src_bnds);

  //target mesh
  coloring_trg.allocate(axis_colors, axis_extents);
  cartmesh::grect trg_bnds;
  trg_bnds[0][0] = 0.0;
  trg_bnds[0][1] = 1.0;
  trg_bnds[1] = trg_bnds[0];
  mesh_trg.allocate(coloring_trg.get(), trg_bnds);

  return 0;
} // init_mesh


