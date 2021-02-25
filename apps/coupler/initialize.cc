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

  coloring.allocate(axis_colors, axis_extents);

  cartmesh::grect geometry;
  geometry[0][0] = 0.0;
  geometry[0][1] = 1.0;
  geometry[1] = geometry[0];

  m.allocate(coloring.get(), geometry);
  return 0;
} // init_mesh


