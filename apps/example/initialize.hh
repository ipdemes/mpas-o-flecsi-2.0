/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <mpasoflecsi/specialization/control.hh>
#include <mpasoflecsi/specialization/mesh.hh>


#include <flecsi/flog.hh>

#include "options.hh"
#include "state.hh"

namespace mpas {
namespace action {

int
init() {
  flog(info) << "Initializing mesh: " << mesh_filename.value() << std::endl;

  coloring.allocate(mesh_filename.value());
  // m.allocate(coloring.get());

  return 0;
} //init_mesh

control::action<init, cp::initialize> init_action;

} // namespace action
} // namespace mpas
