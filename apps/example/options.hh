/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "flecsi/execution.hh"

flecsi::program_option<std::string>
mesh_filename("mesh-file", "Filename for the mesh file.", 1);
