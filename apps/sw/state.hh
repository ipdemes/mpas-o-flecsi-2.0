/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "flecsi/execution.hh"
#include "mpasoflecsi/common/types.hh"
#include "mpasoflecsi/specialization/mesh.hh"

namespace mpas { namespace sw {

inline flecsi::program_option<std::string>
mesh_filename("mesh-file", "Filename for the mesh file.", 1);

inline flecsi::program_option<int>
test_case("Problem", "testcase,t", "Test case to run [1,5].",
          {{flecsi::option_default, 1}},
          [](flecsi::any const & v, std::stringstream &) {
            const int value = flecsi::option_value<int>(v);
            return (value == 1 or value == 5);
          });

inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::cells> h;
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::edges> u;

inline mesh::slot m;
inline mesh::cslot coloring;

}}
