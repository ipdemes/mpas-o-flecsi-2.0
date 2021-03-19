/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "flecsi/execution.hh"
#include "mpasoflecsi/common/types.hh"
#include "mpasoflecsi/specialization/mesh.hh"
#include "mpasoflecsi/io/desc.hh"

namespace mpas { namespace sw {

enum class test { case1, case5, case6 };
inline std::istream & operator>>(std::istream & is, test & tc)
{
  std::string token;
  is >> token;

  if (token == "1")
    tc = test::case1;
  else if (token == "5")
    tc = test::case5;
  else if (token == "6")
    tc = test::case6;
  else
    is.setstate(std::ios_base::failbit);

  return is;
}
inline std::ostream & operator<<(std::ostream & os, test tc)
{
  switch (tc) {
  case test::case1 :
    os << "case 1";
    break;
  case test::case5 :
    os << "case 5";
    break;
  case test::case6 :
    os << "case 6";
    break;
  }
  return os;
}


namespace inputs {

inline flecsi::program_option<std::string>
meshfile("mesh-file", "Filename for the mesh file.", 1);

inline flecsi::program_option<test>
test_case("Problem", "testcase,t", "Test case to run [1,5,6].",
          {{flecsi::option_default, test::case1}});

inline flecsi::program_option<std::size_t>
    nsteps("Timesteps", "nsteps,s", "Number of timesteps to run.",
           {{flecsi::option_default, 10}});

using output_fields = type_list<io::desc::thickness,
                                io::desc::topography,
                                io::desc::vorticity,
                                io::desc::pv_vertex,
                                io::desc::tracers>;

}

enum time_versions { prev, curr=prev, next, provis, num_timeversions };
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::cells> h[num_timeversions];
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::edges> u[num_timeversions];
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::edges> h_edge[num_timeversions];
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::vertices> h_vertex[num_timeversions];
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::vertices> circulation;
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::vertices> vorticity;
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::cells> ke[num_timeversions];
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::edges> pv_edge[num_timeversions];
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::vertices> pv_vertex;
inline const flecsi::field<vltracer>::definition<mesh, mesh::cells> tracers[num_timeversions];
inline const flecsi::field<vlreal>::definition<mesh, mesh::cells> tend_h;
inline const flecsi::field<vlreal>::definition<mesh, mesh::edges> tend_u;
inline const flecsi::field<vltracer>::definition<mesh, mesh::cells> tend_tracers;

inline mesh::slot m;
inline mesh::cslot coloring;

}}
