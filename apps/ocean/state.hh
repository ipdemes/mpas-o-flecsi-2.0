/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "flecsi/execution.hh"
#include "mpasoflecsi/common/types.hh"
#include "mpasoflecsi/specialization/mesh.hh"
#include "mpasoflecsi/io/desc.hh"

namespace mpas { namespace ocean {

enum class test { case1, case2, case3 };
inline std::istream & operator>>(std::istream & is, test & tc)
{
  std::string token;
  is >> token;

  if (token == "1")
    tc = test::case1;
  else if (token == "2")
    tc = test::case2;
  else if (token == "3")
    tc = test::case3;
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
  case test::case2 :
    os << "case 2";
    break;
  case test::case3 :
    os << "case 3";
    break;
  }
  return os;
}


namespace inputs {

inline flecsi::program_option<std::string>
meshfile("mesh-file", "Filename for the mesh file.", 1);

inline flecsi::program_option<test>
test_case("Problem", "testcase,t", "Test case to run [1,2,3].",
          {{flecsi::option_default, test::case1}});

inline flecsi::program_option<std::size_t>
nsteps("Timesteps", "nsteps,s", "Number of timesteps to run.",
       {{flecsi::option_default, 1}});

inline flecsi::program_option<double>
delta_t("Problem", "dt,d", "Timestep size for time integrator.",
        {{flecsi::option_default, 96}});
}

enum time_versions { prev, curr, next, prov, num_timeversions };
inline const flecsi::field<double>::definition<mesh, mesh::edges>    meshScalingDel2;

inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::cells>    layerThickness[num_timeversions];
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::edges>    normalVelocity[num_timeversions];
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::cells>    hProvInv;
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::cells>    hNewInv;
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::edges>    normalThicknessFlux;

inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::edges>    layerThicknessEdge;
inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::edges>    w;

inline const flecsi::field<metensor<double>>::definition<mesh, mesh::cells>    edgeSignOnCell;

inline const flecsi::field<me2tensor<double>>::definition<mesh, mesh::edges>    advCoef;
inline const flecsi::field<me2tensor<double>>::definition<mesh, mesh::edges>    advCoef_3rd;

inline const flecsi::field<vltracer>::definition<mesh, mesh::cells>    flxIn;
inline const flecsi::field<vltracer>::definition<mesh, mesh::cells>    flxOut;
inline const flecsi::field<vltracer>::definition<mesh, mesh::edges>    highOrderFlux;
inline const flecsi::field<vltracer>::definition<mesh, mesh::edges>    lowOrderFlux;
inline const flecsi::field<vltracer>::definition<mesh, mesh::cells>    tracerMax;
inline const flecsi::field<vltracer>::definition<mesh, mesh::cells>    tracerMin;
inline const flecsi::field<vltracer>::definition<mesh, mesh::cells>    tracers[num_timeversions];
inline const flecsi::field<vltracer>::definition<mesh, mesh::cells>    tracers_tend;

inline const flecsi::field<derivtensor>::definition<mesh, mesh::edges>    derivTwo;

inline const flecsi::field<int>::definition<mesh, mesh::edges>    maxLevelEdgeTop;
inline const flecsi::field<int>::definition<mesh, mesh::edges>    nAdvCellsForEdge;

inline const flecsi::field<vltensor<int>>::definition<mesh, mesh::cells>    boundaryCell;
inline const flecsi::field<vltensor<int>>::definition<mesh, mesh::edges>    highOrderAdvMask;

inline const flecsi::field<me2tensor<int>>::definition<mesh, mesh::edges>    advCellsForEdge;

inline mesh::slot m;
inline mesh::cslot coloring;

}}
