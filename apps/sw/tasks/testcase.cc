/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "mpasoflecsi/common/constants.hh"
#include "mpasoflecsi/eqns/metrics.hh"

#include "../state.hh"
#include "testcase.hh"

namespace mpas { namespace sw { namespace task {

using namespace flecsi;

void init_extra_fields(mesh::accessor<ro, ro> m,
                       acc<vltracer, wo, wo> tracers)
{
  for (auto c : m.cells()) {
    for (std::size_t j{0}; j < nVertLevels; j++) {
      for (std::size_t k{0}; k < maxTracers; k++) {
        tracers(c)[j][k] = 0.0;
      }
    }
  }
}

void setup_case_1(mesh::accessor<ro, ro> m,
                  field<double>::accessor<ro, ro> dvEdge,
                  field<double>::accessor<ro, ro> latCell,
                  field<double>::accessor<ro, ro> lonCell,
                  field<double>::accessor<ro, ro> latVertex,
                  field<double>::accessor<ro, ro> lonVertex,
                  field<vltensor<double>>::accessor<wo, na> u,
                  field<vltensor<double>>::accessor<wo, na> h)
{
  using namespace mpas_constants;
  constexpr double u0 = 2.0 * pi * a / (12.0 * 86400.0);
  constexpr double alpha = 0.0;//pi / 4.0;
  constexpr double h0 = 1000.0;
  constexpr double theta_c = 0.0;
  constexpr double lambda_c = 3.0 * pi / 2.0;

  // initialize wind field
  const auto psivertex = [u0, alpha, &latVertex, &lonVertex](auto vtx) {
    return -a * u0 *
           (std::sin(latVertex(vtx)) * std::cos(alpha) -
             std::cos(lonVertex(vtx)) * std::cos(latVertex(vtx)) *
               std::sin(alpha));
  };

  for(auto e : m.edges()) {
    auto verticesOnEdge = m.vertices(e);
    u(e)[0] = -1.0 *
              (psivertex(verticesOnEdge[1]) - psivertex(verticesOnEdge[0])) /
              dvEdge(e);
  }

  // initialize cosine bell at (theta_c, lambda_c)
  for(auto c : m.cells()) {
    auto r =
      eqns::sphere_distance(theta_c, lambda_c, latCell(c), lonCell(c), a);
    if(r < a / 3.0)
      h(c)[0] = (h0 / 2.0) * (1.0 + cos(pi * r * 3.0 / a));
    else
      h(c)[0] = h0 / 2.0;
  }
}

}}}
