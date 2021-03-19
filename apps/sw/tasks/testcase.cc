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


void setup_case_5(mesh::accessor<ro, ro> m,
                  acc<double, ro, na> dvEdge,
                  acc<double, ro, na> latCell,
                  acc<double, rw, na> lonCell,
                  acc<double, ro, na> latVertex,
                  acc<double, ro, na> lonVertex,
                  acc<double, ro, na> latEdge,
                  acc<double, ro, na> lonEdge,
                  acc<double, wo, na> h_s,
                  acc<double, wo, na> fVertex,
                  acc<double, wo, na> fEdge,
                  acc<vlreal, wo, na> u,
                  acc<vlreal, wo, na> h,
                  acc<vltracer, wo, na> tracers)
{
  using namespace mpas_constants;
  constexpr double u0 = 20;
  constexpr double gh0 = 5960.0 * gravity;
  constexpr double hs0 = 2000;
  constexpr double alpha = 0;
  constexpr double theta_c = pi / 6.0;
  constexpr double lambda_c = 3.0 * pi / 2.0;
  constexpr double rr = pi / 9.0;

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

  // generate rotated coriolis field
  for(auto e : m.edges()) {
    fEdge(e) = 2.0 * omega *
               (-std::cos(lonEdge(e)) * std::cos(latEdge(e)) * std::sin(alpha) +
                 std::sin(latEdge(e)) * std::cos(alpha));
  }
  for(auto v : m.vertices()) {
    fVertex(v) =
      2.0 * omega *
      (-std::cos(lonVertex(v)) * std::cos(latVertex(v)) * std::sin(alpha) +
        std::sin(latVertex(v)) * std::cos(alpha));
  }

  // initialize mountain
  for(auto c : m.cells()) {
    if(lonCell(c) < 0)
      lonCell(c) = lonCell(c) + 2.0 * pi;
    auto r = std::sqrt(
      std::min(std::pow(rr, 2), std::pow((lonCell(c) - lambda_c), 2) +
                                  std::pow((latCell(c) - theta_c), 2)));
    h_s(c) = hs0 * (1.0 - r / rr);
  }

  // initialize tracer fields
  for(auto c : m.cells()) {
    auto r = std::sqrt(
      std::min(std::pow(rr, 2), std::pow((lonCell(c) - lambda_c), 2) +
                                  std::pow((latCell(c) - theta_c), 2)));
    tracers(c)[0][0] = 1.0 - r / rr;
  }
  if(maxTracers > 1) {
    for(auto c : m.cells()) {
      auto r = std::sqrt(std::min(
        std::pow(rr, 2), std::pow((lonCell(c) - lambda_c), 2) +
                           std::pow((latCell(c) - theta_c - pi / 6.0), 2)));
      tracers(c)[0][1] = 1.0 - r / rr;
    }
  }

  // initialize height field (actually, fluid thickness field)
  for(auto c : m.cells()) {
    h(c)[0] = (gh0 - (a * omega * u0 + 0.5 * std::pow(u0, 2)) *
                       std::pow((-std::cos(lonCell(c)) * std::cos(latCell(c)) *
                                    std::sin(alpha) +
                                  std::sin(latCell(c)) * std::cos(alpha)),
                         2)) /
              gravity;
    h(c)[0] = h(c)[0] - h_s(c);
  }
}


void setup_case_6(mesh::accessor<ro, ro> m,
                  acc<double, ro, na> dvEdge,
                  acc<double, ro, na> latCell,
                  acc<double, ro, na> lonCell,
                  acc<double, ro, na> latVertex,
                  acc<double, ro, na> lonVertex,
                  acc<vlreal, wo, na> u,
                  acc<vlreal, wo, na> h)
{
  using namespace mpas_constants;
  constexpr double h0 = 8000.0;
  constexpr double w = 7.848e-6;
  constexpr double K = 7.848e-6;
  constexpr double R = 4.0;

  // helper functions
  const auto aa = [](auto theta) {
    return 0.5 * w * (2.0 * omega + w) * std::pow(std::cos(theta), 2) +
           0.25 * std::pow(K, 2) * std::pow(std::cos(theta), (2 * R)) *
             ((R + 1.0) * std::pow(std::cos(theta), 2) + 2.0 * std::pow(R, 2) -
               R - 2.0 - 2.0 * std::pow(R, 2) * std::pow(std::cos(theta), -2));
  };

  const auto bb = [](auto theta) {
    return (2.0 * (omega + w) * K / ((R + 1.0) * (R + 2.0))) *
           std::pow(std::cos(theta), R) *
           ((std::pow(R, 2) + 2.0 * R + 2.0) -
             std::pow((R + 1.0) * std::cos(theta), 2));
  };

  const auto cc = [](auto theta) {
    return 0.25 * std::pow(K, 2) * std::pow(std::cos(theta), (2.0 * R)) *
           ((R + 1.0) * std::pow(std::cos(theta), 2) - R - 2.0);
  };

  // initialize wind field
  const auto psivertex = [w, K, R, &latVertex, &lonVertex](auto vtx) {
    return -a * a * w * std::sin(latVertex(vtx)) +
           a * a * K * std::pow(std::cos(latVertex(vtx)), R) *
             std::sin(latVertex(vtx)) * std::cos(R * lonVertex(vtx));
  };

  for(auto e : m.edges()) {
    auto verticesOnEdge = m.vertices(e);
    u(e)[0] = -1.0 *
              (psivertex(verticesOnEdge[1]) - psivertex(verticesOnEdge[0])) /
              dvEdge(e);
  }

  // initialize height field (actually, fluid thickness field)
  for(auto c : m.cells()) {
    h(c)[0] = (gravity * h0 + a * a * aa(latCell(c)) +
                a * a * bb(latCell(c)) * std::cos(R * lonCell(c)) +
                a * a * cc(latCell(c)) * std::cos(2.0 * R * lonCell(c))) /
              gravity;
  }
}

}}}
