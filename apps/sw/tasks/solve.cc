#include "../state.hh"
#include "solve.hh"


namespace mpas { namespace sw { namespace task {

using namespace flecsi;

void compute_solve_diagnostics(
    mesh::accessor<ro, ro> m,
    field<double>::accessor<ro, ro> dvEdge,
    field<double>::accessor<ro, ro> dcEdge,
    field<double>::accessor<ro, ro> areaCell,
    field<double>::accessor<ro, ro> areaTriangle,
    field<vdtensor<int>>::accessor<ro, ro> kiteAreasOnVertex,
    field<double>::accessor<ro, ro> fVertex,
    field<vltensor<double>>::accessor<ro, ro> h,
    field<vltensor<double>>::accessor<ro, ro> u,
    field<vltensor<double>>::accessor<wo, na> h_edge,
    field<vltensor<double>>::accessor<wo, na> h_vertex,
    field<vltensor<double>>::accessor<wo, na> circulation,
    field<vltensor<double>>::accessor<wo, na> vorticity,
    field<vltensor<double>>::accessor<wo, na> ke,
    field<vltensor<double>>::accessor<wo, na> pv_edge,
    field<vltensor<double>>::accessor<wo, na> pv_vertex)
{
  // compute height on cell edges at velocity locations
  auto nCells = m.cells().size();
  for(auto e : m.edges()) {
    auto cell1 = m.cells(e)[0];
    auto cell2 = m.cells(e)[1];
    // if((cell1 <= nCells - 1) and (cell2 <= nCells - 1)) {
      for(std::size_t k = 0; k < nVertLevels; k++) {
        h_edge(e)[k] = 0.5 * (h(cell1)[k] + h(cell2)[k]);
      }
    // }
  }

  // compute circulation and relative vorticity at each vertex
  for(auto v : m.vertices()) {
    for(std::size_t k = 0; k < nVertLevels; k++) {
      circulation(v)[k] = 0.0;
    }
  }
  for(auto e : m.edges()) {
    auto v1 = m.vertices(e)[0];
    auto v2 = m.vertices(e)[1];
    for(std::size_t k = 0; k < nVertLevels; k++) {
      // TODO: check sign
      circulation(v1)[k] = circulation(v1)[k] - dcEdge(e) * u(e)[k];
      circulation(v2)[k] = circulation(v2)[k] + dcEdge(e) * u(e)[k];
    }
  }
  for(auto v : m.vertices()) {
    for(std::size_t k = 0; k < nVertLevels; k++) {
      vorticity(v)[k] = circulation(v)[k] / areaTriangle(v);
    }
  }

  // compute kinetic energy in each cell
  for(auto c : m.cells()) {
    for(std::size_t k = 0; k < nVertLevels; k++) {
      ke(c)[k] = 0.0;
    }

    for(auto e : m.edges(c)) {
      for(std::size_t k = 0; k < nVertLevels; k++) {
        ke(c)[k] += 0.25 * dcEdge(e) * dvEdge(e) * std::pow(u(e)[k], 2);
      }
    }

    for(std::size_t k = 0; k < nVertLevels; k++) {
      ke(c)[k] /= areaCell(c);
    }
  }

  // compute height at vertices, pv at vertices, and average pv to edge
  // locations
  for(auto v : m.vertices()) {
    for(std::size_t k = 0; k < nVertLevels; k++) {
      h_vertex(v)[k] = 0.0;
      for(std::size_t i = 0; i < vertexDegree; i++) {
        // h_vertex(v)[k] += h(cellsOnVertex(v)[i])[k] * kiteAreasOnVertex(v)[i];
        h_vertex(v)[k] += h(m.cells(v)[i])[k] * kiteAreasOnVertex(v)[i];
      }
      h_vertex(v)[k] /= areaTriangle(v);

      pv_vertex(v)[k] = (fVertex(v) + vorticity(v)[k]) / h_vertex(v)[k];
    }
  }

  // compute pv at the edges
  for(auto e : m.edges()) {
    for(std::size_t k = 0; k < nVertLevels; k++) {
      pv_edge(e)[k] = 0.0;
    }
  }
  for(auto v : m.vertices()) {
    for(auto e : m.edges(v)) {
      for(std::size_t k = 0; k < nVertLevels; k++) {
        pv_edge(e)[k] += 0.5 * pv_vertex(v)[k];
      }
    }
  }
}


void init_provis(mesh::accessor<ro, ro> m,
                 acc<vlreal, ro, na> u_old,
                 acc<vlreal, ro, na> h_old,
                 acc<vltracer, ro, na> tracers_old,
                 acc<vlreal, ro, na> h_edge_old,
                 acc<vlreal, ro, na> ke_old,
                 acc<vlreal, ro, na> pv_edge_old,
                 acc<vlreal, wo, na> u_provis,
                 acc<vlreal, wo, na> h_provis,
                 acc<vltracer, wo, na> tracers_provis,
                 acc<vlreal, wo, na> h_edge_provis,
                 acc<vlreal, wo, na> ke_provis,
                 acc<vlreal, wo, na> pv_edge_provis)
{
  for(auto c : m.cells()) {
    for(std::size_t j = 0; j < nVertLevels; j++) {
      h_provis(c)[j] = h_old(c)[j];
      ke_provis(c)[j] = ke_old(c)[j];
      for(std::size_t k = 0; k < maxTracers; k++) {
        tracers_provis(c)[j][k] = tracers_old(c)[j][k];
      }
    }
  }

  for(auto e : m.edges()) {
    for(std::size_t j = 0; j < nVertLevels; j++) {
      u_provis(e)[j] = u_old(e)[j];
      h_edge_provis(e)[j] = h_edge_old(e)[j];
      pv_edge_provis(e)[j] = pv_edge_old(e)[j];
    }
  }
}


void setup_timestep(mesh::accessor<ro, ro> m,
                    acc<vlreal, ro, na> u_old,
                    acc<vlreal, ro, na> h_old,
                    acc<vltracer, ro, na> tracers_old,
                    acc<vlreal, wo, na> u_new,
                    acc<vlreal, wo, na> h_new,
                    acc<vltracer, wo, na> tracers_new)
{
  for(auto e : m.edges()) {
    for(std::size_t j = 0; j < nVertLevels; j++) {
      u_new(e)[j] = u_old(e)[j];
    }
  }

  for(auto c : m.cells()) {
    for(std::size_t j = 0; j < nVertLevels; j++) {
      h_new(c)[j] = h_old(c)[j];
      for(std::size_t k = 0; k < maxTracers; k++) {
        tracers_new(c)[j][k] = tracers_old(c)[j][k] * h_old(c)[j];
      }
    }
  }
}


void compute_tend(mesh::accessor<ro, ro> m,
                  acc<double, ro, ro> areaCell,
                  acc<double, ro, ro> dvEdge,
                  acc<double, ro, ro> dcEdge,
                  acc<me2tensor<double>, ro, ro> weightsOnEdge,
                  acc<vlreal, ro, ro> h_edge,
                  acc<vlreal, ro, ro> h,
                  acc<double, ro, ro> h_s,
                  acc<vlreal, ro, ro> u,
                  acc<vlreal, ro, ro> ke,
                  acc<vlreal, ro, ro> pv_edge,
                  acc<vlreal, wo, na> tend_h,
                  acc<vlreal, wo, na> tend_u)
{
  using namespace mpas_constants;

  // Compute height tendency for each cell
  for(auto c : m.cells()) {
    for(auto k = 0; k < nVertLevels; k++) {
      tend_h(c)[k] = 0;
    }
  }
  for(auto e : m.edges()) {
    auto cell1 = m.cells(e)[0];
    auto cell2 = m.cells(e)[1];

    // TODO: double check sign of edge
    for(std::size_t k = 0; k < nVertLevels; k++) {
      auto flux = u(e)[k] * dvEdge(e) * h_edge(e)[k];
      tend_h(cell1)[k] -= flux;
      tend_h(cell2)[k] += flux;
    }
  }
  for(auto c : m.cells()) {
    for(std::size_t k = 0; k < nVertLevels; k++) {
      tend_h(c)[k] /= areaCell(c);
    }
  }

  // Compute u (normal) velocity tendency for each edge
  for(auto e : m.edges()) {
    for(std::size_t k = 0; k < nVertLevels; k++) {
      tend_u(e)[k] = 0.0;
    }
  }

  for(auto e : m.edges()) {
    auto cell1 = m.cells(e)[0];
    auto cell2 = m.cells(e)[1];

    for(std::size_t k = 0; k < nVertLevels; k++) {
      real_t q = 0;
      std::size_t j = 0;
      for(auto eoe : m.edges(e)) {
        real_t workpv = 0.5 * (pv_edge(e)[k] + pv_edge(eoe)[k]);
        q = q + weightsOnEdge(e)[j] * u(eoe)[k] * workpv * h_edge(eoe)[k];
        ++j;
      }

      tend_u(e)[k] =
        q - (ke(cell2)[k] - ke(cell1)[k] +
              gravity * (h(cell2)[k] + h_s(cell2) - h(cell1)[k] - h_s(cell1))) /
              dcEdge(e);
    }
  }
}


void compute_scalar_tend(mesh::accessor<ro, ro> m,
                         acc<double, ro, ro> areaCell,
                         acc<double, ro, ro> dvEdge,
                         acc<vlreal, ro, ro> h_edge,
                         acc<vlreal, ro, ro> u,
                         acc<vltracer, ro, ro> tracers,
                         acc<vltracer, wo, na> tracer_tend)
{
  auto nTracers = maxTracers;
  for(auto c : m.cells()) {
    for(std::size_t k = 0; k < nVertLevels; k++) {
      for(std::size_t i = 0; i < nTracers; i++) {
        tracer_tend(c)[k][i] = 0.0;
      }
    }
  }

  // auto nCells = m.cells().size();
  for(auto e : m.edges()) {
    auto cell1 = m.cells(e)[0];
    auto cell2 = m.cells(e)[1];
    // if((cell1.id() <= nCells - 1) and (cell2.id() <= nCells - 1)) {
      for(std::size_t k = 0; k < nVertLevels; k++) {
        for(std::size_t i = 0; i < nTracers; i++) {
          auto tracer_edge =
            0.5 * (tracers(cell1)[k][i] + tracers(cell2)[k][i]);
          auto flux = u(e)[k] * dvEdge(e) * h_edge(e)[k] * tracer_edge;
          tracer_tend(cell1)[k][i] -= flux / areaCell(cell1);
          tracer_tend(cell2)[k][i] += flux / areaCell(cell2);
        }
      }
    // }
  }
}


void compute_substep(mesh::accessor<ro, ro> m,
                     acc<vlreal, wo, na> u_provis,
                     acc<vlreal, ro, na> u_old,
                     acc<vlreal, ro, na> u_tend,
                     acc<vlreal, wo, na> h_provis,
                     acc<vlreal, ro, na> h_old,
                     acc<vlreal, ro, na> h_tend,
                     acc<vltracer, wo, na> tracers_provis,
                     acc<vltracer, ro, na> tracers_old,
                     acc<vltracer, ro, na> tracers_tend,
                     double rkweight)
{
  for(auto e : m.edges()) {
    for(std::size_t j = 0; j < nVertLevels; j++) {
      u_provis(e)[j] = u_old(e)[j] + rkweight * u_tend(e)[j];
    }
  }

  for(auto c : m.cells()) {
    for(std::size_t j = 0; j < nVertLevels; j++) {
      h_provis(c)[j] = h_old(c)[j] + rkweight * h_tend(c)[j];
    }
  }

  for(auto c : m.cells()) {
    for(std::size_t j = 0; j < nVertLevels; j++) {
      for(std::size_t k = 0; k < maxTracers; k++) {
        tracers_provis(c)[j][k] = (h_old(c)[j] * tracers_old(c)[j][k] +
                                    rkweight * tracers_tend(c)[j][k]) /
                                  h_provis(c)[j];
      }
    }
  }

  if(inputs::test_case.value() == test::case1) {
    for(auto e : m.edges()) {
      for(std::size_t j = 0; j < nVertLevels; j++) {
        u_provis(e)[j] = u_old(e)[j];
      }
    }
  }
}


void accumulate_update(mesh::accessor<ro, ro> m,
                       acc<vlreal, rw, na> u_new,
                       acc<vlreal, ro, na> u_tend,
                       acc<vlreal, rw, na> h_new,
                       acc<vlreal, ro, na> h_tend,
                       acc<vltracer, rw, na> tracers_new,
                       acc<vltracer, ro, na> tracers_tend,
                       double rkweight)
{
  for(auto e : m.edges()) {
    for(std::size_t j = 0; j < nVertLevels; j++) {
      u_new(e)[j] += rkweight * u_tend(e)[j];
    }
  }

  for(auto c : m.cells()) {
    for(std::size_t j = 0; j < nVertLevels; j++) {
      h_new(c)[j] += rkweight * h_tend(c)[j];

      for(std::size_t k = 0; k < maxTracers; k++) {
        tracers_new(c)[j][k] += rkweight * tracers_tend(c)[j][k];
      }
    }
  }
}


void finalize_timestep(mesh::accessor<ro, ro> m,
                       acc<vlreal, ro, na> u_old,
                       acc<vlreal, wo, na> u_new,
                       acc<vlreal, ro, na> h_new,
                       acc<vltracer, rw, na> tracers_new)
{
  // decouple new scalar fields
  for(auto c : m.cells()) {
    for(std::size_t j = 0; j < nVertLevels; j++) {
      for(std::size_t k = 0; k < maxTracers; k++) {
        tracers_new(c)[j][k] /= h_new(c)[j];
      }
    }
  }

  // for case 1, wind field should be fixed
  if(inputs::test_case == test::case1) {
    for(auto e : m.edges()) {
      for(std::size_t j = 0; j < nVertLevels; j++) {
        u_new(e)[j] = u_old(e)[j];
      }
    }
  }
}


void timeshift(mesh::accessor<ro, ro> m,
               acc<vlreal, wo, na> u_old,
               acc<vlreal, wo, na> h_old,
               acc<vltracer, wo, na> tracers_old,
               acc<vlreal, wo, na> h_edge_old,
               acc<vlreal, wo, na> ke_old,
               acc<vlreal, wo, na> pv_edge_old,
               acc<vlreal, ro, na> u_new,
               acc<vlreal, ro, na> h_new,
               acc<vltracer, ro, na> tracers_new,
               acc<vlreal, ro, na> h_edge_new,
               acc<vlreal, ro, na> ke_new,
               acc<vlreal, ro, na> pv_edge_new)
{
  for(auto c : m.cells()) {
    for(std::size_t j = 0; j < nVertLevels; j++) {
      h_old(c)[j] = h_new(c)[j];
      ke_old(c)[j] = ke_new(c)[j];
      for(std::size_t k = 0; k < maxTracers; k++) {
        tracers_old(c)[j][k] = tracers_new(c)[j][k];
      }
    }
  }

  for(auto e : m.edges()) {
    for(std::size_t j = 0; j < nVertLevels; j++) {
      u_old(e)[j] = u_new(e)[j];
      h_edge_old(e)[j] = h_edge_new(e)[j];
      pv_edge_old(e)[j] = pv_edge_new(e)[j];
    }
  }
}

}}}
