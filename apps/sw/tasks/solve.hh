/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "mpasoflecsi/specialization/mesh.hh"

namespace mpas { namespace sw { namespace task {

void compute_solve_diagnostics(mesh::accessor<flecsi::ro, flecsi::ro> m,
                               flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> dvEdge,
                               flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> dcEdge,
                               flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> areaCell,
                               flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> areaTriangle,
                               flecsi::field<vdtensor<double>>::accessor<flecsi::ro, flecsi::ro> kiteAreasOnVertex,
                               flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> fVertex,
                               flecsi::field<vltensor<double>>::accessor<flecsi::ro, flecsi::ro> h,
                               flecsi::field<vltensor<double>>::accessor<flecsi::ro, flecsi::ro> u,
                               flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> h_edge,
                               flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> h_vertex,
                               flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> circulation,
                               flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> vorticity,
                               flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> ke,
                               flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> pv_edge,
                               flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> pv_vertex);

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
                 acc<vlreal, wo, na> pv_edge_provis);


/**
   Initialize state for a timestep

   - Initialize new state with old state
   - Initialize first RK state
   - Couple new tracers with h
*/
void setup_timestep(mesh::accessor<ro, ro> m,
                    acc<vlreal, ro, na> u_old,
                    acc<vlreal, ro, na> h_old,
                    acc<vltracer, ro, na> tracers_old,
                    acc<vlreal, wo, na> u_new,
                    acc<vlreal, wo, na> h_new,
                    acc<vltracer, wo, na> tracers_new);


/**
   Compute height and normal wind tendencies.
*/
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
                  acc<vlreal, wo, na> tend_u);


void compute_scalar_tend(mesh::accessor<ro, ro> m,
                         acc<double, ro, ro> areaCell,
                         acc<double, ro, ro> dvEdge,
                         acc<vlreal, ro, ro> h_edge,
                         acc<vlreal, ro, ro> u,
                         acc<vltracer, ro, ro> tracers,
                         acc<vltracer, wo, na> tracer_tend);


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
                     double rkweight);


void accumulate_update(mesh::accessor<ro, ro> m,
                       acc<vlreal, rw, na> u_new,
                       acc<vlreal, ro, na> u_tend,
                       acc<vlreal, rw, na> h_new,
                       acc<vlreal, ro, na> h_tend,
                       acc<vltracer, rw, na> tracers_new,
                       acc<vltracer, ro, na> tracers_tend,
                       double rkweight);


void finalize_timestep(mesh::accessor<ro, ro> m,
                       acc<vlreal, ro, na> u_old,
                       acc<vlreal, rw, na> u_new,
                       acc<vlreal, ro, na> h_new,
                       acc<vltracer, rw, na> tracers_new);


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
               acc<vlreal, ro, na> pv_edge_new);

double compute_error(mesh::accessor<ro, ro> m,
                     acc<vlreal, ro, na> h,
                     acc<double, ro, na> latCell,
                     acc<double, ro, na> lonCell,
                     acc<double, ro, na> areaCell,
                     double elapsed);

}}}
