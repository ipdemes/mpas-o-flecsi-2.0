/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "mpasoflecsi/specialization/mesh.hh"

namespace mpas { namespace ocean { namespace task {

double rotation_velocity(double latE, double angleE);

double initial_gaussian_dist(double xC);

double initial_slotCyl_dist(double xC, double yC, double zC, double latC, double lonC);

double spherical_harmonic_dist(double latC, double lonC, double coef, double elapsed);

void init_extra_fields(mesh::accessor<ro, ro> m,
                       flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> w,
                       flecsi::field<vltensor<int>>::accessor<flecsi::wo, flecsi::na> bc,
                       flecsi::field<double>::accessor<flecsi::wo, flecsi::na> bld,
                       flecsi::field<double>::accessor<flecsi::wo, flecsi::na> ss,
                       flecsi::field<int>::accessor<flecsi::wo, flecsi::na> mlc);

void setup_case_1(mesh::accessor<flecsi::ro, flecsi::ro> m,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> angleEdge,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> latEdge,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> xCell,
                  flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> u,
                  flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> h,
                  flecsi::field<vltracer>::accessor<flecsi::wo, flecsi::na> tracers);

void setup_case_2(mesh::accessor<ro, ro> m,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> angleEdge,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> latEdge,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> latCell,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> lonCell,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> xCell,  
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> yCell,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> zCell,
                  flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> u,
                  flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> h,
                  flecsi::field<vltracer>::accessor<flecsi::wo, flecsi::na> tracers);

void setup_case_3(mesh::accessor<flecsi::ro, flecsi::ro> m,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> angleEdge,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> latEdge,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> latCell,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> lonCell,
                  flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> u,
                  flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> h,
                  flecsi::field<vltracer>::accessor<flecsi::wo, flecsi::na> tracers,
                  double coef);

void calc_error_case_1(mesh::accessor<flecsi::ro, flecsi::ro> m,
                      flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> latCell,
                      flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> lonCell,
                      flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> xCell,
                      flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> yCell,
                      flecsi::field<vltracer>::accessor<flecsi::ro, flecsi::ro> tracers,
                      double timeElapsed);

void calc_error_case_2(mesh::accessor<flecsi::ro, flecsi::ro> m,
                       flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> latCell,
                       flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> lonCell,
                       flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> xCell,  
                       flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> yCell,
                       flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> zCell,
                       flecsi::field<vltracer>::accessor<flecsi::ro, flecsi::ro> tracers,
                       double timeElapsed);

void calc_error_case_3(mesh::accessor<flecsi::ro, flecsi::ro> m,
                       flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> latCell,
                       flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> lonCell,
                       flecsi::field<vltracer>::accessor<flecsi::ro, flecsi::ro> tracers,
                       double hmix_coef,
                       double timeElapsed);

}}}
