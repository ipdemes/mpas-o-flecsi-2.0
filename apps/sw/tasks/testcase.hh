/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "mpasoflecsi/specialization/mesh.hh"

namespace mpas { namespace sw { namespace task {

/**
 * Initialize fields that may be left uninitialized by test cases.
 */
void init_extra_fields(mesh::accessor<flecsi::ro, flecsi::ro> m,
                       acc<vltracer, wo, wo> tracers);

void setup_case_1(mesh::accessor<flecsi::ro, flecsi::ro> m,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> dvEdge,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> latCell,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> lonCell,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> latVertex,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> lonVertex,
                  flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> u,
                  flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> h);

}}}
