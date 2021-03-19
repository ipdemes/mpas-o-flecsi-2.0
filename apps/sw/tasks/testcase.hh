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

////////////////////////////////////////////////////////////////////////////////
//! \brief Setup shallow water test case 1: Advection of Cosine Bell over the
//! Pole
//!
//! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical
//! Approximations to the Shallow Water Equations in Spherical Geometry"
//! J. of Comp. Phys., 102, pp. 211--224
////////////////////////////////////////////////////////////////////////////////
void setup_case_1(mesh::accessor<flecsi::ro, flecsi::ro> m,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> dvEdge,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> latCell,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> lonCell,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> latVertex,
                  flecsi::field<double>::accessor<flecsi::ro, flecsi::ro> lonVertex,
                  flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> u,
                  flecsi::field<vltensor<double>>::accessor<flecsi::wo, flecsi::na> h);

////////////////////////////////////////////////////////////////////////////////
//! \brief Setup shallow water test case 5: Zonal Flow over an Isolated Mountain
//!
//! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical
//! Approximations to the Shallow Water Equations in Spherical Geometry"
//! J. of Comp. Phys., 102, pp. 211--224
////////////////////////////////////////////////////////////////////////////////
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
                  acc<vltracer, wo, na> tracers);


////////////////////////////////////////////////////////////////////////////////
//! \brief Setup shallow water test case 6: Rossby-Haurwitz Wave
//!
//! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical
//! Approximations to the Shallow Water Equations in Spherical Geometry"
//! J. of Comp. Phys., 102, pp. 211--224
////////////////////////////////////////////////////////////////////////////////
void setup_case_6(mesh::accessor<ro, ro> m,
                  acc<double, ro, na> dvEdge,
                  acc<double, ro, na> latCell,
                  acc<double, ro, na> lonCell,
                  acc<double, ro, na> latVertex,
                  acc<double, ro, na> lonVertex,
                  acc<vlreal, wo, na> u,
                  acc<vlreal, wo, na> h);

}}}
