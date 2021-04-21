/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <mpasoflecsi/specialization/mesh.hh>
#include "../state.hh"

namespace mpas { namespace ocean { namespace task {

void init_maxLevelEdgeTop(mesh::accessor<flecsi::ro, flecsi::ro> m,
                          flecsi::field<int>::accessor<flecsi::ro, flecsi::ro>    maxLevelCell,
                          flecsi::field<int>::accessor<flecsi::wo, flecsi::na>    maxLevelEdgeTop);

void init_layerThicknessEdge(mesh::accessor<flecsi::ro, flecsi::ro> m,
                             flecsi::field<vlreal>::accessor<flecsi::ro, flecsi::ro> layerThickness,
                             flecsi::field<int>::accessor<flecsi::ro, flecsi::ro>    maxLevelEdgeTop,
                             flecsi::field<vlreal>::accessor<flecsi::wo, flecsi::na> layerThicknessEdge); 

void init_meshScalingDel2(mesh::accessor<flecsi::ro, flecsi::ro> m,
                          flecsi::field<double>::accessor<flecsi::ro, flecsi::ro>    meshDensity,
                          flecsi::field<double>::accessor<flecsi::wo, flecsi::na>    meshScalingDel2);

void init_edgeSignOnCell(mesh::accessor<flecsi::ro, flecsi::ro> m,
                         flecsi::field<metensor<double>>::accessor<flecsi::wo, flecsi::na>    edgeSignOnCell);

void init_2nd_deriv(mesh::accessor<flecsi::ro, flecsi::ro> m,
                    flecsi::field<double>::accessor<flecsi::ro, flecsi::ro>    angleEdge,
                    flecsi::field<double>::accessor<flecsi::ro, flecsi::ro>    dcEdge,
                    flecsi::field<double>::accessor<flecsi::ro, flecsi::ro>    xCell,
                    flecsi::field<double>::accessor<flecsi::ro, flecsi::ro>    yCell,
                    flecsi::field<double>::accessor<flecsi::ro, flecsi::ro>    zCell,
                    flecsi::field<double>::accessor<flecsi::ro, flecsi::ro>    xVertex,
                    flecsi::field<double>::accessor<flecsi::ro, flecsi::ro>    yVertex,
                    flecsi::field<double>::accessor<flecsi::ro, flecsi::ro>    zVertex,
                    flecsi::field<derivtensor>::accessor<flecsi::wo, flecsi::na>    derivTwo);

void init_tracer_adv_coeff(mesh::accessor<flecsi::ro, flecsi::ro> m,
                           flecsi::field<int>::accessor<flecsi::ro, flecsi::ro>    maxLevelCell,
                           flecsi::field<double>::accessor<flecsi::ro, flecsi::ro>    dcEdge,
                           flecsi::field<double>::accessor<flecsi::ro, flecsi::ro>    dvEdge,
                           flecsi::field<vltensor<int>>::accessor<flecsi::ro, flecsi::ro>    boundaryCell,
                           flecsi::field<derivtensor>::accessor<flecsi::ro, flecsi::ro>    derivTwo,
                           flecsi::field<int>::accessor<flecsi::wo, flecsi::na>    nAdvCellsForEdge,
                           flecsi::field<vltensor<int>>::accessor<flecsi::wo, flecsi::na>    highOrderAdvMask,
                           flecsi::field<me2tensor<int>>::accessor<flecsi::wo, flecsi::na>    advCellsForEdge,
                           flecsi::field<me2tensor<double>>::accessor<flecsi::wo, flecsi::na>    advCoef,
                           flecsi::field<me2tensor<double>>::accessor<flecsi::wo, flecsi::na>    advCoef_3rd);


}}}
