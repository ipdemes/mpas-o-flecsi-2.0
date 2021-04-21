/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "mpasoflecsi/specialization/mesh.hh"

namespace mpas { namespace ocean { namespace task {

void copy_tracers(mesh::accessor<flecsi::ro, flecsi::ro> m,
                  acc<vltracer, ro, ro> tracers,
                  acc<vltracer, wo, na> tracersNew,
                  acc<vltracer, wo, na> tracersPrev);

void reinitialize_tracer_tendencies(mesh::accessor<flecsi::ro, flecsi::ro> m,
                                    acc<vltracer, wo, na> tend);

void calc_normal_thickness_flux(mesh::accessor<flecsi::ro, flecsi::ro> m,
                                acc<vltensor<double>, ro, ro> layerThicknessEdge,
                                acc<vltensor<double>, ro, ro> normalVelocity,
                                acc<vltensor<double>, wo, na> normalThicknessFlux); 

void calc_provisional_layer_thicknesses(mesh::accessor<flecsi::ro, flecsi::ro> m,
                                        acc<int, ro, ro> maxLevelCell,
                                        acc<double, ro, ro> areaCell,
                                        acc<double, ro, ro> dvEdge,
                                        acc<metensor<double>, ro, ro> edgeSignOnCell,
                                        acc<vltensor<double>, ro, ro> layerThickness,
                                        acc<vltensor<double>, ro, ro> normalThicknessFlux,
                                        acc<vltensor<double>, ro, ro> w,
                                        acc<vltensor<double>, wo, na> hProv,
                                        acc<vltensor<double>, wo, na> hProvInv,
                                        acc<vltensor<double>, wo, na> hNewInv,
                                        double dt);

void calc_horiz_tracer_advection_flux_mono_low_order(mesh::accessor<flecsi::ro, flecsi::ro> m,
                                                     acc<double, ro, ro> dvEdge,
                                                     acc<vltensor<double>, ro, ro> normalThicknessFlux,
                                                     acc<vltracer, ro, ro> tracers,
                                                     acc<vltracer, wo, na> lowOrderFlux);

void calc_horiz_tracer_advection_flux_mono_high_order(mesh::accessor<flecsi::ro, flecsi::ro> m,
                                                      acc<vltensor<int>, ro, ro> highOrderAdvMask,
                                                      acc<double, ro, ro> dvEdge,
                                                      acc<vltensor<double>, ro, ro> normalThicknessFlux,
                                                      acc<me2tensor<double>, ro, ro> advCoef,
                                                      acc<me2tensor<double>, ro, ro> advCoef_3rd,
                                                      acc<int, ro, ro> maxLevelCell,
                                                      acc<int, ro, ro> nAdvCellsForEdge,
                                                      acc<me2tensor<int>, ro, ro> advCellsForEdge,
                                                      acc<vltracer, ro, ro> tracers,
                                                      acc<vltracer, wo, na> highOrderFlux);

void correct_and_rescale_high_order_mono_horiz_flux(mesh::accessor<flecsi::ro, flecsi::ro> m,
                                                    acc<metensor<double>, ro, ro> edgeSignOnCell,
                                                    acc<vltensor<double>, ro, ro> hProvInv,
                                                    acc<vltensor<double>, ro, ro> layerThickness,
                                                    acc<vltensor<double>, ro, ro> normalThicknessFlux,
                                                    acc<double, ro, ro> areaCell,
                                                    acc<double, ro, ro> dvEdge,
                                                    acc<vltensor<int>, ro, ro> highOrderAdvMask,
                                                    acc<int, ro, ro> maxLevelCell,
                                                    acc<int, ro, ro> maxLevelEdgeTop,
                                                    acc<vltracer, ro, ro> lowOrderFlux,
                                                    acc<vltracer, ro, ro> tracers,
                                                    acc<vltracer, rw, ro> highOrderFlux,
                                                    acc<vltracer, wo, na> flxIn,
                                                    acc<vltracer, wo, na> flxOut,
                                                    acc<vltracer, wo, na> tend,
                                                    acc<vltracer, wo, na> tracerMax,
                                                    acc<vltracer, wo, na> tracerMin,
                                                    double dt);

void accumulate_horiz_tracer_advection_tend_mono_FCT(mesh::accessor<flecsi::ro, flecsi::ro> m,
                                                     acc<metensor<double>, ro, ro> edgeSignOnCell,
                                                     acc<vltensor<double>, ro, ro> hProvInv,
                                                     acc<vltensor<double>, ro, ro> layerThickness,
                                                     acc<double, ro, ro> areaCell,
                                                     acc<int, ro, ro> maxLevelCell,
                                                     acc<int, ro, ro> maxLevelEdgeTop,
                                                     acc<vltracer, ro, ro> highOrderFlux,
                                                     acc<vltracer, ro, ro> tracers,
                                                     acc<vltracer, wo, na> tracersProv,
                                                     acc<vltracer, rw, ro> tend,
                                                     double dt);

void tracer_hmix_del2_tend(mesh::accessor<flecsi::ro, flecsi::ro> m,
                            acc<double, ro, ro> areaCell,
                            acc<double, ro, ro> dcEdge,
                            acc<double, ro, ro> dvEdge,
                            acc<double, ro, ro> meshScalingDel2,
                            acc<int, ro, ro> maxLevelEdgeTop,
                            acc<metensor<double>, ro, ro> edgeSignOnCell,
                            acc<vltensor<double>, ro, ro> layerThicknessEdge,
                            acc<vltracer, ro, ro> tracers,
                            acc<vltracer, rw, ro> tend,
                            double coef);

void step_tracers_forward(mesh::accessor<flecsi::ro, flecsi::ro> m,
                          acc<vltracer, ro, ro> tend,
                          acc<vltracer, ro, ro> tracers,
                          acc<vltracer, rw, ro> tracersNew,
                          acc<vltracer, wo, na> tracersPrev,
                          double rk_weight,
                          double rk_subweight);

void update_tracers(mesh::accessor<flecsi::ro, flecsi::ro> m,
                          acc<vltracer, ro, ro> tracersNew,
                          acc<vltracer, wo, na> tracers);


}}}
