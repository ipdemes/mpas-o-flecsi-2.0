#include "../state.hh"
#include "solve.hh"

#include <stdio.h>

namespace mpas { namespace ocean { namespace task {

using namespace flecsi;

// 3rd/4th order blending coefficient,  TODO: this should be an input, hard-coded for now.
const double coef3rdOrder = 0.25;

void copy_tracers(mesh::accessor<flecsi::ro, flecsi::ro> m,
                  acc<vltracer, ro, ro> tracers,
                  acc<vltracer, wo, na> tracersNew,
                  acc<vltracer, wo, na> tracersPrev)
{
  for(auto c: m.cells()) {
    for(std::size_t k = 0; k < nVertLevels; k++) {
      for(std::size_t t = 0; t < maxTracers; t++) {
        tracersNew(c)[k][t] = tracers(c)[k][t];
        tracersPrev(c)[k][t] = tracers(c)[k][t];
      }
    }
  }
}

void reinitialize_tracer_tendencies(mesh::accessor<flecsi::ro, flecsi::ro> m,
                                    acc<vltracer, wo, na> tend)
{
  for(auto c: m.cells()) {
    for(std::size_t k = 0; k < nVertLevels; k++) {
      for(std::size_t t = 0; t < maxTracers; t++) {
        tend(c)[k][t] = 0.0;
      }
    }
  }
}

void calc_normal_thickness_flux(mesh::accessor<flecsi::ro, flecsi::ro> m,
                                acc<vltensor<double>, ro, ro> layerThicknessEdge,
                                acc<vltensor<double>, ro, ro> normalVelocity,
                                acc<vltensor<double>, wo, na> normalThicknessFlux)
{

//FILE * outfile;
//outfile = fopen("debug.out","w");

  for(auto e: m.edges()) {
    for(std::size_t k = 0; k < nVertLevels; k++) {
      normalThicknessFlux(e)[k] = normalVelocity(e)[k] * layerThicknessEdge(e)[k];
//      fprintf(outfile,"   %8d   % 18.15e   % 18.15e   % 18.15e\n", e, normalThicknessFlux(e)[k], normalVelocity(e)[k], layerThicknessEdge(e)[k]);
    }
  }

//fclose(outfile);
}


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
                                        double dt)
{

  for(auto c: m.cells()) {
    double dtInvAreaCell = dt / areaCell(c);
    std::size_t kmax = maxLevelCell(c);
    for(std::size_t k = 0; k < kmax; k++) {
      hProv(c)[k] = layerThickness(c)[k];
    }

    std::size_t je = 0;
    for(auto e: m.edges(c)) {
      double signedFactor = dtInvAreaCell * dvEdge(e) * edgeSignOnCell(c)[je];
      // Provisional layer thickness is after horizontal thickness flux only
      for(std::size_t k = 0; k < kmax; k++) {
        hProv(c)[k] += signedFactor * normalThicknessFlux(e)[k];
      }
      je++;
    }

    for(std::size_t k = 0; k < kmax; k++) {
      hProvInv(c)[k] = 1.0 / hProv(c)[k];
      // New layer thickness is after horizontal and vertical thickness flux
      hNewInv(c)[k] = 1.0 / (hProv(c)[k]); //- dt * w(c)[k] + dt * w(c)[k+1]);
    }
  }
}

void calc_horiz_tracer_advection_flux_mono_low_order(mesh::accessor<flecsi::ro, flecsi::ro> m,
                                                     acc<double, ro, ro> dvEdge,
                                                     acc<vltensor<double>, ro, ro> normalThicknessFlux,
                                                     acc<vltracer, ro, ro> tracers,
                                                     acc<vltracer, wo, na> lowOrderFlux)
{
    // compute low order upwind fluxes
    for(auto e: m.edges()) {
    auto c1 = m.cells(e)[0];
    auto c2 = m.cells(e)[1];
    for(std::size_t k = 0; k < nVertLevels; k++) {
      double ntf1 = std::max(0.0,normalThicknessFlux(e)[k]);
      double ntf2 = std::min(0.0,normalThicknessFlux(e)[k]);
      for(std::size_t t = 0; t < maxTracers; t++) {
        lowOrderFlux(e)[k][t] = dvEdge(e) * (ntf1*tracers(c1)[k][t] + ntf2*tracers(c2)[k][t]);
      }
    }
  }

}

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
                                                      acc<vltracer, wo, na> highOrderFlux)
{
  std::vector<double> wgtTmp(nVertLevels);
  std::vector<double> sgnTmp(nVertLevels);

  std::vector<std::vector<double>> flxTmp(nVertLevels, std::vector<double>(maxTracers));

  for(auto e: m.edges()) {
    // compute some common intermediate factors
    for(std::size_t k = 0; k< nVertLevels; k++) {
      wgtTmp[k] = normalThicknessFlux(e)[k] * highOrderAdvMask(e)[k];
      sgnTmp[k] = (normalThicknessFlux(e)[k] < 0.0) ? -1.0 : (normalThicknessFlux(e)[k] > 0.);
      for(std::size_t t = 0; t < maxTracers; t++) {
        flxTmp[k][t] = 0.0;
      }
    }

    // compute 3rd/4th order fluxes where requested
    for(std::size_t i = 0; i < nAdvCellsForEdge(e); i++) {
      std::size_t iCell = advCellsForEdge(e)[i];
      double coef1 = advCoef(e)[i];
      double coef3 = advCoef_3rd(e)[i] * coef3rdOrder;
      for(std::size_t k = 0; k < maxLevelCell(iCell); k++) {
        double weightFactor = wgtTmp[k] * (coef1 + coef3*sgnTmp[k]);
        for(std::size_t t = 0; t < maxTracers; t++) {
          flxTmp[k][t] += tracers(iCell)[k][t] * weightFactor;
        }
      }
    }

    for(std::size_t k = 0; k < nVertLevels; k++) {
      for(std::size_t t = 0; t < maxTracers; t++) {
        highOrderFlux(e)[k][t] = flxTmp[k][t];
      }
    }
  }

}

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
                                                    double dt)
{
//FILE * outfile;
//outfile = fopen("debug1.out","w");

  const double eps = 1.e-10;

  // determine bounds on tracer (Min and Max) from surrounding cells
  for(auto c: m.cells()) {
    for(std::size_t k = 0; k < maxLevelCell(c); k++) {
      for(std::size_t t = 0; t < maxTracers; t++) {
        tracerMax(c)[k][t] = tracers(c)[k][t];
        tracerMin(c)[k][t] = tracers(c)[k][t];
        flxIn(c)[k][t] = 0.0;
        flxOut(c)[k][t] = 0.0;
      }
    }

    for(auto c2: m.cells(c)) {
      std::size_t kmax = std::min(maxLevelCell(c), maxLevelCell(c2));
      for(std::size_t k = 0; k < kmax; k++) {
        for(std::size_t t = 0; t < maxTracers; t++) {
          tracerMax(c)[k][t] = std::max(tracerMax(c)[k][t], tracers(c2)[k][t]);
          tracerMin(c)[k][t] = std::min(tracerMin(c)[k][t], tracers(c2)[k][t]);
        }
      }
    }
  }

  for(auto e: m.edges()) {
    auto c1 = m.cells(e)[0];
    auto c2 = m.cells(e)[1];
    // remove low order flux from high order flux for FCT
    // store left over high order flux in highOrderFlx array
    for(std::size_t k = 0; k < maxLevelEdgeTop(e); k++) {
      double tracerWeight = (highOrderAdvMask(e)[k]+1 & 1) * (dvEdge(e) * 0.5) *
                            normalThicknessFlux(e)[k];
      for(std::size_t t = 0; t < maxTracers; t++) {
        highOrderFlux(e)[k][t] += tracerWeight * (tracers(c1)[k][t] + tracers(c2)[k][t]);
        highOrderFlux(e)[k][t] -= lowOrderFlux(e)[k][t];
      }
    }
  }

  // accumulate remaining high order fluxes
  for(auto c: m.cells()) {
    double invAreaCell1 = 1.0 / areaCell(c);

    std::size_t je = 0;
    for(auto e: m.edges(c)) {
      double signedFactor = edgeSignOnCell(c)[je] * invAreaCell1;
      for(std::size_t k = 0; k < maxLevelEdgeTop(e); k++) {
        for(std::size_t t = 0; t < maxTracers; t++) {
          // accumulate upwind low order fluxes into tend array
          tend(c)[k][t] += signedFactor * lowOrderFlux(e)[k][t];

          // accumulate remaining high order fluxes
          flxIn(c)[k][t]  += std::max(0.0, signedFactor * highOrderFlux(e)[k][t]);
          flxOut(c)[k][t] += std::min(0.0, signedFactor * highOrderFlux(e)[k][t]);
        }
      }
      je++;
    }

    // build the factors for FCT computed using previously computed bounds
    // and the bound on the newly updated value factors are placed in flxIn and flxOut arrays
    for(std::size_t k = 0; k < maxLevelCell(c); k++) {
      for(std::size_t t = 0; t < maxTracers; t++) {
        double tracerUpwindNew = (tracers(c)[k][t] * layerThickness(c)[k] 
                               + dt * tend(c)[k][t]) * hProvInv(c)[k];
        double tracerMaxNew = tracerUpwindNew + dt * flxIn(c)[k][t] * hProvInv(c)[k];
        double tracerMinNew = tracerUpwindNew + dt * flxOut(c)[k][t] * hProvInv(c)[k];
        double scaleFactor = (tracerMax(c)[k][t] - tracerUpwindNew)
                           / (tracerMaxNew - tracerUpwindNew + eps);
        flxIn(c)[k][t] = std::min(1.0, std::max(0.0, scaleFactor));
        scaleFactor = (tracerUpwindNew - tracerMin(c)[k][t])
                    / (tracerUpwindNew - tracerMinNew + eps);
        flxOut(c)[k][t] = std::min(1.0, std::max(0.0, scaleFactor));
      }
    }
  }

  // rescale horizontal fluxes
  for(auto e: m.edges()) {
    auto c1 = m.cells(e)[0];
    auto c2 = m.cells(e)[1];
//    fprintf(outfile, "   %8d   %8d   %8d\n", e, c1, c2);
    for(std::size_t k = 0; k < maxLevelEdgeTop(e); k++) {
      for(std::size_t t = 0; t < maxTracers; t++) {
        highOrderFlux(e)[k][t] = std::max(0.0, highOrderFlux(e)[k][t])
                               * std::min(flxOut(c1)[k][t], flxIn(c2)[k][t])
                               + std::min(0.0, highOrderFlux(e)[k][t])
                               * std::min(flxIn(c1)[k][t], flxOut(c2)[k][t]);
      }
    }
  }
//fclose(outfile);
//exit(0);
} // end correct_and_rescale_high_order_mono_horiz_flux

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
                                                     double dt)
{
  // accumulate the scaled high order horizontal tendencies, and the upwind tendencies
  for(auto c: m.cells()) {
    double invAreaCell1 = 1.0 / areaCell(c);

    std::size_t je = 0;
    for(auto e: m.edges(c)) {
      double signedFactor = edgeSignOnCell(c)[je] * invAreaCell1;

      for(std::size_t k = 0; k < maxLevelEdgeTop(e); k++) {
        for(std::size_t t = 0; t < maxTracers; t++) {
          // tend before is the upwind tendency
          // tend after is the total horizontal advection tendency
          tend(c)[k][t] += signedFactor * highOrderFlux(e)[k][t];
        }
      }
      je++;
    }
    // save provisional tracer after horizontal fluxes only
    for(std::size_t k = 0; k < maxLevelCell(c); k++) {
      for(std::size_t t = 0; t < maxTracers; t++) {
        tracersProv(c)[k][t] = (tracers(c)[k][t] * layerThickness(c)[k]
          + dt * tend(c)[k][t]) * hProvInv(c)[k];
      }
    }
  }

}

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
                            double coef)
{
  for(auto c: m.cells()) {
    double invArea = 1.0 / areaCell(c);

    std::size_t je = 0;
    for(auto e: m.edges(c)) {
      auto c1 = m.cells(e)[0];
      auto c2 = m.cells(e)[1];

      double rho_scaling = meshScalingDel2(e) * coef * dvEdge(e) / dcEdge(e);

      for(std::size_t k = 0; k < maxLevelEdgeTop(e); k++) {
        for(std::size_t t = 0; t < maxTracers; t++) {

          // \kappa_2 \nabla \phi on edge
          double tracer_turb_flux = tracers(c2)[k][t] - tracers(c1)[k][t];

          // div(h \kappa_2 \nabla \phi) at cell center
          double flux = layerThicknessEdge(e)[k] * tracer_turb_flux * rho_scaling;

          tend(c)[k][t] -= edgeSignOnCell(c)[je] * flux * invArea;
        }
      }
      je++;
    }
  }

}

void step_tracers_forward(mesh::accessor<flecsi::ro, flecsi::ro> m,
                          acc<vltracer, ro, ro> tend,
                          acc<vltracer, ro, ro> tracers,
                          acc<vltracer, rw, ro> tracersNew,
                          acc<vltracer, wo, na> tracersPrev,
                          double rk_weight,
                          double rk_subweight)
{

  for(auto c: m.cells()) {
    for(std::size_t j = 0; j < nVertLevels; j++) {
      for(std::size_t k = 0; k < maxTracers; k++) {
        //This calculates the new tracer value that will accumulate.
        tracersNew(c)[j][k] += rk_weight * tend(c)[j][k];
        //This calculates the tracer value to be used for calculating next tend.
        tracersPrev(c)[j][k] = tracers(c)[j][k] + rk_subweight * tend(c)[j][k];
      }
    }
  }

//if(rk_subweight == 0.0)
//{                             
//FILE * outfile;
//outfile = fopen("tracers.out","w");
//
//  for(auto c: m.cells()) {
//    fprintf(outfile, "   %8d   % 18.15e   % 18.15e   % 18.15e\n", c, tracers(c)[0][0], tracersNew(c)[0][0], tracersPrev(c)[0][0]);
//  }
//fclose(outfile);
//}

}

void update_tracers(mesh::accessor<flecsi::ro, flecsi::ro> m,
                          acc<vltracer, ro, ro> tracersNew,
                          acc<vltracer, wo, na> tracers)
{

  for(auto c: m.cells()) {
    for(std::size_t k = 0; k < nVertLevels; k++) {
      for(std::size_t t = 0; t < maxTracers; t++) {
        tracers(c)[k][t] = tracersNew(c)[k][t];
      }
    }
  }

}

}}}
