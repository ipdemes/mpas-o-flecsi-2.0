/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "mpasoflecsi/common/constants.hh"

#include "../state.hh"
#include "testcase.hh"
#include "mpasoflecsi/eqns/metrics.hh"

namespace mpas { namespace ocean { namespace task {

using namespace flecsi;

double rotation_velocity(double latE, double angleE)
{
  using namespace mpas_constants;
  constexpr double u0 = 2.0 * pi * a / 86400.0;  // velocity parameter

  return u0 * std::cos(latE) * std::cos(angleE);
}

double initial_gaussian_dist(double xC)
{
  using namespace mpas_constants;

  constexpr double maxG = 21983.7698959881;      // max value of unnormalized Gaussian
  constexpr double Knorm = 20.0/maxG;            // normalization parameter
  constexpr double kappa = 10.0;                 // concentration parameter
  constexpr double T0 = 50.0;                    // base tracer value

  return T0 + Knorm * std::exp(kappa * xC / a);
}

double initial_slotCyl_dist(double xC, double yC, double zC, double latC, double lonC)
{
  using namespace mpas_constants;

  double aa = a/2.;
  double hw = a/12.;
  double bot = -a/3.;

  double dist = eqns::sphere_arc_length(xC, yC, zC, a, 0., 0.);

  double temp;

  double longC = lonC;
  if(longC > pi) longC -= 2.*pi;

  if(dist > aa){
    temp = 50.;
  } else if( (a*latC > bot) && (a*latC < aa) && (a*std::abs(longC) < hw) ){
    temp = 50.;
  } else {
    temp = 100.;
  }

  return temp;
}

double spherical_harmonic_dist(double latC, double lonC, double coef, double elapsed)
{
  using namespace mpas_constants;

  double cosLat = std::cos(latC);
  double ampl = 20.;
  double T0 = 50.;

  return T0 + ampl * cosLat * cosLat * std::cos(2. * lonC)
         * std::exp(-6.0 * coef * elapsed / (a*a));

}

// init w, boundaryCell, boundaryLevelDepth, surfaceStress, maxLevelCell
void init_extra_fields(mesh::accessor<ro, ro> m,
                              field<vltensor<double>>::accessor<wo, na> w,
                              field<vltensor<int>>::accessor<wo, na> bc,
                              field<double>::accessor<wo, na> bld,
                              field<double>::accessor<wo, na> ss,
                              field<int>::accessor<wo, na> mlc) {
  for(auto c: m.cells()) {
    bld(c) = 0.;
    mlc(c) = 1;
    for(std::size_t k = 0; k < nVertLevels; k++) {
      w(c)[k] = 0.;
      bc(c)[k] = 0;
    }
  }
  for(auto e: m.edges()) {
    ss(e) = 0.;
  }

}

void setup_case_1(mesh::accessor<ro, ro> m,
                  field<double>::accessor<ro, ro> angleEdge,
                  field<double>::accessor<ro, ro> latEdge,
                  field<double>::accessor<ro, ro> xCell,  
                  field<vltensor<double>>::accessor<wo, na> u,
                  field<vltensor<double>>::accessor<wo, na> h,
                  field<vltracer>::accessor<wo, na> tracers)
{
  // initialize normal velocities, all latitudes have rotation period of 1 day
  for(auto e: m.edges()){
    u(e)[0] = rotation_velocity(latEdge(e), angleEdge(e));
  }

  // set layer thickness to unity
  for(auto c: m.cells()){
    h(c)[0] = 1.;
  }

  // initialize tracer as Gaussian hill
  for(auto c: m.cells()){
    tracers(c)[0][0] = initial_gaussian_dist(xCell(c));
  }
}

void setup_case_2(mesh::accessor<ro, ro> m,
                  field<double>::accessor<ro, ro> angleEdge,
                  field<double>::accessor<ro, ro> latEdge,
                  field<double>::accessor<ro, ro> latCell,
                  field<double>::accessor<ro, ro> lonCell,
                  field<double>::accessor<ro, ro> xCell,  
                  field<double>::accessor<ro, ro> yCell,
                  field<double>::accessor<ro, ro> zCell,
                  field<vltensor<double>>::accessor<wo, na> u,
                  field<vltensor<double>>::accessor<wo, na> h,
                  field<vltracer>::accessor<wo, na> tracers)
{
  // initialize normal velocities, all latitudes have rotation period of 1 day
  for(auto e: m.edges()){
    u(e)[0] = rotation_velocity(latEdge(e), angleEdge(e));
  }

  // set layer thickness to unity
  for(auto c: m.cells()){
    h(c)[0] = 1.;
  }

  // initialize tracer as slotted cylinder
  for(auto c: m.cells()){
    tracers(c)[0][0] = initial_slotCyl_dist(xCell(c), yCell(c), zCell(c), latCell(c), lonCell(c));
  }

}

void setup_case_3(mesh::accessor<ro, ro> m,
                  field<double>::accessor<ro, ro> angleEdge,
                  field<double>::accessor<ro, ro> latEdge,
                  field<double>::accessor<ro, ro> latCell,
                  field<double>::accessor<ro, ro> lonCell,
                  field<vltensor<double>>::accessor<wo, na> u,
                  field<vltensor<double>>::accessor<wo, na> h,
                  field<vltracer>::accessor<wo, na> tracers,
                  double coef)
{
  // initialize normal velocities at rest
  for(auto e: m.edges()){
    u(e)[0] = 0.0; 
  }

  // set layer thickness to unity
  for(auto c: m.cells()){
    h(c)[0] = 1.;
  }

  // initialize tracer as spherical harmonic
  for(auto c: m.cells()){
    tracers(c)[0][0] = spherical_harmonic_dist(latCell(c), lonCell(c), coef, 0.0);
  }
}

void calc_error_case_1(mesh::accessor<ro, ro> m,
                      field<double>::accessor<ro, ro> latCell,
                      field<double>::accessor<ro, ro> lonCell,
                      field<double>::accessor<ro, ro> xCell,
                      field<double>::accessor<ro, ro> yCell,
                      field<vltracer>::accessor<ro, ro> tracers,
                      double timeElapsed)
{
  using namespace mpas_constants;
  double dayFrac = timeElapsed / 86400.;
  double rotAng = 2.0 * pi * dayFrac;

  double errSqd = 0.0;
  double errMax = 0.0;

  auto nCells = m.cells().size();

  FILE *fp;
  fp = fopen("calcEnd_gauss.dat", "w");

  for(auto c: m.cells()) {
    double x_prime = xCell(c)*std::cos(-rotAng) - yCell(c)*std::sin(-rotAng);
    double trueAnswer = initial_gaussian_dist(x_prime);

    double error = (tracers(c)[0][0] - trueAnswer) / trueAnswer;
    fprintf(fp, "% 18.15e    % 18.15e    % 18.15e    % 18.15e    % 18.15e\n", latCell(c), lonCell(c), tracers(c)[0][0], trueAnswer, error);

    if(std::abs(error) > errMax) errMax = error;
    errSqd += error * error;
  }

  double totalError = std::sqrt(errSqd) / nCells;
  flog(info) << "Total L2 RMS Norm of the error was: "<< totalError << std::endl;
  flog(info) << " Max error:  "<< std::abs(errMax) << std::endl;

  fclose(fp);

}

void calc_error_case_2(mesh::accessor<ro, ro> m,
                      field<double>::accessor<ro, ro> latCell,
                      field<double>::accessor<ro, ro> lonCell,
                      field<double>::accessor<ro, ro> xCell,
                      field<double>::accessor<ro, ro> yCell,
                      field<double>::accessor<ro, ro> zCell,
                      field<vltracer>::accessor<ro, ro> tracers,
                      double timeElapsed)
{
  using namespace mpas_constants;
  double dayFrac = timeElapsed / 86400.;
  double rotAng = 2.0 * pi * dayFrac;

  double errSqd = 0.0;
  double errMax = 0.0;

  auto nCells = m.cells().size();

  FILE *fp;
  fp = fopen("calcEnd_slotCyl.dat", "w");

  for(auto c: m.cells()) {
    double x_prime = xCell(c)*std::cos(rotAng) + yCell(c)*std::sin(rotAng);
    double y_prime =-xCell(c)*std::sin(rotAng) + yCell(c)*std::cos(rotAng);

    double lon_prime = lonCell(c) - rotAng;

    double trueAnswer = initial_slotCyl_dist(x_prime, y_prime, zCell(c), latCell(c), lon_prime);

    double error = (tracers(c)[0][0] - trueAnswer) / trueAnswer;

    fprintf(fp, "% 18.15e    % 18.15e    % 18.15e    % 18.15e    % 18.15e\n", latCell(c), lonCell(c), tracers(c)[0][0], trueAnswer, error);

    if(std::abs(error) > errMax) errMax = error;
    errSqd += error * error;
  }

  double totalError = std::sqrt(errSqd) / nCells;
  flog(info) << "Total L2 RMS Norm of the error was: "<< totalError << std::endl;
  flog(info) << " Max error:  "<< std::abs(errMax) << std::endl;

  fclose(fp);

}

void calc_error_case_3(mesh::accessor<ro, ro> m,
                      field<double>::accessor<ro, ro> latCell,
                      field<double>::accessor<ro, ro> lonCell,
                      field<vltracer>::accessor<ro, ro> tracers,
                      double hmix_coef,
                      double timeElapsed)
{
  using namespace mpas_constants;

  double errSqd = 0.0;
  double errMax = 0.0;

  auto nCells = m.cells().size();

  FILE *fp;
  fp = fopen("calcEnd_diffus.dat", "w");

  for(auto c: m.cells()) {

    double trueAnswer = spherical_harmonic_dist(latCell(c), lonCell(c), hmix_coef, timeElapsed);

    double error = (tracers(c)[0][0] - trueAnswer) / trueAnswer;

    fprintf(fp, "% 18.15e    % 18.15e    % 18.15e    % 18.15e    % 18.15e\n", latCell(c), lonCell(c), tracers(c)[0][0], trueAnswer, error);

    if(std::abs(error) > errMax) errMax = error;
    errSqd += error * error;
  }

  double totalError = std::sqrt(errSqd) / nCells;
  flog(info) << "Total L2 RMS Norm of the error was: "<< totalError << std::endl;
  flog(info) << " Max error:  "<< std::abs(errMax) << std::endl;

  fclose(fp);

}

}}}
