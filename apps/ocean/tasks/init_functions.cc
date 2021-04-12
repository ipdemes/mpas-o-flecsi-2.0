/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include <stdlib.h>

#include "init_functions.hh"
#include "mpasoflecsi/eqns/matrix.hh"
#include "mpasoflecsi/eqns/metrics.hh"

#include <algorithm>

namespace mpas { namespace ocean { namespace task {

using namespace flecsi;

void init_maxLevelEdgeTop(mesh::accessor<ro, ro> m,
                          field<int>::accessor<ro, ro>    maxLevelCell,
                          field<int>::accessor<wo, na>    maxLevelEdgeTop)
{
  for(auto e: m.edges()) {
    auto cell1 = m.cells(e)[0];
    auto cell2 = m.cells(e)[1];
    maxLevelEdgeTop(e) = std::max(maxLevelCell(cell1), maxLevelCell(cell2));
  }
}

// TODO: doublecheck boundary cells
void init_layerThicknessEdge(mesh::accessor<ro, ro> m,
                             field<vlreal>::accessor<ro, ro>    layerThickness,  
                             field<int>::accessor<ro, ro>    maxLevelEdgeTop,
                             field<vlreal>::accessor<wo, na>    layerThicknessEdge)
{
  for(auto e: m.edges()) {
    for(std::size_t k = 0; k < nVertLevels; k++) {
      layerThicknessEdge(e)[k] = 0.0;
    }
    auto cell1 = m.cells(e)[0];
    auto cell2 = m.cells(e)[1];
    for(std::size_t k = 0; k < nVertLevels; k++) {
      layerThicknessEdge(e)[k] = 0.5 * (layerThickness(cell1)[k] + layerThickness(cell2)[k]);
    }
  }
}

void init_meshScalingDel2(mesh::accessor<ro, ro> m,
                          field<double>::accessor<ro, ro>    meshDensity,
                          field<double>::accessor<wo, na>    meshScalingDel2)
{
  double maxDensity = 0.0;
  for(auto c: m.cells()) {
    if(meshDensity(c) > maxDensity) {
      maxDensity = meshDensity(c);
    }
  }

  double relativeDensity;
  for(auto e: m.edges()) {
    auto cell1 = m.cells()[0];
    auto cell2 = m.cells()[1];
    relativeDensity = (meshDensity(cell1) + meshDensity(cell2))
    / (2.0*maxDensity);
    meshScalingDel2(e) = 1.0 / std::pow(relativeDensity, 0.25);
  }
}

void init_edgeSignOnCell(mesh::accessor<ro, ro> m,
                         field<metensor<double>>::accessor<wo, na>    edgeSignOnCell)
{
  for(auto c: m.cells()) {
    auto nEonC = m.edges(c).size();
    std::size_t j = 0;
    for(auto e: m.edges(c)) {
      if(c == m.cells(e)[0])
        edgeSignOnCell(c)[j] = -1.0;
      else
        edgeSignOnCell(c)[j] = 1.0;
      j++;
    }
  }
}
void init_2nd_deriv(mesh::accessor<ro, ro> m,
                    field<double>::accessor<ro, ro>    angleEdge,
                    field<double>::accessor<ro, ro>    dcEdge,
                    field<double>::accessor<ro, ro>    xCell,
                    field<double>::accessor<ro, ro>    yCell,
                    field<double>::accessor<ro, ro>    zCell,
                    field<double>::accessor<ro, ro>    xVertex,
                    field<double>::accessor<ro, ro>    yVertex,
                    field<double>::accessor<ro, ro>    zVertex,
                    field<derivtensor>::accessor<wo, na>    derivTwo)
{
  using namespace mpas_constants;

  const std::size_t polynomial_order = 2;
  const std::size_t ll = 25;

  for(auto e: m.edges()) {
    for(std::size_t j = 0; j < maxEdges2; j++) {
      for(std::size_t i = 0; i < 2; i++) {
        derivTwo(e)[j][i] = 0.0;
      }
    }
  }

  const int nEdges = m.edges().size();
  const int nCells = m.cells().size();

  std::vector<double> thetae1(nEdges);
  std::vector<double> thetae2(nEdges);

  // loop over every cell in mesh
  for(auto c: m.cells()) {

    // create list of this cell and all neighboring cells needed for constructing polynomial
    std::vector<int> cell_list{};

    const int nEonC = m.edges(c).size();

    cell_list.push_back(c);
//    fprintf(outfile, "   %8d", c);

    for(auto c2: m.cells(c)) {
      cell_list.push_back(c2);
//    fprintf(outfile, "   %8d", c2);
    }
//    fprintf(outfile, "\n");

    int n = nEonC + 1; // TODO: doublecheck - nEonC not necessarily = ncells on cell near boundaries
    bool do_the_cell = true;

    // check if all needed cell are available, otherwise skip the cell
    for(std::size_t i = 1; i < n; i++) {
      if(cell_list[i] < 0 or cell_list[i] >= nCells) {
        do_the_cell = false;
        break;
      }
    }

    if(not do_the_cell) continue;

    // TODO: this should be input, hard-coded for now.
    bool on_a_sphere = true;

    std::vector<double> xc{};
    std::vector<double> yc{};
    std::vector<double> zc{};

    double theta_abs;

    std::vector<double> theta_v{};
    std::vector<double> theta_t{};
    std::vector<double> dl_sphere{};

    std::vector<double> angle_2d{};

    std::vector<double> xp{};
    std::vector<double> yp{};

    if (on_a_sphere) {
      for(std::size_t i = 0; i < n; i++) {
        xc.push_back(xCell(cell_list[i]) / a);
        yc.push_back(yCell(cell_list[i]) / a);
        zc.push_back(zCell(cell_list[i]) / a);
      }

      if (zc[0] == 1.0)
        theta_abs = pi / 2.0;
      else
        theta_abs = pi / 2.0 - eqns::sphere_angle(
                    xc[0], yc[0], zc[0], xc[1], yc[1], zc[1], 0.0, 0.0, 1.0 );

      double length_scale = 1.0;

      // angles from cell center to neighbor centers (theta_v) & arc lengths from cell centers to neighbor centers (dl_sphere)
      for(std::size_t i = 0; i < n-1; i++) {
        std::size_t ip2 = i + 2;
        if (ip2 > n - 1) ip2 = 1;
        theta_v.push_back(eqns::sphere_angle(xc[0], yc[0], zc[0],
                          xc[i+1], yc[i+1], zc[i+1], xc[ip2], yc[ip2], zc[ip2]));
        dl_sphere.push_back(a * eqns::sphere_arc_length(xc[0], yc[0], zc[0],
                            xc[i+1], yc[i+1], zc[i+1]) / length_scale);
      }

      theta_t.push_back(theta_abs);

      for(std::size_t i = 1; i < n-1; i++) {
        theta_t.push_back(theta_t[i-1] + theta_v[i-1]);
      }

      // points projected onto tangent plane
      for(std::size_t i = 0; i < n-1; i++) {
        xp.push_back(std::cos(theta_t[i]) * dl_sphere[i]);
        yp.push_back(std::sin(theta_t[i]) * dl_sphere[i]);
      }

    } else { // on a plane

    }
////////////////////////////////////////////////////////////////////////////////
////////! translated from mpas_poly_fit_2() in MPAS-Model/src/operators/mpas_geometry_utils.F
////////!   solve for bmatrx (polynomial coefficients)
//////////////////////////////////////////////////////////////////////////////////

    std::size_t na;
    std::size_t ma = n;
    std::size_t mw = nEonC;

    double amatrix[ll][ll] = {0}, bmatrix[ll][ll] = {0}, wmatrix[ll][ll] = {0};

    // higher order polynomials than 2nd are not fully implemented in Fortran version, but
    // the guts are there, leaving open the option for future development here
    if (polynomial_order == 2) {
      na = 6;

      amatrix[0][0] = 1.0;
      wmatrix[0][0] = 1.0;

      for(std::size_t i = 1; i < ma; i++) {
        amatrix[i][0] = 1.0;
        amatrix[i][1] = xp[i-1];
        amatrix[i][2] = yp[i-1];
        amatrix[i][3] = xp[i-1]*xp[i-1];
        amatrix[i][4] = xp[i-1]*yp[i-1];
        amatrix[i][5] = yp[i-1]*yp[i-1];

        wmatrix[i][i] = 1.0;
      }
    }

    std::vector<std::vector<double>> aa, ww;

    for(std::size_t j = 0; j < ma; j++){
      std::vector<double> atmp;
      for(std::size_t i = 0; i < na; i++){
        atmp.push_back(amatrix[j][i]);
      }
      aa.push_back(atmp);
    }

    for(std::size_t j = 0; j < ma; j++){
      std::vector<double> wtmp;
      for(std::size_t i = 0; i < ma; i++){
        wtmp.push_back(wmatrix[j][i]);
      }
      ww.push_back(wtmp);
    }

    auto at = eqns::transpose(aa);
    auto wt = eqns::transpose(ww);

    auto h = eqns::matMul(wt,ww);

    auto ath = eqns::matMul(at,h);

    auto atha = eqns::matMul(ath,aa);

    auto atha_inv = eqns::migs(atha);

    auto bb = eqns::matMul(atha_inv,ath);

    for(std::size_t j = 0; j < na; j++){
      for(std::size_t i = 0; i < ma; i++){
        bmatrix[j][i] = bb[j][i];
      }
    }
////////////////////////////////////////////////////////////////////////////////
////////!  end mpas_poly_fit_2
////////////////////////////////////////////////////////////////////////////////

    for(std::size_t i=0; i < nEonC; i++) {
      double xv1, yv1, zv1, xv2, yv2, zv2;
      auto e1 = m.edges(c)[i];
      if( on_a_sphere ) {
        xv1 = xVertex(m.vertices(e1)[0]) / a;
        yv1 = yVertex(m.vertices(e1)[0]) / a;
        zv1 = zVertex(m.vertices(e1)[0]) / a;
        xv2 = xVertex(m.vertices(e1)[1]) / a;
        yv2 = yVertex(m.vertices(e1)[1]) / a;
        zv2 = zVertex(m.vertices(e1)[1]) / a;
      } else { // on a plane
      }

      double xec,yec,zec;
      double cos2t, sin2t, costsint;
      if( on_a_sphere ) {
        eqns::arc_bisect(xv1, yv1, zv1,
                               xv2, yv2, zv2,
                               xec, yec, zec);


        double thetae_tmp = eqns::sphere_angle( 
                            xc[0], yc[0], zc[0], xc[i+1], yc[i+1], zc[i+1],
                            xec, yec, zec);

        thetae_tmp += theta_t[i];

        cos2t = cos(thetae_tmp);
        sin2t = sin(thetae_tmp);
        costsint = cos2t*sin2t;
        cos2t = cos2t*cos2t;
        sin2t = sin2t*sin2t;

        if( c == m.cells(e1)[0] ){
          for(std::size_t j=0; j < n; j++){
            thetae1[e1] = thetae_tmp;
            derivTwo(e1)[j][0] = 2.*cos2t*bmatrix[3][j]
                               + 2.*costsint*bmatrix[4][j]
                               + 2.*sin2t*bmatrix[5][j];
          }
        } else {
          for(std::size_t j=0; j < n; j++){
            thetae2[e1] = thetae_tmp;
            derivTwo(e1)[j][1] = 2.*cos2t*bmatrix[3][j]
                               + 2.*costsint*bmatrix[4][j]
                               + 2.*sin2t*bmatrix[5][j];
          }
        }

      } else { // on a plane
      } // end if

    } // end of loop over edges on cells

  } // end of loop over cells

} // end init_2nd_deriv

void init_tracer_adv_coeff(mesh::accessor<ro, ro> m,
                           field<int>::accessor<ro, ro>    maxLevelCell,
                           field<double>::accessor<ro, ro>    dcEdge,
                           field<double>::accessor<ro, ro>    dvEdge,
                           field<vltensor<int>>::accessor<ro, ro>    boundaryCell,
                           field<derivtensor>::accessor<ro, ro>    derivTwo,
                           field<int>::accessor<wo, na>    nAdvCellsForEdge,
                           field<vltensor<int>>::accessor<wo, na>    highOrderAdvMask,
                           field<me2tensor<int>>::accessor<wo, na>    advCellsForEdge,
                           field<me2tensor<double>>::accessor<wo, na>    advCoef,
                           field<me2tensor<double>>::accessor<wo, na>    advCoef_3rd)
{
  for(auto e: m.edges()) {
    for(std::size_t i = 0; i < maxEdges2; i++) {
      advCoef(e)[i] = 0.;
      advCoef_3rd(e)[i] = 0.;
      advCellsForEdge(e)[i] = 0;
    }

    for(std::size_t k = 0; k < nVertLevels; k++) {
      highOrderAdvMask(e)[k] = 1;
      for(auto c: m.cells(e)) {
        if(boundaryCell(c)[k] == 1)
          highOrderAdvMask(e)[k] = 0;
      }
    }
  }

  for(auto e: m.edges()) {
    nAdvCellsForEdge(e) = 0;

    auto cell1 = m.cells(e)[0];
    auto cell2 = m.cells(e)[1];

    // Build unique list of cells used for advection on this edge, starting with cells touching edge
    std::vector<std::size_t> cell_list{cell1, cell2};

    // add new cells touching cell1
    std::size_t nEonC1 = m.edges(cell1).size();
    for(auto c: m.cells(cell1)) {
      if(std::find(cell_list.begin(), cell_list.end(), c) == cell_list.end())
        cell_list.push_back(c);
    }

    // add new cells touching cell2
    std::size_t nEonC2 = m.edges(cell2).size();
    for(auto c: m.cells(cell2)) {
      if(std::find(cell_list.begin(), cell_list.end(), c) == cell_list.end())
        cell_list.push_back(c);
    }

    nAdvCellsForEdge(e) = cell_list.size();
    for(std::size_t iCell = 0; iCell < nAdvCellsForEdge(e) ; iCell++) {
      advCellsForEdge(e)[iCell] = cell_list[iCell];
    }

    std::vector<std::size_t>::iterator it;
    std::size_t idx;

    // compute advection coefficients
    it = std::find(cell_list.begin(), cell_list.end(), cell1);
    idx = std::distance(cell_list.begin(), it);
    if(idx < nAdvCellsForEdge(e)) {
      advCoef(e)[idx] = advCoef(e)[idx] + derivTwo(e)[0][0];
      advCoef_3rd(e)[idx] = advCoef_3rd(e)[idx] + derivTwo(e)[0][0];
    }

    for(std::size_t iCell = 0; iCell < nEonC1; iCell++) {
      auto c = m.cells(cell1)[iCell];
      it = std::find(cell_list.begin(), cell_list.end(), c);
      idx = std::distance(cell_list.begin(), it);
      if(idx < nAdvCellsForEdge(e)) {
        advCoef(e)[idx] = advCoef(e)[idx] + derivTwo(e)[iCell+1][0];
        advCoef_3rd(e)[idx] = advCoef_3rd(e)[idx] + derivTwo(e)[iCell+1][0];
      }
    }

    it = std::find(cell_list.begin(), cell_list.end(), cell2);
    idx = std::distance(cell_list.begin(), it);
    if(idx < nAdvCellsForEdge(e)) {
      advCoef(e)[idx] = advCoef(e)[idx] + derivTwo(e)[0][1];
      advCoef_3rd(e)[idx] = advCoef_3rd(e)[idx] - derivTwo(e)[0][1];
    }

    for(std::size_t iCell = 0; iCell < nEonC2; iCell++) {
      auto c = m.cells(cell2)[iCell];
      it = std::find(cell_list.begin(), cell_list.end(), c);
      idx = std::distance(cell_list.begin(), it);
      if(idx < nAdvCellsForEdge(e)) {
        advCoef(e)[idx] = advCoef(e)[idx] + derivTwo(e)[iCell+1][1];
        advCoef_3rd(e)[idx] = advCoef_3rd(e)[idx] - derivTwo(e)[iCell+1][1];
      }
    }

    for(std::size_t iCell = 0; iCell < nAdvCellsForEdge(e); iCell++) {
      advCoef(e)[iCell]     = -(dcEdge(e)*dcEdge(e)) * advCoef(e)[iCell]     / 12.0;
      advCoef_3rd(e)[iCell] = -(dcEdge(e)*dcEdge(e)) * advCoef_3rd(e)[iCell] / 12.0;
    }

    it = std::find(cell_list.begin(), cell_list.end(), cell1);
    idx = std::distance(cell_list.begin(), it);
    if(idx < nAdvCellsForEdge(e)) {
      advCoef(e)[idx] = advCoef(e)[idx] + 0.5;
    }
    
    it = std::find(cell_list.begin(), cell_list.end(), cell2);
    idx = std::distance(cell_list.begin(), it);
    if(idx < nAdvCellsForEdge(e)) {
      advCoef(e)[idx] = advCoef(e)[idx] + 0.5;
    }

    for(std::size_t iCell = 0; iCell < nAdvCellsForEdge(e); iCell++) {    
      advCoef(e)[iCell]     = dvEdge(e) * advCoef(e)[iCell];
      advCoef_3rd(e)[iCell] = dvEdge(e) * advCoef_3rd(e)[iCell];
    }
  }

}

}}}
