/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "io.hh"

namespace mpas { namespace task {

using namespace flecsi;

void read_mesh_fields(mesh::accessor<ro, ro> m,
                      hid_t file,
                      io::acc<double> latCell,
                      io::acc<double> lonCell,
                      io::acc<double> latVertex,
                      io::acc<double> lonVertex,
                      io::acc<double> latEdge,
                      io::acc<double> lonEdge,
                      io::acc<double> xCell,
                      io::acc<double> yCell,
                      io::acc<double> zCell,
                      io::acc<double> areaCell,
                      io::acc<double> areaTriangle,
                      io::acc<vdtensor<int>> kiteAreasOnVertex,
                      io::acc<double> dvEdge,
                      io::acc<double> dcEdge,
                      io::acc<double> meshDensity,
                      io::acc<double> fVertex,
                      io::acc<metensor<double>> weightsOnEdge)
{
  using namespace mpas::io;

  constexpr auto cell_space = mesh::index_space::cells;
  constexpr auto edge_space = mesh::index_space::edges;
  constexpr auto vert_space = mesh::index_space::vertices;
  { // read 1D float datasets
    helper1d<double> reader(file);
    reader.read<vert_space>("fVertex", fVertex, m);
    reader.read<cell_space>("areaCell", areaCell, m);
    reader.read<vert_space>("areaTriangle", areaTriangle, m);
    reader.read<edge_space>("dvEdge", dvEdge, m);
    reader.read<edge_space>("dcEdge", dcEdge, m);
    reader.read<cell_space>("meshDensity", meshDensity, m);
    reader.read<cell_space>("xCell", xCell, m);
    reader.read<cell_space>("yCell", yCell, m);
    reader.read<cell_space>("zCell", zCell, m);
    reader.read<cell_space>("latCell", latCell, m);
    reader.read<cell_space>("lonCell", lonCell, m);
    reader.read<vert_space>("latVertex", latVertex, m);
    reader.read<vert_space>("lonVertex", lonVertex, m);
    reader.read<vert_space>("latEdge", latEdge, m);
    reader.read<vert_space>("lonEdge", lonEdge, m);
  }

  { // read 2D float datasets
    helper2d<double> reader(file);
    reader.read<edge_space>("weightsOnEdge", weightsOnEdge, m);
    reader.read<vert_space>("kiteAreasOnVertex", kiteAreasOnVertex, m);
  }
}

}}
