/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "mpasoflecsi/io/io_hdf5.hh"
#include "mpasoflecsi/io/state.hh"
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
                      io::acc<vdtensor<double>> kiteAreasOnVertex,
                      io::acc<double> dvEdge,
                      io::acc<double> dcEdge,
                      io::acc<double> meshDensity,
                      io::acc<double> fVertex,
                      io::acc<me2tensor<double>> weightsOnEdge)
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
    reader.read<edge_space>("latEdge", latEdge, m);
    reader.read<edge_space>("lonEdge", lonEdge, m);
  }

  { // read 2D float datasets
    helper2d<double> reader(file);
    reader.read<edge_space>("weightsOnEdge", weightsOnEdge, m);
    reader.read<vert_space>("kiteAreasOnVertex", kiteAreasOnVertex, m);
  }
}


void mesh_output_init(mesh::accessor<ro, ro> m,
                      acc<double, ro, na> areaCell,
                      acc<double, ro, na> dvEdge,
                      acc<double, ro, na> dcEdge,
                      acc<double, ro, na> xCell,
                      acc<double, ro, na> yCell,
                      acc<double, ro, na> zCell,
                      acc<double, ro, na> latCell,
                      acc<double, ro, na> lonCell,
                      acc<double, ro, na> latVertex,
                      acc<double, ro, na> lonVertex,
                      acc<double, ro, na> meshDensity,
                      const char * mfile_name,
                      hid_t * mfile)
{
  using namespace mpas::io;

  constexpr auto cell_space = mesh::index_space::cells;
  constexpr auto edge_space = mesh::index_space::edges;
  constexpr auto vert_space = mesh::index_space::vertices;

  hid_t file = mpi_ph5::create_file(mfile_name);

  { // write attributes
    auto root = H5Gopen(file, "/", H5P_DEFAULT);
    hid_t dsid = H5Screate(H5S_SCALAR);
    auto atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, 3);
    H5Tset_strpad(atype, H5T_STR_NULLTERM);
    auto attrid =
      H5Acreate2(root, "on_a_sphere", atype, dsid, H5P_DEFAULT, H5P_DEFAULT);
    auto ret = H5Awrite(attrid, atype, "YES");
    ret = H5Sclose(dsid);
    ret = H5Aclose(attrid);
    ret = H5Tclose(atype);
    ret = H5Gclose(root);
  } // end write attributes

  { // Write dimension scales
    helper1d<int> writer(file);
    writer.write_dimscale<cell_space>(m);
    writer.write_dimscale<edge_space>(m);
    writer.write_dimscale<vert_space>(m);
    writer.write_dimscale("vertexDegree", vertexDegree);
    writer.write_dimscale("nVertLevels", nVertLevels);
    writer.write_dimscale("Time", io::inputs::num_output_times.value());
    writer.write_dimscale("nTracers", maxTracers);
  } // end write dimension scales
  { // Write 1D float arrays
    helper1d<double> writer(file);
    writer.writeds_dim<cell_space>("areaCell", areaCell, m);
    writer.writeds_dim<cell_space>("xCell", xCell, m);
    writer.writeds_dim<cell_space>("yCell", yCell, m);
    writer.writeds_dim<cell_space>("zCell", zCell, m);
    writer.writeds_dim<cell_space>("latCell", latCell, m);
    writer.writeds_dim<cell_space>("lonCell", lonCell, m);
    writer.writeds_dim<vert_space>("latVertex", latVertex, m);
    writer.writeds_dim<vert_space>("lonVertex", lonVertex, m);
    writer.writeds_dim<cell_space>("mDensity", meshDensity, m);
    writer.writeds_dim<edge_space>("dvEdge", dvEdge, m);
    writer.writeds_dim<edge_space>("dcEdge", dcEdge, m);
  } // End 2D float arrays
  { // Write 2D int arrays
    std::vector<int> tmp;
    for (auto v : m.vertices()) {
      for (auto c : m.cells(v)) {
        tmp.push_back(c + 1);
      }
    }
    std::size_t gsize = tmp.size();
    std::vector<hsize_t> coord;
    coord.resize(tmp.size() * 2);
    std::size_t count{0};
    for (auto v : m.vertices()) {
      std::size_t i{0};
      for (auto c : m.cells(v)) {
        coord[2*count] = v;
        coord[2*count + 1] = i;
        ++i; ++count;
      }
    }
    mpi_ph5::write_buffer<2>(file, "cellsOnVertex", {gsize, 3}, tmp.size(), coord.data(), tmp.data());
    h5::attach_scale(file, "cellsOnVertex", "nVertices", 0);
    h5::attach_scale(file, "cellsOnVertex", "vertexDegree", 1);
  } // End 2D int arrays

  *mfile = file;
}


void mesh_output_finalize(hid_t mfile)
{
  io::mpi_ph5::close_file(mfile);

}

}}
