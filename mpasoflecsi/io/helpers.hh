/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <mpasoflecsi/specialization/mesh.hh>
#include <mpasoflecsi/specialization/definition.hh>

#include "state.hh"
#include "io_hdf5.hh"

namespace mpas { namespace io {

template<class T>
using acc = typename flecsi::field<T>::template accessor<flecsi::wo, flecsi::na>;

template<mesh::index_space ispace> inline auto it_from_ispace(mesh::accessor<flecsi::ro, flecsi::ro> mesh)
{
  if constexpr (ispace == mesh::index_space::cells) return mesh.cells();
  else if constexpr (ispace == mesh::index_space::edges) return mesh.edges();
  else return mesh.vertices();
}

template<mesh::index_space ispace> constexpr inline auto ispace_size_label()
{
  if constexpr (ispace == mesh::index_space::cells) return "nCells";
  else if constexpr (ispace == mesh::index_space::edges) return "nEdges";
  else return "nVertices";
}

template<class T>
class helper1d
{
public:
  helper1d(hid_t infile) : file(infile) {}

  template<mesh::index_space ispace>
  void read(const char * label,
            acc<T> & arr,
            mesh::accessor<flecsi::ro, flecsi::ro> mesh) {
    readds(label, arr, it_from_ispace<ispace>(mesh));
  }

  template<class ITERATOR>
  void readds(const char * label,
              acc<T> & arr,
              ITERATOR it) {
    h5::read_dataset_1D(file, label, temp);
    for (auto ind : it) {
      arr(ind) = temp[ind];
    }
    temp.clear();
  }

  template<mesh::index_space ispace>
  void writeds_dim(const char * label,
                   mpas::acc<T, ro, na> & arr,
                   mesh::accessor<ro, ro> m) {
    using namespace h5;

    auto iter = it_from_ispace<ispace>(m);

    std::size_t gsize = iter.size();

    std::vector<hsize_t> coord;
    temp.resize(iter.size());
    coord.resize(iter.size());
    std::size_t count = 0;
    for(auto ind : iter) {
      temp[count] = arr(ind);
      coord[count] = ind;
      count++;
    }
    mpi_ph5::write_buffer<1>(
      file, label, {gsize}, coord.size(), coord.data(), temp.data());

    attach_scale(file, label, ispace_size_label<ispace>(), 0);

    temp.clear();
  }

  template<mesh::index_space ispace>
  void write_dimscale(mesh::accessor<ro, ro> m) {
    using namespace h5;
    std::size_t gsize = 0;
    std::size_t displ = 0;

    gsize = it_from_ispace<ispace>(m).size();

    auto lsize = gsize;
    constexpr auto label = ispace_size_label<ispace>();

    temp.resize(lsize, 0);
    mpi_ph5::write_buffer<1>(
      file, label, {lsize}, {gsize}, {displ}, temp.data());
    hid_t dset = open_dataset(file, label);
    H5DSset_scale(dset, label);
    close_dataset(dset);
  }

  void write_dimscale(const char * label, std::size_t nelem) {
    using namespace h5;

    hid_t dtype = type_equiv<T>::h5_type();

    temp.resize(nelem, 0);
    hid_t dset = create_dataset<T, 1>(file, {temp.size()}, label);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // if(color == 0)
      H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, plist_id, temp.data());
    // else {
    //   hid_t dset_space = H5Dget_space(dset);
    //   H5Sselect_none(dset_space);
    //   H5Dwrite(dset, dtype, dset_space, dset_space, plist_id, NULL);
    //   H5Sclose(dset_space);
    // }

    H5DSset_scale(dset, label);
    close_dataset(dset);
  }

  template<mesh::index_space ispace>
  void write_timestep(const char * label,
                      mpas::acc<T, ro, na> & arr,
                      mesh::accessor<ro, ro> m,
                      std::size_t timestep)
  {
    using namespace h5;
    auto iter = it_from_ispace<ispace>(m);
    std::size_t gsize = iter.size();
    std::vector<hsize_t> coord;
    temp.resize(iter.size());
    coord.resize(iter.size() * 2);
    std::size_t count = 0;
    for (auto ind : iter) {
      temp[count] = arr(ind);
      coord[2 * count + 0] = timestep;
      coord[2 * count + 1] = ind;
      ++count;
    }

    mpi_ph5::write_buffer<2>(file, label, {inputs::num_output_times.value(), gsize}, temp.size(),
                             coord.data(), temp.data());

    temp.clear();
  }

protected:
  hid_t file;
  std::vector<T> temp;
};


template<class T>
class helper2d
{
public:
  helper2d(hid_t infile) : file(infile), dsreader(h5::read_dataset_2D<T>) {}

  template<mesh::index_space ispace, class ACC>
  void read(const char * label,
            ACC & arr,
            mesh::accessor<flecsi::ro, flecsi::ro> mesh) {
    readds(label, arr, it_from_ispace<ispace>(mesh));
  }

  template<class ACC, class ITERATOR_TYPE>
  void readds(const char * label,
              ACC & arr,
              ITERATOR_TYPE it)
  {
    auto stride = dsreader(file, label, temp);
    for (auto ind : it) {
      for (auto i = 0; i < arr(i).size(); i++) {
        arr(ind)[i] = temp[ind];
      }
    }
    temp.clear();
  }

  template<mesh::index_space ispace, class ACC>
  void write_timestep(const char * label,
                      ACC & arr,
                      mesh::accessor<ro, ro> m,
                      std::size_t timestep)
  {
    auto iter = it_from_ispace<ispace>(m);
    std::size_t gsize = iter.size();
    std::vector<hsize_t> coord;
    temp.resize(iter.size() * arr(0).size());
    coord.resize(iter.size() * arr(0).size() * 3);
    std::size_t count = 0;
    for (auto ind : iter) {
      for (int i = 0; i < arr(i).size(); i++) {
        temp[count] = arr(ind)[i];
        coord[3 * count + 0] = timestep;
        coord[3 * count + 1] = ind;
        coord[3 * count + 2] = i;
        ++count;
      }
    }

    mpi_ph5::write_buffer<3>(file, label, {inputs::num_output_times.value(), gsize, arr(0).size()},
                             temp.size(), coord.data(), temp.data());

    temp.clear();
  }

protected:
  hid_t file;
  std::vector<T> temp;
  std::function<std::size_t(hid_t, const std::string &, std::vector<T> &)> dsreader;
};


template<unsigned short ND, class T>
struct helper;

template<class T>
struct helper<1, T> { using type = helper1d<T>; };
template<class T>
struct helper<2, T> { using type = helper2d<T>; };


}}
