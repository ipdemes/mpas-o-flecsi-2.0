/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <mpasoflecsi/specialization/mesh.hh>
#include <mpasoflecsi/specialization/definition.hh>

namespace mpas { namespace io {

template<class T>
using acc = typename flecsi::field<T>::template accessor<flecsi::wo, flecsi::na>;

template<mesh::index_space ispace> inline auto it_from_ispace(mesh::accessor<flecsi::ro, flecsi::ro> mesh)
{
  if constexpr (ispace == mesh::index_space::cells) return mesh.cells();
  else if constexpr (ispace == mesh::index_space::edges) return mesh.edges();
  else return mesh.vertices();
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

protected:
  hid_t file;
  std::vector<T> temp;
  std::function<std::size_t(hid_t, const std::string &, std::vector<T> &)> dsreader;
};

}}
