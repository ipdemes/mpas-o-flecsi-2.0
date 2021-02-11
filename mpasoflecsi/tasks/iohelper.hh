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

template<mesh::index_space ispace> inline auto it_from_ispace(mesh::accessor<flecsi::ro, flecsi::ro> mesh);

template<>
inline auto
it_from_ispace<mesh::index_space::cells>(mesh::accessor<flecsi::ro, flecsi::ro> mesh) {
  return mesh.cells();
}

template<>
inline auto
it_from_ispace<mesh::index_space::edges>(mesh::accessor<flecsi::ro, flecsi::ro> mesh) {
  return mesh.edges();
}

template<>
inline auto
it_from_ispace<mesh::index_space::vertices>(mesh::accessor<flecsi::ro, flecsi::ro> mesh) {
  return mesh.vertices();
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

}}
