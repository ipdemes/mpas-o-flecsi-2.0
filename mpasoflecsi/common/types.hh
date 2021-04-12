/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <array>

#include "flecsi/data/accessor.hh"

#include "constants.hh"

namespace mpas {

template<class T, int... II>
struct tensor;

template<class T, int... II>
using tensor_t = typename tensor<T, II...>::type;

template<class T>
struct tensor<T> {
  using type = T;
};
template<class T, int I, int... II>
struct tensor<T, I, II...> {
  using type = std::array<tensor_t<T, II...>, I>;
};

template<class T>
using vdtensor = std::array<T, vertexDegree>;

template<class T>
using vltensor = std::array<T, nVertLevels>;
using vlreal = vltensor<double>;

using vltracer = tensor_t<double, nVertLevels, maxTracers>;

template<class T>
using metensor = std::array<T, maxEdges>;

template<class T>
using me2tensor = std::array<T, maxEdges2>;

using derivtensor = tensor_t<double, maxEdges2, 2>;

template <class T, flecsi::partition_privilege_t P1, flecsi::partition_privilege_t P2>
using acc = typename flecsi::template field<T>::template accessor<P1, P2>;

using flecsi::ro;
using flecsi::rw;
using flecsi::na;
using flecsi::wo;

template<class... T>
struct type_list {};

template<class ta, class tb>
struct type_cat;

template<typename... a, typename... b>
struct type_cat<type_list<a...>, type_list<b...>> {
  typedef type_list<a..., b...> type;
};

}
