/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <flecsi/execution.hh>

#include "mpasoflecsi/util/part.hh"

namespace mpas { namespace io {
namespace inputs {

inline flecsi::program_option<std::size_t>
num_output_times("Output", "num_output", "Number of timesteps to output",
                 {{flecsi::option_default, 0}});

}

class output_flagger
{
public:
  output_flagger(int num_timesteps) :
    part(num_timesteps, inputs::num_output_times.value() - 1) {}

  inline operator bool() const {
    return inputs::num_output_times.value() > 0;
  }

  inline bool output_step(std::size_t i) const {
    return part.high(part.owner(i)) == i;
  }

protected:
  util::block_partition part;
};

}}
