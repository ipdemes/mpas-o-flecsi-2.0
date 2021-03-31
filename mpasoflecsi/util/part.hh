/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

namespace mpas { namespace util {

class partition
{
public:
  using len_t = std::size_t;

partition(len_t n, len_t num_procs) : num_procs(num_procs), n(n) {}
	virtual ~partition() {};
	virtual len_t low(len_t id) const = 0;
	virtual len_t high(len_t id) const = 0;
	virtual len_t size(len_t id) const = 0;
	virtual len_t owner(len_t index) const  = 0;

protected:
	len_t num_procs;
	len_t n;
};


class block_partition: public partition
{
public:
block_partition(len_t n, len_t num_procs):
	partition(n, num_procs) {};

	inline len_t low(len_t id) const override
	{
		return id * n / num_procs;
	};

	inline len_t high(len_t id) const override
	{
		return low(id+1) - 1;
	};

	inline len_t size(len_t id) const override
	{
		return high(id) - low(id) + 1;
	};

	inline len_t owner(len_t index) const override
	{
		return (num_procs*(index+1)-1)/n;
	};
};

}}
