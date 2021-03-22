/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#pragma once

/// \file

#include <mpi.h>

extern "C" {
#include <hdf5.h>
#include <hdf5_hl.h>
}

// system includes
#include <algorithm>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "flecsi/flog.hh"
#include "flecsi/data.hh"


namespace mpas {
namespace io {

namespace detail {
  struct crs : flecsi::topo::unstructured_impl::crs {
    std::size_t size() const {
      if (offsets.empty())
        return 0;
      else
        return offsets.size() - 1;
    }

    void clear() {
      offsets.clear();
      indices.clear();
    }
  };
} // namespace detail

/* Notes
 *
 * The use of HDF5 here is fairly primitive. Any API failure will cause a
 * program abort, but there are likely ways to recover gracefully in many cases.
 *
 * It might be worth trying to manage what resources are open and need to be
 * closed with good C++ idioms, but I don't know if that really makes sense in a
 * multi-node context.  I'll have to answer that once I've got a better handle
 * on the parallel use case.
 *
 * It's worth considering that the HDF5 C++ API does seems to do good job of
 * resource management, but is not considered thread safe, occasionally making
 * multiple calls into the (optionally threadsafe) C API at a time.
 *
 *
 * # Datatypes and conversions
 *
 * The scalar type of the data in any HDF5 dataset is fixed, but when calling
 * H5DRead, the datatype argument does not need to match this.  The library will
 * attempt to convert types appropriately so long as the given datatype and the
 * "real" datatype are in the same type class.  Some types cannot be converted
 * without e.g. loss of precision or overflow, in which case I *think* an error
 * is signaled, and (if properly registered) a callback function will be called.
 *
 */

/*!
 * \brief Simple wrappers for some HDF5 functions that abort
 * on error
 */
namespace h5 {

/* Static inferrence of the appropriate HDF5 datatype for
 * type T */
template<typename T>
struct type_equiv {};

#define make_type_equiv(cxx_t, h5_t)                                           \
  template<>                                                                   \
  struct type_equiv<cxx_t> {                                                   \
    static hid_t h5_type() {                                                   \
      return h5_t;                                                             \
    }                                                                          \
  }

/* HDF5 native datatypes are aliases for width-specific types set depending on
 * architecture.  It's likely better to use the more precise types if the HDF5
 * schema is known. */
make_type_equiv(char, H5T_NATIVE_CHAR);
make_type_equiv(signed char, H5T_NATIVE_SCHAR);
make_type_equiv(unsigned char, H5T_NATIVE_UCHAR);
make_type_equiv(short, H5T_NATIVE_SHORT);
make_type_equiv(unsigned short, H5T_NATIVE_USHORT);
make_type_equiv(int, H5T_NATIVE_INT);
make_type_equiv(unsigned, H5T_NATIVE_UINT);
make_type_equiv(long, H5T_NATIVE_LONG);
make_type_equiv(unsigned long, H5T_NATIVE_ULONG);
make_type_equiv(long long, H5T_NATIVE_LLONG);
make_type_equiv(unsigned long long, H5T_NATIVE_ULLONG);
make_type_equiv(float, H5T_NATIVE_FLOAT);
make_type_equiv(double, H5T_NATIVE_DOUBLE);
make_type_equiv(long double, H5T_NATIVE_LDOUBLE);

#undef make_type_equiv

/*!
 * \brief Callback function for HDF5 datatype conversion
 * errors
 *
 * This currently just warns the client about the type of
 * conversion error, and aborts the type conversion.
 *
 * See
 * https://portal.hdfgroup.org/display/HDF5/H5P_SET_TYPE_CONV_CB
 */
inline H5T_conv_ret_t
handle_conversion_err(H5T_conv_except_t except_type,
  hid_t src_type_id,
  hid_t dst_type_id,
  void * src_buf,
  void * dst_buf,
  void * op_data) {
  // TODO: Maybe try to recover from the error and perform a
  // manual conversion? It seems likely that if conversion
  // errors happen, there's something else that needs to be
  // fixed (e.g. the code reading datasets does not match
  // the actual schema), so it's probably ok to just bail
  // out.
  switch(except_type) {
    case H5T_CONV_EXCEPT_RANGE_HI:
      flog(warn) << "Conversion failure: Overflow\n"
                 << "Source value is positive and "
                 << "its magnitude is too big for the destination." << std::endl;
      break;

    case H5T_CONV_EXCEPT_RANGE_LOW:
      flog(warn) << "Conversion failure: Overflow\n"
                 << "Source value is negative and "
                 << "its magnitude is too big for the destination." << std::endl;
      break;

    case H5T_CONV_EXCEPT_TRUNCATE:
      flog(warn) << "Conversion failure: Truncation\n"
                 << "Source is floating-point "
                 << "type and destination is integer.\n The floating-point number "
                 << "has fractional part." << std::endl;
      break;

    case H5T_CONV_EXCEPT_PRECISION:
      flog(warn) << "Conversion failure: Precision\n"
                 << "Source is integer and destination is floating-point type. The "
                 << "mantissa of floating-point type is not big enough to hold all "
                 << "the digits of the integer." << std::endl;
      break;

    case H5T_CONV_EXCEPT_PINF:
      flog(warn) << "Conversion failure: Infinity\n"
                 << "Source is floating-point type "
                 <<"and the value is positive infinity." << std::endl;
      break;

    case H5T_CONV_EXCEPT_NINF:
      flog(warn) << "Conversion failure: Infinity\n"
                 << "Source is floating-point type "
                 << "and the value is negative infinity." << std::endl;
      break;

    case H5T_CONV_EXCEPT_NAN:
      flog(warn) << "Conversion failure: NAN\n"
                 << "Source is floating-point type and "
                 << "the value is NaN (not a number, including QNaN and SNaN).";
      break;
    default:; // should be unreachable - by API guarantee
  }

  // For now, just do the safe thing and always signal an
  // aborted conversion. This should print an hdf5 stack
  // trace before returning to where the conversion was
  // intiated.  Note that it's explicitly advised against
  // using C++ exceptions here because they will bypass
  // cleanup routines in the HDF5 library.
  return H5T_CONV_ABORT;
}

inline hid_t
create_file(const std::string name, unsigned flags, hid_t fapl_id) {
  hid_t file = H5Fcreate(name.c_str(), flags, fapl_id, H5P_DEFAULT);
  if(file < 0) {
    flog_fatal("could not open file: " << name);
  }
  return file;
}


/*!
 * \brief Open an HDF5 file with the given flags and access properties
 *
 * This will fail hard if there are problems opening the file (e.g. not
 * readable/writeable, does not exist, is not an HDF5 file)
 *
 * \param[in] name file name to open
 *
 * \param[in] flags read/write/create flags to open the file with
 *
 * \param[in] fapl_id handle to the file access property list to use
 *
 * \return a handle to the open file
 */
inline hid_t
open_file(const std::string & name, unsigned flags, hid_t fapl_id) {
  // Note: Third param (fapl_id) is the id for file access props For parallel
  // access, the fapl_id will hold the communicator
  hid_t file = H5Fopen(name.c_str(), flags, H5P_DEFAULT);

  if(file < 0) {
    // Try to give some meaningful feedback if open fails
    htri_t res = H5Fis_hdf5(name.c_str());

    if(res == 0) { // res == 0 ==> exists but not in HDF5 format
      flog_fatal(name << " is not an HDF5 file");
    } // if
    else {
      // TODO: At least check if the file exists
      flog_fatal("Couldn't open " << name << " ... not sure why");
    } // else
  } // if
  return file;
} // open_file

/*!
 * \brief wrapper for H5Fclose
 *
 * Note: API-level failure causes program abort.  What effect this may have on
 * the file being closed is unclear.
 *
 * \param[in] file Handle to an open file
 */
inline void
close_file(hid_t file_id) {
  /* If there are still objects from the given file open, the actual closing of
   * this file is delayed until all of those objects are closed as well.  This
   * is a problem when using the MPI-IO driver, so care must be taken to close
   * all objects before closing the file in an MPI context.  See
   * https://portal.hdfgroup.org/display/HDF5/H5F_CLOSE. */
  if(H5Fclose(file_id) < 0) {
    flog_fatal("Closing file with handle " << file_id << " failed!");
  } // if
} // close_file

/*!
 * \brief Open the specified dataset, failing hard if the open is unsuccessful
 *
 * \param[in] loc the handle to an open location identifier
 *
 * \param[in] name name of the dataset to access
 *
 * \return a handle to the dataset
 */
inline hid_t
open_dataset(hid_t loc, const std::string & name) {
  hid_t handle = H5Dopen(loc, name.c_str(), H5P_DEFAULT);
  flog_assert(handle >= 0, "Failed to open dataset: " << name);
  return handle;
} // open_dataset


/*!
 * \brief Create the specified dataset, failing hard if the open is unsuccessful
 *
 * \tparam T C++ type of dataset
 * \tparam ND number of dimensions
 * \param[in] loc the handle to an open location identifier
 * \param[in] dims dimensions of the dataset
 * \param[in] name name of the dataset to create
 * \return a handle to the dataset
 */
template<class T, unsigned short ND>
inline hid_t
create_dataset(hid_t loc,
  std::array<hsize_t, ND> dims,
  const std::string & name) {
  flog_assert((ND == 1) or (ND == 2),
    "Not implemented: create_dataset for dimension = " << ND);
  auto dataspace_id = H5Screate_simple(ND, dims.data(), NULL);
  auto dtype = type_equiv<T>::h5_type();
  hid_t dataset_id;
  if(ND == 1) {
    dataset_id =
      H5Dcreate1(loc, name.c_str(), dtype, dataspace_id, H5P_DEFAULT);
  }
  else {
    dataset_id = H5Dcreate2(loc, name.c_str(), dtype, dataspace_id, H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);
  }
  flog_assert(dataset_id >= 0, "Failed to create dataset: " << name);

  return dataset_id;
} // create_dataset

/*!
 * \brief wrapper for H5Dclose
 *
 * Note: API-level failure causes program abort.  What effect this may have on
 * the underlying file is unclear.
 *
 * \param[in] dset_id Handle to an open dataset
 */
inline void
close_dataset(hid_t dset_id) {
  if(H5Dclose(dset_id) < 0) {
    flog_fatal("Closing dataset with handle " << dset_id << " failed!");
  } // if
} // close_dataset

/*!
 * \brief Get a handle to a copy of a dataset's dataspace, failing hard if the
 * open is unsuccessful
 *
 * Note that callers are responsible for closing the copied dataset using
 * close_dataspace (or H5Sclose)
 *
 * \param[in] dset_id handle to the open dataspace
 *
 * \return a handle to the copied dataspace
 */
inline hid_t
get_space(hid_t dset_id) {
  hid_t handle = H5Dget_space(dset_id);

  // TODO: use H5Iget_name() to give a more meaningful error
  // message
  flog_assert(handle >= 0, "Failed to make a copy of dataspace");

  return handle;
} // get_space

/*!
 * \brief wrapper for H5Sclose
 *
 * Note: API-level failure causes program abort.  What effect this may have on
 * the underlying file is unclear.
 *
 * \param[in] dspace_id Handle to an open dataspace
 */
inline void
close_dataspace(hid_t dspace_id) {
  if(H5Sclose(dspace_id) < 0) {
    flog_fatal("Closing dataspace with handle " << dspace_id << " failed!");
  } // if
} // close_dataspace


/*!
 * \brief fails hard if the given dataset does not hold the expected class
 *
 * \param[in] dset handle to the dataset
 *
 * \param[in] expect_class expected class of data
 */
inline void
assert_dataset_typeclass(hid_t dset, H5T_class_t expect_class) {
  hid_t dtype_id = H5Dget_type(dset);
  H5T_class_t real_class = H5Tget_class(dtype_id);
  flog_assert(real_class == expect_class,
    "Unexpected dataset type class " << real_class << " != " << expect_class);
  H5Tclose(dtype_id);
}

/*!
 * \brief fails hard if the given dataset does not hold the expected type
 *
 * Example usage: assert_dataset_type(dset_id, H5T_NATIVE_INT)
 *
 * \param[in] dset handle to the dataset
 *
 * \param[in] expect_type handle to hid_t of expected type of data
 */
inline void
assert_dataset_type(hid_t dset, hid_t expect_type) {
  hid_t dtype_id = H5Dget_type(dset);

  if(!H5Tequal(dtype_id, expect_type)) {
    const size_t max = 101;
    char real_name[max - 1], expect_name[max - 1], dset_name[max - 1];

    // TODO: Even for well-known types, these calls will sometimes (often?)
    // fail.  If debugability of this sort of thing becomes very important,
    // it'll likely be worth wrapping this function with one that can infer and
    // stringify names of common types (e.g. H5T_NATIVE_DOUBLE)
    H5Iget_name(dset, dset_name, max);
    H5Iget_name(dtype_id, real_name, max);
    H5Iget_name(expect_type, expect_name, max);
    flog_fatal("Dataset " << dset_name << " has type: " << real_name
               << " but I'm expecting " << expect_name);
  }
  H5Tclose(dtype_id);
} // assert_dataset_type

/*! \brief Get the rank of a dataset, failing hard if the HDF5 API call is
 * unsuccessful
 *
 *  \param[in] dset_id handle to the dataset
 *
 *  \return the rank of the dataset (number of dimensions)
 */
inline int
get_rank(hid_t dset_id) {
  hid_t space = get_space(dset_id);
  int rank = H5Sget_simple_extent_ndims(space);

  // TODO: use H5Iget_name() to give a more meaningful error message
  flog_assert(rank >= 0, "Failed to retrieve dataset rank");

  H5Sclose(space);
  return rank;
} // get_rank

/*!
 * \brief Wrapper for H5Sget_simple_extent_dims to retrieve dimension
 * information of a dataset
 *
 * The caller must allocate appropriate space to hold the dimension extents
 *
 * \param[in] dset_id dataset handle
 *
 * \param[out] dims pointer to array
 *
 * \return the number of dimesions (rank) of the dataset
 */
inline int
get_simple_dims(hid_t dset_id, hsize_t * dims) {
  hid_t space = get_space(dset_id);
  int rank = H5Sget_simple_extent_dims(space, dims, NULL);
  H5Sclose(space);
  return rank;
} // get_simple_dims

/*!
 * \brief Wrapper for H5Pcreate
 *
 * \param[in] plist_type Identifier of the type of property list to create
 *
 */

inline hid_t
create_plist(hid_t plist_type) {
  hid_t plist = H5Pcreate(plist_type);
  flog_assert(plist >= 0, "Failed to create property list");
  return plist;
}

/* !
 * \brief read data from a dataset in a given *open* file
 *
 * \param[in] file_handle to open HDF5 file
 *
 * \param[in] name name of the dataset to read from
 *
 * \param[out] data vector to add coordinates to
 *
 * \return number of entities read  */
template<typename T>
inline size_t
read_dataset_1D(hid_t file_handle,
  const std::string & dset_name,
  std::vector<T> & data) {

  hid_t transfer_plist = create_plist(H5P_DATASET_XFER);

  if(H5Pset_type_conv_cb(transfer_plist, handle_conversion_err, nullptr) < 0) {
    flog_fatal("Failed to register type conversion error callback");
  }

  hid_t dset = open_dataset(file_handle, dset_name);
  hid_t dtype = type_equiv<T>::h5_type();

  /* Check that the type we're trying to read shares the same typeclass as
   * what's actually in the dataset. */
  H5T_class_t typeclass = H5Tget_class(dtype);
  assert_dataset_typeclass(dset, typeclass);

  // Expecting two-dimensional data
  flog_assert(
    get_rank(dset) == 1, "Expected two dimensions for coordinate data");

  hsize_t dset_size;
  get_simple_dims(dset, &dset_size);

  T * buf = new T[dset_size];
  herr_t status = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, transfer_plist, buf);
  flog_assert(status >= 0, "Failed to read dataset \"" << dset_name << "\"");

  data.insert(data.end(), buf, &buf[dset_size]);

  delete[] buf;
  return dset_size;
}

/* !
 * \brief read data from a dataset in a given *open* file
 *
 * \param[in] file_handle hdf5 file handle for file to read
 * \param[in] dset_name name of dataset to read
 * \param[in] data flat vector to read data to (stride for first dimension
 * returned) \return stride for first dimension
 */
template<typename T>
inline size_t
read_dataset_2D(hid_t file_handle,
  const std::string & dset_name,
  std::vector<T> & data) {

  hid_t transfer_plist = create_plist(H5P_DATASET_XFER);

  if(H5Pset_type_conv_cb(transfer_plist, handle_conversion_err, nullptr) < 0) {
    flog_fatal("Failed to register type conversion error callback");
  }

  hid_t dset = open_dataset(file_handle, dset_name);
  hid_t dtype = type_equiv<T>::h5_type();

  /* Check that the type we're trying to read shares the same typeclass as
   * what's actually in the dataset. */
  H5T_class_t typeclass = H5Tget_class(dtype);
  assert_dataset_typeclass(dset, typeclass);

  // Expecting two-dimensional data
  flog_assert(
    get_rank(dset) == 2, "Expected two dimensions for coordinate data");

  hsize_t dset_sizes[2];
  get_simple_dims(dset, dset_sizes);

  size_t dim1_size = dset_sizes[0];
  size_t dim2_size = dset_sizes[1];

  data.resize(dim1_size * dim2_size);

  herr_t status =
    H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
  flog_assert(status >= 0, "Failed to read dataset \"" << dset_name << "\"");

  return dim2_size;
}

/* !
 * \brief 2D connectivity data from a dataset in a given *open* file
 *
 *
 */
template<typename T>
inline size_t
read_connectivity_2D(hid_t file_handle,
                     const std::string & dset_name,
                     detail::crs & data,
                     std::vector<std::size_t> & local_to_global) {

  hid_t transfer_plist = create_plist(H5P_DATASET_XFER);

  if(H5Pset_type_conv_cb(transfer_plist, handle_conversion_err, nullptr) < 0) {
    flog_fatal("Failed to register type conversion error callback");
  }

  hid_t dset = open_dataset(file_handle, dset_name);
  hid_t dtype = type_equiv<T>::h5_type();

  /* Check that the type we're trying to read shares the same typeclass as
   * what's actually in the dataset. */
  H5T_class_t typeclass = H5Tget_class(dtype);
  assert_dataset_typeclass(dset, typeclass);

  // Expecting two-dimensional data
  flog_assert(
    get_rank(dset) == 2, "Expected two dimensions for coordinate data");

  hsize_t dset_sizes[2];
  get_simple_dims(dset, dset_sizes);

  size_t dim1_size = dset_sizes[0];
  size_t dim2_size = dset_sizes[1];

  T * buf = new T[dim1_size * dim2_size];

  herr_t status = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  flog_assert(status >= 0, "Failed to read dataset \"" << dset_name << "\"");

  data.clear();
  data.offsets.emplace_back(0);
  size_t counter{0};

  for(size_t gid : local_to_global) {
    for(size_t j = 0; j < dim2_size; j++) {
      if(buf[gid * dim2_size + j] != 0) {
        data.indices.push_back(buf[gid * dim2_size + j] - 1);
        counter++;
      }
    }
    data.offsets.push_back(counter);
  }

  delete[] buf;

  return dim1_size * dim2_size;
}

template<typename T>
inline size_t
read_connectivity_2D(hid_t file_handle,
                     const std::string & dset_name,
                     detail::crs & data,
                     size_t start,
                     size_t end) {

  hid_t transfer_plist = create_plist(H5P_DATASET_XFER);

  if(H5Pset_type_conv_cb(transfer_plist, handle_conversion_err, nullptr) < 0) {
    flog_fatal("Failed to register type conversion error callback");
  }

  hid_t dset = open_dataset(file_handle, dset_name);
  hid_t dtype = type_equiv<T>::h5_type();

  /* Check that the type we're trying to read shares the same typeclass as
   * what's actually in the dataset. */
  H5T_class_t typeclass = H5Tget_class(dtype);
  assert_dataset_typeclass(dset, typeclass);

  // Expecting two-dimensional data
  flog_assert(
    get_rank(dset) == 2, "Expected two dimensions for coordinate data");

  hsize_t dset_sizes[2];
  get_simple_dims(dset, dset_sizes);

  size_t dim1_size = dset_sizes[0];
  size_t dim2_size = dset_sizes[1];

  T * buf = new T[dim1_size * dim2_size];

  herr_t status = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  flog_assert(status >= 0, "Failed to read dataset \"" << dset_name << "\"");

  data.clear();
  data.offsets.emplace_back(0);
  size_t counter{0};

  for(size_t i = start; i <= end; i++) {
    for(size_t j = 0; j < dim2_size; j++) {
      if(buf[i * dim2_size + j] != 0) {
        data.indices.push_back(buf[i * dim2_size + j] - 1);
        counter++;
      }
    }
    data.offsets.push_back(counter);
  }

  delete[] buf;
  return dim1_size * dim2_size;
}

template<typename T>
inline size_t
read_connectivity_2D(hid_t file_handle,
                     const std::string & dset_name,
                     detail::crs & data) {
  hid_t transfer_plist = create_plist(H5P_DATASET_XFER);

  if(H5Pset_type_conv_cb(transfer_plist, handle_conversion_err, nullptr) < 0) {
    flog_fatal("Failed to register type conversion error callback");
  }

  hid_t dset = open_dataset(file_handle, dset_name);
  hid_t dtype = type_equiv<T>::h5_type();

  /* Check that the type we're trying to read shares the same typeclass as
   * what's actually in the dataset. */
  H5T_class_t typeclass = H5Tget_class(dtype);
  assert_dataset_typeclass(dset, typeclass);

  // Expecting two-dimensional data
  flog_assert(
    get_rank(dset) == 2, "Expected two dimensions for coordinate data");

  hsize_t dset_sizes[2];
  get_simple_dims(dset, dset_sizes);

  size_t dim1_size = dset_sizes[0];
  size_t dim2_size = dset_sizes[1];

  T * buf = new T[dim1_size * dim2_size];

  herr_t status = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  flog_assert(status >= 0, "Failed to read dataset \"" << dset_name << "\"");

  data.clear();
  data.offsets.emplace_back(0);
  size_t counter{0};
  for(size_t i = 0; i < dim1_size; i++) {
    for(size_t j = 0; j < dim2_size; j++) {
      if(buf[i * dim2_size + j] != 0) {
        data.indices.push_back(buf[i * dim2_size + j] - 1);
        counter++;
      }
    }
    data.offsets.push_back(counter);
  }
  delete[] buf;
  return dim1_size * dim2_size;
}

template<typename T>
inline size_t
read_entity_to_entity(hid_t file_handle,
                      const std::string & dset_name,
                      detail::crs & data,
                      std::vector<std::size_t> & local_to_global) {

  hid_t transfer_plist = create_plist(H5P_DATASET_XFER);

  if(H5Pset_type_conv_cb(transfer_plist, handle_conversion_err, nullptr) < 0) {
    flog_fatal("Failed to register type conversion error callback");
  }

  hid_t dset = open_dataset(file_handle, dset_name);
  hid_t dtype = type_equiv<T>::h5_type();

  /* Check that the type we're trying to read shares the same typeclass as
   * what's actually in the dataset. */
  H5T_class_t typeclass = H5Tget_class(dtype);
  assert_dataset_typeclass(dset, typeclass);

  // Expecting two-dimensional data
  flog_assert(
    get_rank(dset) == 2, "Expected two dimensions for coordinate data");

  hsize_t dset_sizes[2];
  get_simple_dims(dset, dset_sizes);

  size_t dim1_size = dset_sizes[0];
  size_t dim2_size = dset_sizes[1];

  T * buf = new T[dim1_size * dim2_size];

  herr_t status = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  flog_assert(status >= 0, "Failed to read dataset \"" << dset_name << "\"");

  data.clear();
  data.offsets.emplace_back(0);
  size_t counter{0};

  for(size_t gid : local_to_global) {
    for(size_t j = 0; j < dim2_size; j++) {
      if(buf[gid * dim2_size + j] != 0 &&
         buf[gid * dim2_size + j] <= dim1_size) {
        data.indices.push_back(buf[gid * dim2_size + j] - 1);
        counter++;
      }
    }
    data.offsets.push_back(counter);
  }

  delete[] buf;
  return dim1_size * dim2_size;
}

template<typename T>
inline size_t
read_entity_to_entity(hid_t file_handle,
                      const std::string & dset_name,
                      detail::crs & data) {

  hid_t transfer_plist = create_plist(H5P_DATASET_XFER);

  if(H5Pset_type_conv_cb(transfer_plist, handle_conversion_err, nullptr) < 0) {
    flog_fatal("Failed to register type conversion error callback");
  }

  hid_t dset = open_dataset(file_handle, dset_name);
  hid_t dtype = type_equiv<T>::h5_type();

  /* Check that the type we're trying to read shares the same typeclass as
   * what's actually in the dataset. */
  H5T_class_t typeclass = H5Tget_class(dtype);
  assert_dataset_typeclass(dset, typeclass);

  // Expecting two-dimensional data
  flog_assert(
    get_rank(dset) == 2, "Expected two dimensions for coordinate data");

  hsize_t dset_sizes[2];
  get_simple_dims(dset, dset_sizes);

  size_t dim1_size = dset_sizes[0];
  size_t dim2_size = dset_sizes[1];

  T * buf = new T[dim1_size * dim2_size];

  herr_t status = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  flog_assert(status >= 0, "Failed to read dataset \"" << dset_name << "\"");

  data.clear();
  data.offsets.emplace_back(0);
  size_t counter{0};

  for(size_t gid = 0; gid < dim1_size; gid++) {
    for(size_t j = 0; j < dim2_size; j++) {
      if(buf[gid * dim2_size + j] != 0 &&
         buf[gid * dim2_size + j] <= dim1_size) {
        data.indices.push_back(buf[gid * dim2_size + j] - 1);
        counter++;
      }
    }
    data.offsets.push_back(counter);
  }
  delete[] buf;
  return dim1_size * dim2_size;
}

/* !
 * \brief read data from a dataset in a given *open* file
 *
 * \param[in] file_handle hdf5 file handle for file to read
 * \param[in] dset_name name of dataset to read
 * \param[in] data flat vector to read data to (stride for first dimension
 * returned) \return stride for first dimension
 */
template<typename T>
inline size_t
read_connectivity_v_2D(hid_t file_handle,
  const std::string & dset_name,
  std::vector<T> & data) {

  hid_t transfer_plist = create_plist(H5P_DATASET_XFER);

  if(H5Pset_type_conv_cb(transfer_plist, handle_conversion_err, nullptr) < 0) {
    flog_fatal("Failed to register type conversion error callback");
  }

  hid_t dset = open_dataset(file_handle, dset_name);
  hid_t dtype = type_equiv<T>::h5_type();

  /* Check that the type we're trying to read shares the same typeclass as
   * what's actually in the dataset. */
  H5T_class_t typeclass = H5Tget_class(dtype);
  assert_dataset_typeclass(dset, typeclass);

  // Expecting two-dimensional data
  flog_assert(
    get_rank(dset) == 2, "Expected two dimensions for coordinate data");

  hsize_t dset_sizes[2];
  get_simple_dims(dset, dset_sizes);

  size_t dim1_size = dset_sizes[0];
  size_t dim2_size = dset_sizes[1];

  data.resize(dim1_size * dim2_size);

  herr_t status =
    H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
  // MPAS meshes use 1 based indexing, flecsi wants 0 based for connectivity
  for(auto & v : data)
    --v;

  flog_assert(status >= 0, "Failed to read dataset \"" << dset_name << "\"");

  return dim2_size;
}

/**
 * \brief Write 1D dataset to a given *open* file
 *
 * \param[in] file_handle file to write dataset
 * \param[in] dset_name name of dataset to write
 * \param[in] data 1D data vector to write to dataset
 */
template<class T>
void
write_dataset_1D(hid_t file_handle,
  const std::string & dset_name,
  std::vector<T> & data) {
  hid_t dset = create_dataset<T, 1>(file_handle, {data.size()}, dset_name);
  hid_t dtype = type_equiv<T>::h5_type();
  herr_t status =
    H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
  // clog_assert(status >= 0, "Failed to write dataset \"" << dset_name << "\"");
  status = H5Dclose(dset);
  // clog_assert(status >= 0, "Failed to close dataset \"" << dset_name << "\"");
}

/**
 * \brief Write 2D dataset to a given *open* file
 *
 * \param[in] file_handle file to write dataset
 * \param[in] dset_name name of dataset to write
 * \param[in] dim shape of array to write
 * \param[in] data contiguous buffer to write (the first dimension is stride
 * dim[1])
 */
template<class T>
void
write_dataset_2D(hid_t file_handle,
  const std::string & dset_name,
  std::array<hsize_t, 2> dim,
  const T * data) {
  hid_t dset = create_dataset<T, 2>(file_handle, dim, dset_name);
  hid_t dtype = type_equiv<T>::h5_type();
  T * buf = new T[dim[0] * dim[1]];
  herr_t status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  // clog_assert(status >= 0, "Failed to write dataset \"" << dset_name << "\"");
  status = H5Dclose(dset);
  // clog_assert(status >= 0, "Failed to close dataset \"" << dset_name << "\"");
  delete[] buf;
}

/**
 * \brief attach dimension scale to a dataset
 *
 * This attaches an existing dimension scale identified by dscale_name to
 * a specified dimension for a given dataset.
 * \param[in] file_handle file to use
 * \param[in] dset_name name of dataset to attach dimension scale to
 * \param[in] dscale_name name of dimension scale to attach to dataset
 * \param[in] dim dimension index indicating which dimension to attach the scale
 */
inline void
attach_scale(hid_t file_handle,
  const std::string & dset_name,
  const std::string & dscale_name,
  unsigned int dim) {
  hid_t dset = open_dataset(file_handle, dset_name);
  hid_t ds_dset = open_dataset(file_handle, dscale_name);
  H5DSattach_scale(dset, ds_dset, dim);

  close_dataset(dset);
  close_dataset(ds_dset);
}
} // namespace h5

namespace detail {

// reading entity coordinates
template<typename T>
size_t
read_coordinates(hid_t file,
  const std::string & x_dataset_name,
  const std::string & y_dataset_name,
  std::vector<T> & entity) {

  // read the coordinate datasets individually, adding the data to the entity
  // vector
  size_t xsize = h5::read_dataset_1D(file, x_dataset_name, entity);
  size_t ysize = h5::read_dataset_1D(file, y_dataset_name, entity);

  flog_assert(xsize == ysize, "Expected x and y coordinate datasets "
                              "to have the same size");
  return xsize;
} // read_coordinates

inline void
dump_connectivity(detail::crs & connectivity) {
  #if 0
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) {
    std::cout << "dump connectivity";
    for(size_t i = 0; i < connectivity.size(); i++) {
      std::cout << "conn[" << i << "] = " << std::endl;
      const auto & tmp = connectivity[i];
      for(size_t j = 0; j < tmp.size(); j++)
        std::cout << tmp[j] << "   ";
      std::cout << std::endl;
    } // fo

  } // if
  #endif
} // dump_connectivity

} // namespace detail


////////////////////////////////////////////////////////////////////////////////
/// \brief This is the two-dimensional mesh reader and writer based on the MPAS
/// HDF5 file format.
////////////////////////////////////////////////////////////////////////////////
template<typename T>
class definition
{
public:
  //============================================================================
  // Typedefs
  //============================================================================

  //! the number of dimensions
  static constexpr std::size_t dimension() { return 2; }

  //! the byte type
  using byte_t = uint8_t;

  //! the size type
  using size_t = std::size_t;
  //! \brief the floating point type
  using real_t = T;
  //! \brief the type used for indexing arrays
  using index_t = std::size_t;

  //! \brief an alias for the vector class
  template<typename U>
  using vector = typename std::vector<U>;

  // sparse data type used from connectivity
  using crs_t = detail::crs;
  //============================================================================
  // Constructors
  //============================================================================

  //! \brief Constructor with filename
  //! \param [in] filename  The name of the file to load
  //
  // An definition should not be valid if there is no backing file, so we
  // require it to build the object.
  definition(const std::string & filename) {
    flog(info) << "Reading mesh from: " << filename << std::endl;
    read_entities(filename);
  }

  /// Default constructor (disabled)
  definition() = delete;

  /// Copy constructor (disabled)
  definition(const definition &) = delete;

  /// Assignment operator (disabled)
  definition & operator=(const definition &) = delete;

  /// Destructor
  ~definition() = default;

  //============================================================================
  //! \brief Implementation of mpas mesh read for burton specialization.
  //
  // \param[in] name Read burton mesh \e m from \e name.
  // \param[out] m Populate burton mesh \e m with contents
  // of \e name.
  //============================================================================
  void read_entities(const std::string & filename) {
    size_t cell_min, cell_max;

    // store input file information
    filename_ = filename;

    // this is a serial reading from not-partitioned HDF5 file

    //! \brief handle on the HDF5 file
    //
    // This is owned/managed solely by the definition object
    hid_t file_handle =
      h5::open_file(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    //--------------------------------------------------------------------------
    // read coordinates

    { // cells
      const std::string xCell_str("xCell");
      const std::string yCell_str("yCell");
      std::vector<real_t> cells_all;
      num_cells_ =
        detail::read_coordinates(file_handle, xCell_str, yCell_str, cells_);
    } // scope

    std::vector<real_t> coordinates;
    { // vertices
      const std::string x_dataset_name("xVertex");
      const std::string y_dataset_name("yVertex");
      num_vertices_ = detail::read_coordinates(
        file_handle, x_dataset_name, y_dataset_name, vertices_);
      num_vertices_ = num_vertices_ / dimension();

      // read connectivity info
      const std::string verticesOnCell_str{"verticesOnCell"};
      const std::string edgesOnCell_str{"edgesOnCell"};
      const std::string cellsOnCell_str{"cellsOnCell"};
      const std::string verticesOnEdge_str{"verticesOnEdge"};
      const std::string edgesOnEdge_str{"edgesOnEdge"};
      h5::read_connectivity_2D<size_t>(file_handle, verticesOnCell_str,
                                       local_connectivity_[2][0]);
      h5::read_connectivity_2D<size_t>(file_handle, edgesOnCell_str,
                                       local_connectivity_[2][1]);
      h5::read_connectivity_2D<size_t>(file_handle, cellsOnCell_str,
                                       local_connectivity_[2][2]);
      h5::read_connectivity_2D<size_t>(file_handle, verticesOnEdge_str,
                                       local_connectivity_[1][0]);
      h5::read_connectivity_2D<size_t>(file_handle, edgesOnEdge_str,
                                       local_connectivity_[1][1]);
    } // scope

    h5::close_file(file_handle);
  }


  /// Return the number of entities of a particular dimension
  /// \param [in] dim
  /// The entity dimension to query.
  size_t num_entities(size_t dim) const  {
    switch(dim) {
      case 0: {
        return vertices_.size() / dimension();
      }
    case 1: {
      return local_connectivity_.at(dim).at(0).size();
    }
      case 2: {
        return local_connectivity_.at(dim).at(0).size();
      }
      default:
        flog_fatal(
          "Dimension out of range: 0 < " << dim << " </ " << dimension());
        return 0;
    }
  }


  /// Return the set of vertices of a particular entity.
  /// \param [in] dimension  The entity dimension to query.
  /// \param [in] entity_id  The id of the entity in
  /// question.
  const std::vector<std::vector<size_t>> & entities(size_t from_dim,
    size_t to_dim) const  {
    assert(false);
    return empty_connectivity_;
  }

  const crs_t & entities_crs(size_t from_dim, size_t to_dim) const  {
    return local_connectivity_.at(from_dim).at(to_dim);
  }

  /// return the set of vertices of a particular entity.
  /// \param [in] dimension  the entity dimension to query.
  /// \param [in] entity_id  the id of the entity in
  /// question.
  std::vector<size_t>
  entities(size_t from_dim, size_t to_dim, size_t from_id) const  {
    //    assert(false);
    const auto & connectivity = local_connectivity_.at(from_dim).at(to_dim);
    auto start = connectivity.offsets[from_id];
    auto end = connectivity.offsets[from_id + 1];
    std::vector<size_t> result(
      &connectivity.indices[start], &connectivity.indices[end]);
    return result;
  }

  /// Return the vertex coordinates for a certain id.
  /// \param [in] vertex_id The
  /// id of the vertex to query.
  template<typename POINT_TYPE>
  auto vertex(size_t vertex_id) const {
    auto num_vertices = vertices_.size() / dimension();
    POINT_TYPE p;
    for(int i = 0; i < dimension(); ++i)
      p[i] = vertices_[i * num_vertices + vertex_id];
    return p;
  } // vertex

  /// Return the vertex coordinates for a certain id.
  /// \param [in] vertex_id  The id of the vertex to query.
  void vertex(size_t vertex_id, real_t * coords) const // 
  {
    constexpr auto num_dims = dimension();
    for(int i = 0; i < num_dims; ++i) {
      coords[i] = vertices_[vertex_id * num_dims + i];
    }
  } // vertex


private:
  //============================================================================
  // Private data
  //============================================================================

  size_t num_vertices_ = 0;
  size_t num_cells_ = 0;
  //  size_t num_edges_ = 0;

  //! \brief storage for cells coordinates
  vector<real_t> cells_;

  //! \brief storage for vertex coordinates
  vector<real_t> vertices_;

  //! regions
  std::vector<index_t> cell_block_id_;

  std::map<index_t, std::map<index_t, crs_t>> local_connectivity_;

  //! \brief need for now (but not used)
  std::vector<std::vector<size_t>> empty_connectivity_;

  std::vector<size_t> empty_vector_;

  std::string filename_;
};

} // namespace io
} // namespace mpas
