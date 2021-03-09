/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#pragma once

#include <mpasoflecsi/specialization/definition.hh>


namespace mpas {
namespace io {

struct mpi_ph5 {
  inline static hid_t create_file(const std::string name) {
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t ret = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
    return ret;
  }

  /**
   * \brief Write buffer (in parallel) to file.
   *
   * @param file_id the hdf5 file id
   * @param dname name of dataset to write
   * @param ldims local dimensions of buffer to write
   * @param gdims global dimensions of file to write
   * @param displ local displacement of buffer in file
   * @param buffer local buffer to write
   */
  template<unsigned short ND, class T>
  inline static void write_buffer(const hid_t file_id,
    const std::string dname,
    std::array<hsize_t, ND> ldims,
    std::array<hsize_t, ND> gdims,
    std::array<hsize_t, ND> displ,
    const T * buffer) {
    auto dtype = h5::type_equiv<T>::h5_type();
    hid_t filespace = H5Screate_simple(ND, gdims.data(), NULL);
    hid_t dset_id = H5Dcreate(file_id, dname.c_str(), dtype, filespace,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    hid_t memspace = H5Screate_simple(ND, ldims.data(), NULL);
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(
      filespace, H5S_SELECT_SET, displ.data(), NULL, ldims.data(), NULL);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    herr_t status =
      H5Dwrite(dset_id, dtype, memspace, filespace, plist_id, buffer);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
  }

  /**
   * \brief Write buffer (in parallel) to file.
   *
   * @param file_id the hdf5 file id
   * @param dname name of dataset to write
   * @param gdims global dimensions of file to write
   * @param num_elements number of local elements to write
   * @param coord global coordinates for local elements
   * @param buffer local buffer to write
   */
  template<unsigned short ND, class T>
  inline static void write_buffer(const hid_t file_id,
    const std::string dname,
    std::array<hsize_t, ND> gdims,
    std::size_t num_elements,
    const hsize_t * coord,
    const T * buffer) {
    auto dtype = h5::type_equiv<T>::h5_type();
    hid_t filespace, dset_id;

    // Open dataset (create if it does not exist)
    htri_t dset_exists = H5Lexists(file_id, dname.c_str(), H5P_DEFAULT);
    if(dset_exists > 0) {
      dset_id = H5Dopen(file_id, dname.c_str(), H5P_DEFAULT);
    }
    else {
      filespace = H5Screate_simple(ND, gdims.data(), NULL);
      dset_id = H5Dcreate(file_id, dname.c_str(), dtype, filespace, H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(filespace);
    }

    hsize_t marray[] = {num_elements};
    hid_t memspace = H5Screate_simple(1, marray, NULL);
    filespace = H5Dget_space(dset_id);
    H5Sselect_elements(filespace, H5S_SELECT_SET, num_elements, coord);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    herr_t status =
      H5Dwrite(dset_id, dtype, memspace, filespace, plist_id, buffer);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
  }

  inline static void close_file(hid_t file) {
    H5Fclose(file);
  }
};


}}
