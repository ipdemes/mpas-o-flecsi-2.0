# MPAS-O-FleCSI-2.0

This is a FleCSI 2.0 version of the MPAS-Ocean code
The project is under development!

# Requirements

The primary requirement for building FleCSI is that you have
a C++17-capable compiler.

## Tools

You'll need the following tools in order to build FleCSI:

* Boost >= 1.70
* CMake >= 3.12
* GCC >= 8
* HDF5 with HL support
* Parmetis

 Building

MPAS-Ocean has additional FleCSI-related dependencies beyond those mentioned
above.

- [FleCSI](https://github.com/flecsi/flecsi)

# Installing the MPASO-FLeCSI dependencies with Spack

First install and initialize spack.
```
cd ${SOMEWHERE}
git clone https://github.com/spack/spack.git
cd spack
source share/spack/setup-env.sh

# load your compilers of choice and tell spack to find them
module load gcc/9.3.0
spack compiler find 
```

Clone MPAS-O-FLeCSI2 repo
```
git clone git@gitlab.lanl.gov:mpas/mpas-o-flecsi-2.0.git
```

Add the MPAS-O-FleCSI 2 spack repo (includes spackage for mpasflecsi-deps), and tell spack to install the dependencies.
```
spack repo add mpas-o-flecsi-2.0/spack-repo
spack install mpasoflecsi-deps%gcc@9.3.0 backend=legion ^openmpi@3.1.6
```

## Building MPAS-O-FleCSI:

Note, in addition to the dependencies installed by spack, you need to have cmake and mpich on your system.

```
cd mpas-o-flecsi-2.0

mkdir build

cd build

cmake ..

make -j 8 
```


