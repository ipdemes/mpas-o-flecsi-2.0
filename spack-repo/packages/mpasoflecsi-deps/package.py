# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *
import os

class MpasoflecsiDeps(BundlePackage):
    '''TODO
    '''
    homepage = 'www.canga-scidac.org'
    git      = 'https://gitlab.lanl.gov/mpas/mpas-o-flecsi-2.0'

    version('master', branch='master', submodules=False, preferred=True)

    variant('build_type', default='RelWithDebInfo',
            values=('Debug', 'Release', 'RelWithDebInfo', 'MinSizeRel'),
            description='The build type to build', multi=False)
    variant('backend', default='legion', values=('serial', 'mpi', 'legion', 'hpx'),
            description='Backend to use for distributed memory', multi=False)
#    variant('debug_backend', default=False,
#            description='Build Backend with Debug Mode')
    variant('hdf5', default=True,
            description='Enable HDF5 Support')
    variant('graphviz', default=False,
            description='Enable GraphViz Support')
    variant('flog', default=False,
            description='Enable FLOG Logging Utility')
    variant('kokkos', default=True,
            description='Enable Kokkos Support')
    variant('yamlcpp', default=True,
            description='Enable Yaml-CPP Support') 
     

    for b in ['Debug', 'Release', 'RelWithDebInfo', 'MinSizeRel']:
        depends_on("flecsi build_type=%s" % b,
            when="build_type=%s" % b)
    for b in ['serial', 'mpi', 'legion', 'hpx', 'charmpp']:
        depends_on("flecsi backend=%s" % b,
            when="backend=%s" % b)
    for v in [ 'hdf5',  'graphviz', 'kokkos', 'flog']:
        depends_on("flecsi +%s" % v, when="+%s" % v)
        depends_on("flecsi ~%s" % v, when="~%s" % v)

    depends_on('cmake@3.12:3.18.4')
    depends_on('hdf5+hl+mpi', when='+hdf5')
    depends_on('yaml-cpp', when='+yamlcpp')

    def setup_run_environment(self, env):
        if '+hdf5' in self.spec:
            env.set('HDF5_ROOT', self.spec['hdf5'].prefix)

