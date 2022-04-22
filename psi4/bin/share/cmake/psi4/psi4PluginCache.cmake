# psi4PluginCache.cmake
# ---------------------
#
# This module sets some likely variable values to initialize the CMake cache in your plugin.
# See ``psi4 --plugin-compile`` for use.
#

set(CMAKE_C_COMPILER          "/usr/local/pace-apps/spack/packages/0.13/linux-rhel7-cascadelake/intel-19.0.5/mvapich2-2.3.2-hpgbkqoytbjh35qn2t63rdorepxcezek/bin/mpicc" CACHE STRING "")
set(CMAKE_C_FLAGS             " -xHost" CACHE STRING "")
set(CMAKE_CXX_COMPILER        "/usr/local/pace-apps/spack/packages/0.13/linux-rhel7-cascadelake/intel-19.0.5/mvapich2-2.3.2-hpgbkqoytbjh35qn2t63rdorepxcezek/bin/mpicxx" CACHE STRING "")
set(CMAKE_CXX_FLAGS           " -xHost" CACHE STRING "")
set(CMAKE_Fortran_COMPILER    "" CACHE STRING "")
set(CMAKE_Fortran_FLAGS       "" CACHE STRING "")

#set(CMAKE_INSTALL_PREFIX      "/storage/home/hcoda1/9/jpederson6/p-jmcdaniel43-0/rich_project_chem-mcdaniel/experimental/psi4/build/stage" CACHE PATH "")
set(CMAKE_INSTALL_LIBDIR      "lib" CACHE STRING "")
set(CMAKE_INSTALL_BINDIR      "bin" CACHE STRING "")
set(CMAKE_INSTALL_DATADIR     "share" CACHE STRING "")
set(CMAKE_INSTALL_INCLUDEDIR  "include" CACHE STRING "")
set(PYMOD_INSTALL_LIBDIR      "/" CACHE STRING "")

set(CMAKE_INSTALL_MESSAGE     "LAZY" CACHE STRING "")
set(pybind11_DIR              "/storage/home/hcoda1/9/jpederson6/.conda/envs/ase_qm/share/cmake/pybind11" CACHE PATH "")

set(PYTHON_VERSION_MAJORMINOR "3.8" CACHE STRING "")
set(Python_VERSION_MAJORMINOR "3.8" CACHE STRING "")
set(PYTHON_EXECUTABLE         "/storage/home/hcoda1/9/jpederson6/.conda/envs/ase_qm/bin/python3.8" CACHE STRING "")
set(Python_EXECUTABLE         "/storage/home/hcoda1/9/jpederson6/.conda/envs/ase_qm/bin/python3.8" CACHE STRING "")

