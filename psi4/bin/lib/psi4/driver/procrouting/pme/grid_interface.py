#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2022 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
import copy
import itertools
import math
import sys

import numpy as np
import scipy.interpolate

from psi4 import core
from psi4.driver import constants


# ** here we access block data to set vext on grid points
#  see method _core_vbase_get_np_xyzw in ~/driver/p4util/python_helpers.py
#  to see access of Vpot.get_np_xyzw(), and we set vext points the same way
def set_blocks_vext(Vpot, extd_pot):
    """Function selecting the algorithm for a MP2 energy call
    and directing to specified or best-performance default modules.

    """
    index_start=0
    for b in range(Vpot.nblocks()):
        block = Vpot.get_block(b)
        npoints_ = block.npoints()
        extd_block = np.array(extd_pot[index_start:index_start+npoints_])
        extd_vector = core.Vector.from_array(extd_block)
        block.set_extd_pot(extd_vector)
        index_start = index_start + npoints_
    print( 'Extended potential set.\n')
    return 1


#**********************************
# this subroutine is in charge of getting/passing gridpoints/vext for
# the one electron potential in the DFT machinery.  Input **kwargs determine choice of method:
# Method 2:  Project quadrature grid to PME grid, and interpolate on PME grid
#
#  Units:  Work in Bohr within this module to avoid unit errors.
#  scf_wfn.molecule().xyz(i) returns atomic coordinates in bohr
#  scf_wfn.V_potential().get_np_xyzw()  returns grid coordinates in bohr
#**********************************
def pass_quadrature_extd_pot(
        wfn,
        qmmm_pme_gridnumber,
        extd_pot,
        interp_method,
        box,
    ):
    x, y, z, w = wfn.V_potential().get_np_xyzw()
    quadrature_grid=[]
    for i in range(len(x)):
        quadrature_grid.append( [ x[i] , y[i] , z[i] ] )
    quadrature_grid = np.array( quadrature_grid )
    print('Interpolating extended potential into quadature grid.\n')
    inverse_box = np.linalg.inv(box)
    quadrature_grid_project = project_to_PME_grid(
        quadrature_grid,
        inverse_box,
        qmmm_pme_gridnumber,
    )
    extd_pot = interpolate_PME_grid(
        qmmm_pme_gridnumber,
        extd_pot,
        quadrature_grid_project,
        interp_method,
    )
    extd_pot = list(extd_pot)
    flag = set_blocks_vext(wfn.V_potential(), extd_pot)


#**********************************
# this subroutine is in charge of getting/passing gridpoints/vext for
# the nuclear repulsion energy   Input **kwargs determine choice of method:
#
# Method 1:  Interpolate in real space using PME_grid_positions (limitation, no PBC, only cubic box ...)
def pass_nuclei_extd_pot(
        wfn,
        qmmm_pme_gridnumber,
        extd_pot,
        interp_method,
        box,
    ):
    x=[]; y=[]; z=[]
    for i in range( wfn.molecule().natom() ):
        xyz = wfn.molecule().xyz(i)
        x.append( xyz[0] )
        y.append( xyz[1] )
        z.append( xyz[2] )
    nuclei_grid=[]
    for i in range(len(x)):
        nuclei_grid.append([x[i], y[i], z[i]])
    nuclei_grid = np.array(nuclei_grid)
    print('Interpolating extended potential into nuclei grid.\n')
    inverse_box = np.linalg.inv(box)
    nuclei_grid_project = project_to_PME_grid(
        nuclei_grid,
        inverse_box,
        qmmm_pme_gridnumber,
    )
    extd_pot = interpolate_PME_grid(
        qmmm_pme_gridnumber,
        extd_pot,
        nuclei_grid_project,
        interp_method,
    )
    extd_vector = core.Vector.from_array(extd_pot)
    flag = wfn.molecule().set_extd_pot(extd_vector)


def pass_nuclei_extd_grad(wfn, qmmm_pme_gridnumber, extd_pot, box):
    nuclei_grid=[]
    for i in range(wfn.molecule().natom()):
        xyz = wfn.molecule().xyz(i)
        nuclei_grid.append([xyz[0], xyz[1], xyz[2]])
    nuclei_grid = np.array( nuclei_grid )
    print('Calculating interpolation derivatives for nuclei grid.\n')
    inverse_box = np.linalg.inv(box)
    nuclei_grid_project = project_to_PME_grid(
        nuclei_grid,
        inverse_box,
        qmmm_pme_gridnumber,
    )
    extd_pot_pad = pad_extd_pot_grid(extd_pot, qmmm_pme_gridnumber)
    xdim = np.array([i for i in range(-1,qmmm_pme_gridnumber+1)])
    grid = (xdim, xdim, xdim )
    extd_pot_3d = np.reshape(
        extd_pot_pad,
        (
            qmmm_pme_gridnumber + 2,
            qmmm_pme_gridnumber + 2,
            qmmm_pme_gridnumber + 2,
        )
    )
    # This code is largely based on
    # scipy.interpolate.RegularGridInterpolator._evaluate linear.
    interp_function = scipy.interpolate.RegularGridInterpolator(
        grid,
        extd_pot_3d,
    )
    indices, norm_dist, out_of_bounds = interp_function._find_indices(
        nuclei_grid_project.T,
    )
    edges = itertools.product(*[[i, i + 1] for i in indices])
    vslice = (slice(None),) + (None,)*(interp_function.values.ndim - len(indices))
    extd_du=[0 for i in range(3)]
    for edge_indices in edges:
        weight = 1.
        weight_du = [1. for i in range(3)]
        j=0
        for ei, i, yi in zip(edge_indices, indices, norm_dist):
            weight *= np.where(ei == i, 1 - yi, yi)
            for k in range(3):
                if j == k:
                    weight_du[k] *= np.where(ei == i, -1.0, 1.0)
                else:
                    weight_du[k] *= np.where(ei == i, 1 - yi, yi)
            j+=1
        extd_du[0] += np.asarray(interp_function.values[edge_indices]) * weight_du[0][vslice]
        extd_du[1] += np.asarray(interp_function.values[edge_indices]) * weight_du[1][vslice]
        extd_du[2] += np.asarray(interp_function.values[edge_indices]) * weight_du[2][vslice]
    extd_dr = qmmm_pme_gridnumber * np.matmul(inverse_box, extd_du)
    print('Completed interpolation derivatives for nuclei grid.\n')
    extd_grad_x = core.Vector.from_array(extd_dr[:,0])
    extd_grad_y = core.Vector.from_array(extd_dr[:,1])
    extd_grad_z = core.Vector.from_array(extd_dr[:,2])
    flag = wfn.molecule().set_extd_grad_x(extd_grad_x)
    flag = wfn.molecule().set_extd_grad_y(extd_grad_y)
    flag = wfn.molecule().set_extd_grad_z(extd_grad_z)
    return 1


#*******************************************
# this is a wrapper that calls scipy routines to
# interpolate from vext values evaluated at PME grid points.
# input "interpolation_method" should be set to either :
#       "interpn"  :: calls  scipy.interpolate.interpn()
#       "griddata" :: calls  scipy.interpolate.griddata()
def interpolate_PME_grid(
        qmmm_pme_gridnumber,
        extd_pot,
        interpolate_coords,
        interp_method,
        pad_grid=True,
        **kwargs,
    ):
    if pad_grid:
        extd_pot_pad = pad_extd_pot_grid(extd_pot, qmmm_pme_gridnumber)
        x_dim = np.array([i for i in range(-1,qmmm_pme_gridnumber+1)])
        grid = (x_dim, x_dim ,x_dim)
        extd_pot_3d = np.reshape(
            extd_pot_pad,
            (
                qmmm_pme_gridnumber + 2,
                qmmm_pme_gridnumber + 2,
                qmmm_pme_gridnumber + 2,
            ),
        )
    else:
        x_dim = np.array([i for i in range(qmmm_pme_gridnumber)])
        grid = (x_dim, x_dim, x_dim)
        extd_pot_3d = np.reshape(
            extd_pot,
            (
                qmmm_pme_gridnumber,
                qmmm_pme_gridnumber,
                qmmm_pme_gridnumber,
            ),
        )
    output_interp = scipy.interpolate.interpn(
        grid,
        extd_pot_3d,
        interpolate_coords,
        method='linear',
    )
    return output_interp


#***********************************
# this method pads the boundaries of vext_tot for interpolation.
# vext_tot satisfies the PBC of the system, but this is done so that we can avoid PBC in the interpolation.
# We might have an interpolation point at the edge of the grid, if we're
# using linear interpolation then we just need one more point to pad the edge
#    so if vext_tot     is dimension pme_grid_size x pme_grid_size x pme_grid_size,
#    then  vext_tot_pad is dimension (pme_grid_size + 2) x (pme_grid_size + 2) x (pme_grid_size + 2)
#    where the padding on either side of the grid is done by using the PBC
def pad_extd_pot_grid(extd_pot, qmmm_pme_gridnumber):
    extd_pot_pad = []
    for i in range(-1, qmmm_pme_gridnumber + 1):
        i_pbc = (i + qmmm_pme_gridnumber) % qmmm_pme_gridnumber
        for j in range(-1, qmmm_pme_gridnumber + 1):
            j_pbc = (j + qmmm_pme_gridnumber) % qmmm_pme_gridnumber
            for k in range(-1, qmmm_pme_gridnumber + 1):
                k_pbc = (k + qmmm_pme_gridnumber) % qmmm_pme_gridnumber
                index = (i_pbc * qmmm_pme_gridnumber**2
                         + j_pbc * qmmm_pme_gridnumber
                         + k_pbc)
                extd_pot_pad.append(extd_pot[index])
    extd_pot_pad = np.array(extd_pot_pad)
    return extd_pot_pad


#*********************************
# project real space points to PME grid
# this algorithm is identical to that used in method
# 'pme_update_grid_index_and_fraction' in OpenMM source code,
# ReferencePME.cpp.  See comments in ReferencePME.cpp
# about how algorithm works ...
#*********************************
def project_to_PME_grid(real_grid_points, inverse_box, qmmm_pme_gridnumber):
    scaled_grid_points = np.matmul(real_grid_points, inverse_box)
    scaled_grid_points = ((scaled_grid_points - np.floor(scaled_grid_points))
                          * qmmm_pme_gridnumber)
    scaled_grid_points = (np.mod(scaled_grid_points.astype(int), qmmm_pme_gridnumber) 
                          + (scaled_grid_points - scaled_grid_points.astype(int)))
    return scaled_grid_points
