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
import numpy as np
from psi4 import core
from psi4.driver import constants


def set_blocks_extd_pot(Vpot, extd_pot):
    """Set extd_pot on grid points in functional blocks.

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
