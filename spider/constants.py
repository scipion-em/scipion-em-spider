# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

# ----------------- Constants values --------------------------------------
SPIDER_HOME = 'SPIDER_HOME'
SPPROC_DIR = 'SPPROC_DIR'
SPMAN_DIR = 'SPMAN_DIR'
SPBIN_DIR = 'SPBIN_DIR'
SPIDER = 'SPIDER'
SPIDER_MPI = 'SPIDER_MPI'

# spider documentation url
SPIDER_DOCS = 'https://spider.wadsworth.org/spider_doc/spider/docs/man/'

# convert.py constants
SHIFTX = 'shiftx'
SHIFTY = 'shifty'

ANGLE_PSI = 'psi'  # in-plane, xmipp psi
ANGLE_THE = 'the'  # tilt in xmipp
ANGLE_PHI = 'phi'  # rot in xmipp

FLIP = 'flip'

# Filter types
FILTER_TOPHAT = 0
FILTER_SPACE_REAL = 1
FILTER_FERMI = 2
FILTER_BUTTERWORTH = 3
FILTER_RAISEDCOS = 4

FILTER_LOWPASS = 0
FILTER_HIGHPASS = 1

# CA-PCA protocol
CA = 0
PCA = 1
IPCA = 2

# Center of gravity values
CG_NONE = 0
CG_PH = 1
CG_RT180 = 2

# BP type (reconstruct protocol)
BP_32F = 0
BP_3F = 1

# Protocol type (projmatch)
DEF_GROUPS = 0
GOLD_STD = 1

# Backprojection method (projmatch protocol)
BP_CG = 0
# BP_3F = 1
BP_RP = 2
BP_3N = 3

#################################################
# viewer constants
ITER_LAST = 0
ITER_SELECTION = 1

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

VOL = 0
VOL_HALF1 = 1
VOL_HALF2 = 2
VOL_FILTERED = 3
VOL_CENTERED = 4

# Template volume names depending on the iteration
VOLNAMES_GOLDSTD = {
    VOL: 'vol_%02d_unfilt',
    VOL_HALF1: 'vol_%02d_s1',
    VOL_HALF2: 'vol_%02d_s2',
    VOL_FILTERED: 'vol_%02d',
    VOL_CENTERED: 'vol_%02d_cent'
}

VOLNAMES_DEFGROUPS = {
    VOL: 'vol%02d',
    VOL_HALF1: 'vol%02d_sub1',
    VOL_HALF2: 'vol%02d_sub2',
}
