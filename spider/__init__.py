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

import os

import pwem
from pyworkflow.utils import Environ

from .constants import *

__version__ = '3.1.5'
_logo = "spider_logo.png"
_references = ['Shaikh2008', 'Frank1996b']


class Plugin(pwem.Plugin):
    _homeVar = SPIDER_HOME
    _pathVars = [SPIDER_HOME]
    _supportedVersions = ['26.06']
    _url = "https://github.com/scipion-em/scipion-em-spider"

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(SPIDER_HOME, 'spider-26.06')
        cls._defineVar(SPIDER, 'spider_linux_mp_intel64')
        cls._defineVar(SPIDER_MPI, 'spider_linux_mpi_opt64')

    @classmethod
    def getEnviron(cls):
        """ Load the environment variables needed for Spider.
        If SPIDER_HOME is defined, the bin, man and proc folders will be
        defined from it. If not, each of them should be defined separately.
        """
        env = Environ(os.environ)
        if cls.getHome():
            env.update(  # Spider needs this extra slash at the end
                {SPBIN_DIR: cls.getHome('spider/bin') + '/',
                 SPMAN_DIR: cls.getHome('spider/man') + '/',
                 SPPROC_DIR: cls.getHome('spider/proc') + '/'})
        else:
            errors = ''
            for var in [SPBIN_DIR, SPMAN_DIR, SPPROC_DIR]:
                if var not in env:
                    errors += "\n   Missing SPIDER variable: '%s'" % var
            if len(errors):
                print("ERRORS: " + errors)

        env.set('PATH', env[SPBIN_DIR], env.END)
        return env

    @classmethod
    def getScript(cls, *paths):
        """ Return the script that will be used. """
        cmd = os.path.join(__path__[0], 'scripts', cls.getActiveVersion(), *paths)
        return str(cmd)

    @classmethod
    def getProgram(cls, mpi=False):
        if mpi:
            program = os.path.basename(cls.getVar(SPIDER_MPI))
        else:
            program = os.path.basename(cls.getVar(SPIDER))

        cmd = os.path.abspath(os.path.join(cls.getEnviron()[SPBIN_DIR], program))

        return str(cmd)

    @classmethod
    def defineBinaries(cls, env):
        env.addPackage('spider', version='26.06',
                       url='https://github.com/spider-em/SPIDER/releases/'
                           'download/v26.06/spiderweb.26.06.tar.gz',
                       createBuildDir=True,
                       buildDir='spider',
                       target="spider/spider",
                       neededProgs=['csh'],
                       default=True)
