# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
from os.path import abspath

import pyworkflow.em
from pyworkflow.utils import Environ, join
import objects

_logo = "spider_logo.png"
_references = ['Shaikh2008', 'Frank1996b']

SPIDER_HOME_VAR = 'SPIDER_HOME'
SPPROC_DIR = 'SPPROC_DIR'
SPMAN_DIR = 'SPMAN_DIR'
SPBIN_DIR = 'SPBIN_DIR'
SPIDER = 'spider_linux_mp_intel64'
SPIDER_MPI = 'spider_linux_mpi_opt64'


# The following class is required for Scipion to detect this Python module
# as a Scipion Plugin. It needs to specify the PluginMeta __metaclass__
# Some function related to the underlying package binaries need to be
# implemented
class Plugin:
    #__metaclass__ = pyworkflow.em.PluginMeta

    @classmethod
    def getEnviron(cls):
        """ Load the environment variables needed for Spider.
        If SPIDER_HOME is defined, the bin, man and proc folders will be
        defined from it. If not, each of them should be defined separately.
        """
        env = Environ(os.environ)
        SPIDER_HOME = os.environ[('%s' % SPIDER_HOME_VAR)]
        if SPIDER_HOME is None:
            errors = ''
            for var in [SPBIN_DIR, SPMAN_DIR, SPPROC_DIR]:
                if not var in env:
                    errors += "\n   Missing SPIDER variable: '%s'" % var
            if len(errors):
                print "ERRORS: " + errors
        else:
            env.update({SPBIN_DIR: join(SPIDER_HOME, 'bin') + '/',  # Spider needs this extra slash at the end
                        SPMAN_DIR: join(SPIDER_HOME, 'man') + '/',
                        SPPROC_DIR: join(SPIDER_HOME, 'proc') + '/'
                        })

        env.set('PATH', env[SPBIN_DIR], env.END)

        return env

    @classmethod
    def getVersion(cls):
        path = os.environ[SPIDER_HOME_VAR]
        for v in cls.getSupportedVersions():
            if v in path or v in os.path.realpath(path):
                return v
        return ''

    @classmethod
    def getSupportedVersions(cls):
        """ Return the list of supported binary versions. """
        return ['24.03']

    @classmethod
    def validateInstallation(cls):
        """ This function will be used to check if package is properly installed. """
        environ = cls.getEnviron()
        missingPaths = ["%s: %s" % (var, environ[var])
                        for var in [SPBIN_DIR, SPMAN_DIR, SPPROC_DIR]
                        if not os.path.exists(environ[var])]

        return (["Missing variables:"] + missingPaths) if missingPaths else []

    @classmethod
    def getScript(cls, *paths):
        """ Return the script that will be used. """
        cmd = join(__path__[0], 'scripts', cls.getVersion(), *paths)
        return str(cmd)

    @classmethod
    def getProgram(cls, mpi=False):
        program = SPIDER if mpi is False else SPIDER_MPI
        cmd = abspath(join(cls.getEnviron()[SPBIN_DIR], program))

        return str(cmd)

pyworkflow.em.Domain.registerPlugin(__name__)