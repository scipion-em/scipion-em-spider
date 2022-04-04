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

from pwem.protocols import EMProtocol

from .. import Plugin
from ..utils import runTemplate
from ..convert import writeSetOfImages


class SpiderProtocol(EMProtocol):
    """ Base protocol for SPIDER utils. """
    _label = None
    _params = None
            
    def convertInput(self, attrName, stackFn, selFn):
        """ Convert from an input pointer of SetOfImages to Spider.
        Params:
            attrName: the attribute name of the input pointer
            stackFn: the name of the stack for converted images
            selFn: the name of the selection file.
        """
        imgSetPointer = getattr(self, attrName)
        writeSetOfImages(imgSetPointer.get(), stackFn, selFn)
        
    def _getFileName(self, key, *args):
        """ Give a key, append the extension
        and prefix the protocol working dir. 
        """
        template = '%(' + key + ')s.%(ext)s'
        
        return self._getPath(template % self._params)
    
    def getExt(self):
        """ Return the extension used in the script,
        stored in a dictionary called self._params. 
        """
        return self._params['ext']
    
    def getScript(self):
        return getattr(self, '_script', None)
    
    def runTemplate(self, inputScript, ext, paramsDict, nummpis=1):
        """ This function will create a valid Spider script
        by copying the template and replacing the values in dictionary.
        After the new file is read, the Spider interpreter is invoked.
        """
        self._enterWorkingDir()

        log = getattr(self, '_log', None)
        mpiFlag = True if nummpis > 1 else False
        program = Plugin.getProgram(mpiFlag)
        runTemplate(inputScript, ext, paramsDict, nummpis=nummpis,
                    program=program, log=log)
        self._leaveWorkingDir()
    
        f = open(self.getLogPaths()[0], 'r')
        for line in f.readlines():
            if 'FATAL ERROR ENCOUNTERED IN BATCH MODE' in line:
                raise RuntimeError('Spider script error!')
        f.close()
