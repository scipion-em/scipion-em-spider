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
from os.path import join, dirname, abspath
import datetime
from collections import OrderedDict
import subprocess
import re

from pyworkflow.utils import runJob
from pyworkflow.utils.path import replaceBaseExt, removeBaseExt

from . import Plugin


END_HEADER = 'END BATCH HEADER'
PATH = abspath(dirname(__file__))
TEMPLATE_DIR = 'templates'

# Regular expressions for parsing vars in scripts header

# Match two string types:
# [key] = value ; some comment
# GLO [key] = 'value' ; some comment
REGEX_KEYVALUE = re.compile("(?P<prefix>\s+GLO\s|[^[a-zA-Z0-9_-]*)(?P<var>\[?[a-zA-Z0-9_-]+\]?)(?P<s1>\s*)=(?P<s2>\s*)(?P<value>\S+)(?P<suffix>\s+.*)")

# Match strings of the type [key]value
# just before a 'fr l' line
REGEX_KEYFRL = re.compile("(?P<var>\[?[a-zA-Z0-9_-]+\]?)(?P<value>\S+)(?P<suffix>\s+.*)")

HEADER_COLUMNS = ['ANGLE_PSI2', 'ANGLE_THE',
                  'ANGLE_PHI', 'REF', 'EXP', 'ANGLE_PSI', 'SHIFTX',
                  'SHIFTY', 'NPROJ', 'DIFF', 'CCROT', 'ROT',
                  'SX', 'SY', 'MIR-CC']


def _getFile(*paths):
    return join(PATH, *paths)


def __substituteVar(match, paramsDict, lineTemplate):
    if match and match.groupdict()['var'] in paramsDict:
        d = match.groupdict()
        d['value'] = paramsDict[d['var']]
        return lineTemplate % d
    return None
    
    
def writeScript(inputScript, outputScript, paramsDict):
    """ Create a new Spider script by substituting 
    params in the input 'paramsDict'.
    """
    fIn = open(Plugin.getScript(inputScript), 'r', encoding='utf-8')
    fOut = open(outputScript, 'w', encoding='utf-8')
    inHeader = True  # After the end of header, no more value replacement
    inFrL = False

    for i, line in enumerate(fIn):
        if END_HEADER in line:
            inHeader = False
        if inHeader:
            try:
                newLine = __substituteVar(REGEX_KEYVALUE.match(line), paramsDict,
                                          "%(prefix)s%(var)s%(s1)s=%(s2)s%(value)s%(suffix)s\n")
                if newLine is None and inFrL:
                    newLine = __substituteVar(REGEX_KEYFRL.match(line), paramsDict,
                                              "%(var)s%(value)s%(suffix)s\n")
                if newLine:
                    line = newLine
            except Exception as ex:
                print(ex, "in line (%d): %s" % (i+1, line))
            inFrL = line.lower().startswith("fr ")
        fOut.write(line)
    fIn.close()
    fOut.close()    
     
    
def runTemplate(inputScript, ext, paramsDict, nummpis=1,
                program=None, log=None, cwd=None):
    """ This function will create a valid Spider script
    by copying the template and replacing the values in dictionary.
    After the new file is read, the Spider interpreter is invoked.
    Usually the execution should be done where the results will
    be left.
    """
    if program is None:
        program = Plugin.getProgram()

    outputScript = replaceBaseExt(inputScript, ext)
    
    if cwd is not None:
        outputScript = join(cwd, outputScript)
        
    # First write the script from the template with the substitutions
    writeScript(inputScript, outputScript, paramsDict)
    # Then proceed to run the script
    runScript(outputScript, ext, program, nummpis, log, cwd)
    

def runScript(inputScript, ext, program, nummpis, log=None, cwd=None):
    scriptName = removeBaseExt(inputScript)
    args = " %s @%s" % (ext, scriptName)
    runJob(log, program, args, numberOfMpi=nummpis,
           env=Plugin.getEnviron(), cwd=cwd)
    

def runCustomMaskScript(filterRadius1, sdFactor,
                        filterRadius2, maskThreshold,
                        workingDir, ext='stk',
                        inputImage='input_image',
                        outputMask='stkmask'):
    """ Utility function to run the custommask.msa script.
    This function will be called from the custom mask protocol
    and from the wizards to create the mask.
    """
    params = {'[filter-radius1]': filterRadius1,
              '[sd-factor]': sdFactor,
              '[filter-radius2]': filterRadius2,
              '[mask-threshold2]': maskThreshold,
              '[input_image]': inputImage,
              '[output_mask]': outputMask,
              } 
    # Run the script with the given parameters
    runTemplate('mda/custommask.msa', ext, params, cwd=workingDir)
    
    
class SpiderShell(object):
    """ This class will open a child process running Spider interpreter
    and will keep connection to send commands. 
    """
    def __init__(self, ext='spi', **kwargs):
        self._debug = kwargs.get('debug', True)
        self._log = kwargs.get('log', None)
        cwd = kwargs.get('cwd', None)
        
        FNULL = open(os.devnull, 'w')
        cmd = Plugin.getProgram().split()
        self._proc = subprocess.Popen(cmd,
                                      stdin=subprocess.PIPE,
                                      stdout=FNULL, stderr=FNULL,
                                      env=Plugin.getEnviron(),
                                      cwd=cwd,
                                      universal_newlines=True)
        if self._debug and self._log:
            self._log = open(self._log, 'w+')
            
        self.runCmd(ext)
        
    def runFunction(self, funcName, *args):
        cmd = funcName
        for a in args:
            cmd += '\n' + str(a)
        self.runCmd(cmd)
    
    def runCmd(self, cmd):
        if self._debug:
            print(cmd, file=self._log)
        self._proc.stdin.write(str(cmd) + '\n')
        self._proc.stdin.flush()
        
    def close(self, end=True):
        if end:
            self.runCmd("end")
        self._proc.wait()
        # self._proc.kill() # TODO: Check if necessary
        

class SpiderDocFile(object):
    """ Handler class to read/write spider docfile. """
    def __init__(self, filename, mode='r'):
        self._file = open(filename, mode)
        self._count = 0

    def nowisthetime(self, dt=None, fmt='%d-%b-%Y AT %H:%M:%S'):
        if dt is None:
            dt = datetime.datetime.now()
        return dt.strftime(fmt).upper()

    def fixHeaders(self, headers):
        """make all headers 11 characters in width; return doc string"""
        w = 11
        docstr = " ; /    "
        for h in headers:
            d = len(h)
            if d > w:
                h = h[:w]
            docstr += h.rjust(w + 1)
        docstr += "\n"
        return docstr

    def makeDocfileHeader(self, filename, batext=None):
        """create the comment line used at the top of SPIDER document files"""
        filename = os.path.basename(filename)
        fn, ext = os.path.splitext(filename)
        ext = ext[1:]
        if batext is None:
            batext = 'spl'  # Spider Python Library
        timestr = self.nowisthetime()
        h = " ;%s/%s   %s   %s\n" % (batext, ext, timestr, filename)
        return h
        
    def writeComment(self, filename, batext='spi'):
        line = self.makeDocfileHeader(filename, batext)

        print(line, file=self._file)

    def writeHeader(self, fields):
        line = self.fixHeaders(fields)

        print(line, file=self._file)
        
    def writeValues(self, *values):
        """ Write values in spider docfile. """
        self._count += 1
        # write data lines
        line = "%5d %2d" % (self._count, len(values))
        for v in values:
            line += " %11g" % float(v)
            
        print(line, file=self._file)
        
    def iterValues(self):
        for line in self._file:
            line = line.strip()
            if not line.startswith(';'):
                values = [float(s) for s in line.split()[2:]]
                yield values
                
    def __iter__(self):
        return self.iterValues()

    def close(self):
        self._file.close()


class SpiderDocAliFile(object):
    """ Handler class to read Spider alignment metadata."""
    def __init__(self, filename, mode='r'):
        self._file = open(filename, mode)
        self._count = 0

    def __iter__(self):
        """PSI, THE, PHI, REF#, EXP#, CUM.{ROT, SX, SY}, NPROJ, DIFF, CCROT, ROT, SX, SY, MIR-CC
        """
        for line in self._file:
            line = line.strip()
            if not line.startswith(';'):
                row = OrderedDict(zip(HEADER_COLUMNS, line.split()[2:]))
                yield row

    def close(self):
        self._file.close()
        
     
def getDocsLink(op, label):
    from .constants import SPIDER_DOCS
    """ Return a label for documentation url of a given command. """
    return '[[%(SPIDER_DOCS)s/%(op)s.html][%(label)s]]' % locals()
