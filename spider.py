# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package will contains Spider protocols
"""
import os
from os.path import join, dirname, abspath, exists, basename
from pyworkflow.object import String
from pyworkflow.em.data import EMObject
from pyworkflow.em import EMProtocol
from pyworkflow.utils.path import copyFile, removeExt, replaceExt
from pyworkflow.utils import runJob
import subprocess

END_HEADER = 'END BATCH HEADER'

SPIDER = 'spider_linux_mp_intel64'

def loadEnvironment():
    """ Load the environment variables needed for Spider.
    If SPIDER_DIR is defined, the bin, man and proc folders will be 
    defined from it. If not, each of them should be defined separately. 
    """
    global SPIDER
    SPIDER_DIR = os.environ.get('SPIDER_DIR', None) # Scipion definition
    
    if SPIDER_DIR is None:
        errors = ''
        for var in ['SPBIN_DIR', 'SPMAN_DIR', 'SPPROC_DIR']:
            if not var in os.environ:
                errors += "\n   Missing SPIDER variable: '%s'" % var
        if len(errors):
            print "ERRORS: " + errors
    else: 
        os.environ['SPBIN_DIR'] = join(SPIDER_DIR, 'bin', '')
        os.environ['SPMAN_DIR'] = join(SPIDER_DIR, 'man', '')
        os.environ['SPPROC_DIR'] = join(SPIDER_DIR, 'proc', '')
    
    # Get the executable or 'spider' by default
    SPIDER = join(os.environ['SPBIN_DIR'], os.environ.get('SPIDER', 'spider_linux_mp_intel64'))
    # expand ~ and vars
    SPIDER = abspath(os.path.expanduser(os.path.expandvars(SPIDER)))
    # Check that executable exists
    if not os.path.exists(SPIDER):
        msg = "SPIDER executable not found at:\n   '%s'" % SPIDER
        msg += "\nPlease create a link inside the bin folder: \n   '%s'" % os.environ['SPBIN_DIR']
        msg += "\n named 'spider' or define the SPIDER environment variable"
        raise Exception(msg)
        
    #
    #TODO: maybe validate that the 
    os.environ['PATH'] = os.environ['PATH'] + os.pathsep + os.environ['SPBIN_DIR']
    
    

TEMPLATE_DIR = abspath(join(dirname(__file__), 'templates'))
        
def getTemplate(templateName):
    """ Return the path to the template file given its name. """
    templateFile = join(TEMPLATE_DIR, templateName)
    if not exists(templateFile):
        raise Exception("getTemplate: template '%s' was not found in templates directory" % templateName)
    
    return templateFile

def copyTemplate(templateName, destDir):
    """ Copy a template file to a diretory """
    template = getTemplate(templateName)
    templateDest = join(destDir, basename(template))
    copyFile(template, templateDest)
    
def runSpiderTemplate(templateName, ext, paramsDict):
    """ This function will create a valid Spider script
    by copying the template and replacing the values in dictionary.
    After the new file is read, the Spider interpreter is invoked.
    """
    loadEnvironment()
    copyTemplate(templateName, '.')
    scriptName = replaceExt(templateName, ext)

    fIn = open(templateName, 'r')
    fOut = open(scriptName, 'w')
    replace = True # After the end of header, not more value replacement
    
    for i, line in enumerate(fIn):
        if END_HEADER in line:
            replace = False
        if replace:
            try:
                line = line % paramsDict
            except Exception, ex:
                print ex, "on line (%d): %s" % (i+1, line)
                raise ex
        fOut.write(line)
    fIn.close()
    fOut.close()    

    scriptName = removeExt(scriptName)  
    runJob(None, SPIDER, "%(ext)s @%(scriptName)s" % locals())


class SpiderShell(object):
    """ This class will open a child process running Spider interpreter
    and will keep conection to send commands. 
    """
    def __init__(self, ext='spi', **args):
        self._debug = args.get('debug', True)
        self._log = args.get('log', None)
        
        loadEnvironment()
        FNULL = open(os.devnull, 'w')
        self._proc = subprocess.Popen(SPIDER, shell=True, 
                                      stdin=subprocess.PIPE,
                                      stdout=FNULL, stderr=FNULL)
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
            #print "SPIDER: ", cmd
            print >> self._log, cmd
        print >> self._proc.stdin, cmd
        self._proc.stdin.flush()
        
    def runScript(self, templateName, paramsDict):
        """ Run all lines in the template script after replace the params
        with their values.
        """
        templateFile = getTemplate(templateName)        
        f = open(templateFile, 'r')
        replace = True # After the end of header, not more value replacement
        
        for line in f:
            line = line.strip()
            
            if END_HEADER in line:
                replace = False
            
            if not line.startswith(';'): # Skip comment lines
                try:
                    if replace:
                        line = line % paramsDict
                except Exception, ex:
                    print ex, "on line: ", line
            self.runCmd(line)
        
        f.close()
        if self._debug and self._log:
            self._log.close()
        
    def close(self, end=True):
        if end:
            self.runCmd("end")
        self._proc.wait()
        # self._proc.kill() TODO: Check if necesary
        

class SpiderDocFile(object):
    """ Handler class to read/write spider docfile. """
    def __init__(self, filename, mode='r'):
        self._file = open(filename, mode)
        self._count = 0
        
    def writeValues(self, *values):
        """ Write values in spider docfile. """
        self._count += 1
            # write data lines
        line = "%5d %2d" % (self._count, len(values))
        for v in values:
            line += " %11g" % float(v)
            
        print >> self._file, line
        
    def iterValues(self):
        for line in self._file:
            line = line.strip()
            if not line.startswith(';'):
                values = [float(s) for s in line.split()[2:]]
                yield values

    def close(self):
        self._file.close()
        
        
class PcaFile(EMObject):
    """ This is a container of files produced by CA PCA Spider protocol.
    It is possible to use the cas_IMC or cas_SEQ files.
    """
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        
        self.filename = String()
        
     
class SpiderProtocol(EMProtocol):
    """ Sub-class of EMProtocol to group some common Spider utils. """
            
    def convertInput(self, attrName, stackFn, selFn):
        """ Convert from an input pointer of SetOfImages to Spider.
        Params:
            attrName: the attribute name of the input pointer
            stackFn: the name of the stack for converted images
            selFn: the name of the selection file.
        """
        imgSetPointer = getattr(self, attrName)
        from convert import writeSetOfImages
        writeSetOfImages(imgSetPointer.get(), stackFn, selFn)
        
        
    def _getFileName(self, key):
        """ Give a key, append the extension
        and prefix the protocol working dir. 
        """
        template = '%(' + key + ')s.%(ext)s'
        
        return self._getPath(template % self._params)
    

    
    