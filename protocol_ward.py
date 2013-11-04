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
This sub-package contains Spider protocol for PCA.
"""


from pyworkflow.em import *  
from pyworkflow.utils import removeExt, removeBaseExt, makePath, moveFile, copyFile, basename
from constants import *
from spider import SpiderShell, SpiderDocFile, SpiderProtocol
from convert import locationToSpider
from glob import glob

      
# TODO: Remove from ProtAlign, and put in other category     
class SpiderProtClassifyWard(ProtClassify, SpiderProtocol):
    """ Ward's method, using 'CL HC' 
    """
    def __init__(self):
        ProtClassify.__init__(self)
        SpiderProtocol.__init__(self)
        
        self._params = {'ext': 'stk',
                        'particles': 'particles',
                        'particlesSel': 'particles_sel',
                        'dendroPs': 'dendrogram',
                        'dendroDoc': 'docdendro',
                        'averages': 'averages'
                        }
    
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True, 
                      pointerClass='SetOfParticles',
                      help='Input images to perform PCA')
        
        #form.addParam('maskType', )
              
        form.addParam('pcaFilePointer', PointerParam, pointerClass='PcaFile',
                      label="PCA file", 
                      help='IMC or SEQ file generated in CA-PCA')        
        form.addParam('numberOfFactors', IntParam, default=10,
                      label='Number of factors',
                      help='After running, examine the eigenimages and decide which ones to use.\n'
                           'Typically all but the first few are noisy.')
        
        
    def _getFileName(self, key):
        #TODO: Move to a base Spider protocol
        template = '%(' + key + ')s.%(ext)s'
        return self._getPath(template % self._params)
    
    def _defineSteps(self):
        pcaFile = self.pcaFilePointer.get().filename.get()
        
        self._insertFunctionStep('convertInput', 'inputParticles',
                                 self._getFileName('particles'), self._getFileName('particlesSel'))
        self._insertFunctionStep('classifyWard', pcaFile, self.numberOfFactors.get())
        self._insertFunctionStep('buildDendroStep')
            
    def classifyWard(self, imcFile, numberOfFactors):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        self._params.update(locals()) # Store input params in dict
        
        # Copy file to working directory, it could be also a link
        imcLocalFile = basename(imcFile)
        copyFile(imcFile, self._getPath(imcLocalFile))
        print "copy from '%s' to '%s' " % (imcFile, imcLocalFile)
        imcLocalFile = removeExt(imcLocalFile)

        self._enterWorkingDir() # Do operations inside the run working dir

        spi = SpiderShell(ext=self._params['ext'], log='script.stk') # Create the Spider process to send commands 
        spi.runFunction('CL HC', imcLocalFile, '1-%d' % numberOfFactors, 0, 5, 
                        'Y', self._params['dendroPs'], 'Y', self._params['dendroDoc'])
        spi.close()
        
        self._leaveWorkingDir() # Go back to project dir
        
        
    def buildDendroStep(self):
        self.buildDendrogram(True)
        
    def buildDendrogram(self, writeAverages=False):
        """ Parse Spider docfile with the information to build the dendogram.
        Params:
            dendroFile: docfile with a row per image. 
                 Each row contains the image id and the height.
        """ 
        dendroFile = self._getFileName('dendroDoc')
        # Dendrofile is a docfile with at least 3 data colums (class, height, id)
        doc = SpiderDocFile(dendroFile)
        values = []
        indexes = []
        for c, h, _ in doc.iterValues(): 
            indexes.append(c)
            values.append(h)
        doc.close()
        
        self.dendroValues = values
        self.dendroIndexes = indexes
        self.dendroImages = self._getFileName('particles')
        self.dendroAverages = self._getFileName('averages')
        
        return self._buildDendrogram(0, len(values)-1, 1, writeAverages)
    
    def _buildDendrogram(self, leftIndex, rightIndex, index, writeAverages=False):
        """ This function is recursively called to create the dendogram graph(binary tree)
        and also to write the average image files.
        Params:
            leftIndex, rightIndex: the indinxes within the list where to search.
            index: the index of the class average.
            writeImages: flag to select when to write averages.
        From self:
            self.dendroValues: the list with the heights of each node
            self.dendroImages: image stack filename to read particles
            self.dendroAverages: stack name where to write averages
        It will search for the max in values list (between minIndex and maxIndex).
        Nodes to the left of the max are left childs and the other right childs.
        """
        maxValue = self.dendroValues[leftIndex]
        maxIndex = 0
        for i, v in enumerate(self.dendroValues[leftIndex+1:rightIndex]):
            if v > maxValue:
                maxValue = v
                maxIndex = i+1
        
        m = maxIndex + leftIndex
        node = {'height': maxValue, 'childs': [], 
                'length': 1, 'index': index}#len(self.dendroValues[leftIndex:rightIndex])}
        
        ih = ImageHandler()

        if writeAverages:
            particleNumber = self.dendroIndexes[m+1]
            node['image'] = ih.read((particleNumber, self.dendroImages))
            node['imageList'] = [particleNumber]
            
        def addChildNode(left, right, index):
            if right > left:
                child = self._buildDendrogram(left, right, index, writeAverages)
                node['childs'].append(child)
                node['length'] += child['length'] 
                if writeAverages:
                    node['image'] += child['image']
                    node['imageList'] += child['imageList']
                    del child['image']
                
        if rightIndex > leftIndex + 1:
            addChildNode(leftIndex, m, 2*index)
            addChildNode(m+1, rightIndex, 2*index+1)
            if writeAverages:
                #TODO: node['image'] /= float(node['length'])
                ih.write(node['image'], (index, self.dendroAverages))
                fn = self._getTmpPath('doc_class%03d.stk' % index)
                f = open(fn, 'w+')
                for i in node['imageList']:
                    f.write('%s\n' % i)
                f.close()
        return node

            
    def _summary(self):
        summary = []
        return summary
    

