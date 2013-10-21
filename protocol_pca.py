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
from pyworkflow.utils import removeExt, removeBaseExt, makePath, moveFile
from constants import *
from spider import SpiderShell
from convert import locationToSpider
from glob import glob

      
class SpiderDefPCA(Form):        
    """Definition of parameters for Spider protocol: Align APSR"""
    def __init__(self):
        Form.__init__(self)

        self.addSection(label='Input')
        self.addParam('inputImages', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfParticles',
                      help='Input images to perform PCA')
        
        #self.addParam('maskType', )
              
        self.addParam('filterRadius1', FloatParam, default=0.6,
                      label='Fourier radius for input image',
                      help='Fourier radius for input image.')        
        self.addParam('sdFactor', FloatParam, default=0.1,
                      label='First threshold',
                      help='first threshold == image average plus this number * s.d.')
        self.addParam('filterRadius2', FloatParam, default=0.1,
                      label='Fourier radius for initial binary mask',
                      help='Fourier radius for initial binary mask.')       
        self.addParam('maskThreshold', FloatParam, default=0.01,
                      label='Mask threshold',
                      help='Threshold for filtered mask.')

# TODO: Remove from ProtAlign, and put in other category     
class SpiderProtPCA(ProtAlign):
    """ Reference-free alignment shift and rotational alignment of an image series. 
    Uses Spider AP SR command.
    """
    _definition = SpiderDefPCA()
    
    def __init__(self):
        ProtAlign.__init__(self)
        self._params = {'ext': 'stk',
                        'inputImage': 'input_image',
                        'outputMask': 'output_mask'
                        }
    
    def _defineSteps(self):
        # Define some names
        # Insert processing steps
        self.outFn = self._getPath('%(inputImage)s.%(ext)s' % self._params)
        self.inputImg = self.inputImage.get()
        index, filename = self.inputImg.getLocation()
        self._insertFunctionStep('convertInput', index, filename)
        self._insertFunctionStep('createMask', 
                                 self.filterRadius1.get(), self.sdFactor.get(),
                                 self.filterRadius2.get(), self.maskThreshold.get())
        #self._insertFunctionStep('createOutput')
        
    def convertInput(self, index, filename):
        """ Convert the input image to a Spider (with stk extension). """
        ImageHandler().convert((index, filename), (1, self.outFn))
        
    def createMask(self, filterRadius1, sdFactor, filterRadius2, maskThreshold):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        self._params.update(locals()) # Store input params in dict
        

        self._enterWorkingDir() # Do operations inside the run working dir

        spi = SpiderShell(ext=self._params['ext'], log='script.stk') # Create the Spider process to send commands 
        spi.runScript('PCA.txt', self._params)
        spi.close(end=False)
        
        self._leaveWorkingDir() # Go back to project dir

        maskFn = self._getPath('%(outputMask)s.%(ext)s' % self._params )
        maskSet = self._createSetOfParticles()
        maskSet.copyInfo(self.inputImg)
        
        for i in range(1, 8): # 7 images in mask stack
            img = Particle()
            img.setLocation(i, maskFn)
            maskSet.append(img)
            
        maskSet.write()
        self._defineOutputs(outputMask=maskSet)
            
    def _summary(self):
        summary = []
        return summary
    

