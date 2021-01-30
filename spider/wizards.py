# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se)
# *              Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
import tkinter as tk
from tkinter import ttk

from pyworkflow.utils.path import cleanPath
import pyworkflow.gui.dialog as dialog
from pyworkflow.gui.widgets import LabelSlider, HotButton
from pwem.constants import UNIT_PIXEL
from pwem.emlib.image import ImageHandler
from pwem.wizards import (EmWizard, ParticleMaskRadiusWizard,
                          ParticlesMaskRadiiWizard,
                          FilterParticlesWizard, DownsampleDialog,
                          ImagePreviewDialog, ListTreeProvider,
                          MaskRadiiPreviewDialog)

from . import Plugin
from .utils import SpiderShell, runCustomMaskScript
from .constants import FILTER_FERMI
from .convert import locationToSpider
from .protocols import (SpiderProtCAPCA, SpiderProtAlignAPSR,
                        SpiderProtAlignPairwise, SpiderProtFilter,
                        SpiderProtCustomMask, SpiderProtRefinement)


# =============================================================================
# MASKS
# =============================================================================

class SpiderProtMaskWizard(ParticleMaskRadiusWizard):
    _targets = [(SpiderProtCAPCA, ['radius']),
                (SpiderProtRefinement, ['radius'])]

    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        protParams = dict()
        if protocol.getClassName().endswith('Refinement'):
            protParams['input'] = protocol.input3DReference
        else:
            protParams['input'] = protocol.inputParticles
        protParams['label'] = label
        protParams['value'] = value
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input'] 
        return ParticleMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form, *args):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        ParticleMaskRadiusWizard.show(self, form, _value, _label)
        
    
class SpiderParticlesMaskRadiiWizard(ParticlesMaskRadiiWizard):
    _targets = [(SpiderProtAlignAPSR, ['innerRadius', 'outerRadius']),
                (SpiderProtAlignPairwise, ['innerRadius', 'outerRadius'])]        
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = dict()
        protParams['input'] = protocol.inputParticles
        protParams['label'] = label
        protParams['value'] = value
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return ParticlesMaskRadiiWizard._getListProvider(self, _objs)
    
    def show(self, form, *args):
        params = self._getParameters(form.protocol)
        value = params['value']
        label = params['label']
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            d = MaskRadiiPreviewDialog(form.root,
                                       provider,
                                       innerRadius=value[0],
                                       outerRadius=value[1],
                                       unit=UNIT_PIXEL)
            if d.resultYes():
                form.setVar(label[0], int(d.getRadius(d.radiusSliderIn)))
                form.setVar(label[1], int(d.getRadius(d.radiusSliderOut)))
        else:
            dialog.showWarning("Empty input", "Select elements first", form.root)
    
# =============================================================================
# FILTERS
# =============================================================================


class SpiderFilterParticlesWizard(FilterParticlesWizard):    
    _targets = [(SpiderProtFilter, ['filterRadius', 'lowFreq',
                                    'highFreq', 'temperature'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = dict()
        protParams['input'] = protocol.inputParticles
        protParams['label'] = label
        protParams['value'] = value
        protParams['mode'] = [protocol.filterType.get(),
                              protocol.filterMode.get(),
                              protocol.usePadding.get()]

        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return FilterParticlesWizard._getListProvider(self, _objs)
    
    def show(self, form, *args):
        protocol = form.protocol
        provider = self._getProvider(protocol)

        installErrors = Plugin.validateInstallation()

        if installErrors:
            dialog.showError("SPIDER not properly installed.",
                             "\n".join(installErrors),  form.root)
            return

        if provider is not None:
            d = SpiderFilterDialog(form.root, provider, 
                                   protocolParent=protocol)
            if d.resultYes():
                if protocol.filterType <= FILTER_FERMI:
                    form.setVar('filterRadius', d.getRadius())
                else:
                    form.setVar('lowFreq', d.getLowFreq())
                    form.setVar('highFreq', d.getHighFreq())
                if protocol.filterType == FILTER_FERMI:
                    form.setVar('temperature', d.getTemperature())
        else:
            dialog.showWarning("Input particles",
                               "Select particles first", form.root)

# =============================================================================
# UTILS
# =============================================================================


class SpiderFilterDialog(DownsampleDialog):
    
    def _beforePreview(self):
        ImagePreviewDialog._beforePreview(self)
        self.lastObj = None
        self.rightPreviewLabel = "Filtered particle"
        self.message = "Filtering particle..."
        self.previewLabel = "Particle"
        self.rightImage = ImageHandler()._img
        
    def _createControls(self, frame):
        self.freqFrame = ttk.LabelFrame(frame, text="Frequencies", padding="5 5 5 5")
        self.freqFrame.grid(row=0, column=0)
        if self.protocolParent.filterType <= FILTER_FERMI:
            self.radiusSlider = self.addFreqSlider('Radius',
                                                   self.protocolParent.filterRadius.get(),
                                                   col=0)
        else:
            self.lfSlider = self.addFreqSlider('Low freq',
                                               self.protocolParent.lowFreq.get(),
                                               col=0)
            self.hfSlider = self.addFreqSlider('High freq',
                                               self.protocolParent.highFreq.get(),
                                               col=1)
        if self.protocolParent.filterType == FILTER_FERMI:
            self.tempSlider = self.addFreqSlider('Temperature',
                                                 self.protocolParent.temperature.get(),
                                                 col=2)
        radiusButton = tk.Button(self.freqFrame, text='Preview',
                                 command=self._doPreview)
        radiusButton.grid(row=0, column=3, padx=5, pady=5)
        
    def _doPreview(self, e=None):
        if self.lastObj is None:
            dialog.showError("Empty selection",
                             "Select an item first before preview", self)
        else:
            self._computeRightPreview()
            
    def getRadius(self):
        return self.radiusSlider.get()
    
    def addFreqSlider(self, label, value, col):
        slider = LabelSlider(self.freqFrame, label, to=0.5,
                             value=value)
        slider.grid(row=0, column=col, padx=5, pady=5)
        return slider
    
    def getLowFreq(self):
        return self.lfSlider.get()
        
    def getHighFreq(self):
        return self.hfSlider.get()

    def getTemperature(self):
        return self.tempSlider.get()
    
    def updateFilteredImage(self):
        self.rightPreview.updateData(self.rightImage.getData())
        
    def _computeRightPreview(self):
        """ This function should compute the right preview
        using the self.lastObj that was selected
        """
        # Copy image to filter to Tmp project folder
        outputName = os.path.join("Tmp", "filtered_particle")
        outputPath = outputName + ".spi"
        cleanPath(outputPath)

        outputLoc = (1, outputPath)
        ih = ImageHandler()
        ih.convert(self.lastObj.getLocation(), outputLoc) 
                
        outputLocSpiStr = locationToSpider(1, outputName)
        
        pars = dict()
        pars["filterType"] = self.protocolParent.filterType.get()
        pars["filterMode"] = self.protocolParent.filterMode.get()
        pars["usePadding"] = self.protocolParent.usePadding.get()
        pars["op"] = "FQ"
        
        if self.protocolParent.filterType <= FILTER_FERMI:
            pars['filterRadius'] = self.getRadius()
        else:
            pars['lowFreq'] = self.getLowFreq()
            pars['highFreq'] = self.getHighFreq()
            
        if self.protocolParent.filterType == FILTER_FERMI:
            pars['temperature'] = self.getTemperature()

        filter_spider(outputLocSpiStr, outputLocSpiStr, **pars)
        
        # Get output image and update filtered image
        img = ImageHandler()._img
        locXmippStr = ImageHandler.locationToXmipp((1, outputPath))
        img.read(locXmippStr)
        self.rightImage = img
        self.updateFilteredImage()


# TODO: Refactor this function to be used also by method filterParticles
def filter_spider(inputLocStr, outputLocStr, **pars):
    """ Function to filter an image located on inputLocStr and
    write it to outputLocStr. """
     
    spi = SpiderShell()  # Create the Spider process to send commands
    filterNumber = pars["filterType"] * 2 + 1
    
    # Consider low-pass or high-pass
    filterNumber += pars["filterMode"]
    
    OP = pars["op"]
    if not pars["usePadding"]:
        OP += ' NP'
        
    args = []
    
    if pars["filterType"] <= FILTER_FERMI:
        args.append(pars['filterRadius'])
    else:
        args.append('%f %f' % (pars['lowFreq'], pars['highFreq']))
        
    if pars["filterType"] == FILTER_FERMI:
        args.append(pars['temperature'])
        
    spi.runFunction(OP, inputLocStr, outputLocStr, filterNumber, *args)
    spi.close()
    
    
# -------------- Custom mask Wizard -------------------------------------------

CUSTOMMASK_VARS = {'filterRadius1': 'First radius', 'sdFactor': 'First threshold',
                   'filterRadius2': 'Second radius', 'maskThreshold': 'Max threshold'}

# and the label to the given image
MASKRESULT_LABELS = ['input image',
                     'filtered',
                     'thresholded',
                     'filtered mask',
                     'final mask',
                     'mask * image',
                     'inverted mask',
                     'inverted \nmask * image',
                     ]


class SpiderCustomMaskWizard(EmWizard):    
    _targets = [(SpiderProtCustomMask, CUSTOMMASK_VARS.keys())]
    
    def _getParameters(self, protocol):
        protParams = dict()
        protParams['input'] = protocol.inputImage
        protParams['label'] = CUSTOMMASK_VARS.keys()
        protParams['labelText'] = MASKRESULT_LABELS
        protParams['value'] = [protocol.getAttributeValue(a) for a in protParams['label']]
        
        return protParams
    
    def _getProvider(self, protocol):
        return ListTreeProvider([self._getParameters(protocol)['input'].get()])
    
    def show(self, form, *args):
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if protocol.inputImage.get():
            d = CustomMaskDialog(form.root, provider,
                                 protocolParent=protocol)
            if d.resultYes():
                for varName in CUSTOMMASK_VARS:
                    form.setVar(varName, d.getVarValue(varName))
        else:
            dialog.showWarning("Input error",
                               "Select the input image first", form.root)


class CustomMaskDialog(ImagePreviewDialog):
        
    def _beforePreview(self):
        imgLocation = self.protocolParent.inputImage.get().getLocation()
        self.dim = ImageHandler().getDimensions(imgLocation)[0]
        self.lastObj = None
        self.rightPreviewLabel = "Final mask"
        self.message = "Generating mask..."
        self.ih = ImageHandler()
        self.rightImage = self.ih.createImage()
        
    def _createPreview(self, frame):
        """ Should be implemented by subclasses to 
        create the items preview. 
        """
        self._previews = []
        for i, label in enumerate(MASKRESULT_LABELS):
            self.previewLabel = label
            previewFrame = tk.Frame(frame)
            ImagePreviewDialog._createPreview(self, previewFrame)
            self._previews.append(self.preview)  # store all previews created
            previewFrame.grid(row=i//4, column=i % 4)
            
    def _itemSelected(self, obj):
        self.lastObj = obj
        dialog.FlashMessage(self, self.message,
                            func=self._computeRightPreview)
           
    def _createVarWidgets(self, parent, varName, varLabel, row, col):
        var = tk.StringVar()
        self._vars[varName] = var
        var.set(self.protocolParent.getAttributeValue(varName))
        varLabel = tk.Label(parent, text=varLabel)
        varLabel.grid(row=row, column=col*2, padx=5, pady=5)
        varEntry = tk.Entry(parent, width=10, textvariable=var)
        varEntry.grid(row=row, column=col*2+1, padx=5, pady=5)
    
    def _createControls(self, frame):
        self._vars = {}
        inputFrame = tk.Frame(frame)
        inputFrame.grid(row=0, column=0)
        
        for i, varName in enumerate(CUSTOMMASK_VARS):
            self._createVarWidgets(inputFrame, varName,
                                   CUSTOMMASK_VARS[varName], i % 2, i//2)
            
        previewBtn = HotButton(frame, text='Preview',
                               command=self._computeRightPreview)
        previewBtn.grid(row=1, column=1, padx=5, pady=5)
            
    def getVarValue(self, varName):
        return self._vars[varName].get()
    
    def _computeRightPreview(self, e=None):
        """ This function should compute the right preview
        using the self.lastObj that was selected
        """
        prot = self.protocolParent  # short notation
        tmp = prot.getProject().getTmpPath()
        ext = prot.getExt()
        # Convert input image to spider
        imgPrefix = 'inputImage'
        imgName = '%s.%s' % (imgPrefix, ext)
        imgFn = os.path.join(tmp, imgName)
        self.ih.convert(self.lastObj, (1, imgFn))
        
        runCustomMaskScript(self.getVarValue('filterRadius1'), 
                            self.getVarValue('sdFactor'), 
                            self.getVarValue('filterRadius2'), 
                            self.getVarValue('maskThreshold'), 
                            workingDir=tmp, ext=ext,
                            inputImage=imgPrefix+'@1')
        
        for i, preview in enumerate(self._previews):
            if i == 0:
                self.rightImage.read(imgFn)
            else:
                self.rightImage.read('%d@%s/stkmask.%s' % (i, tmp, ext))
            preview.updateData(self.rightImage.getData())
