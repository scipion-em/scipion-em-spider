# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from pyworkflow.em.protocol import ProtImportParticles
from pyworkflow.tests import setupTestProject, DataSet, unittest, BaseTest
from pyworkflow.tests.em.workflows.test_workflow import TestWorkflow

from spider.convert import writeSetOfImages
from spider.protocols import *

  

class TestSpiderConvert(TestWorkflow):
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('mda')
        cls.particlesFn = cls.dataset.getFile('particles')
    
    def test_convert(self):
        """ Run an Import particles protocol. """
        protImport = self.newProtocol(ProtImportParticles, filesPath=self.particlesFn, samplingRate=3.5)
        self.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if getattr(protImport, 'outputParticles', None) is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % self.particlesFn)
        
        stackFn = self.getOutputPath('stack.stk')
        selFn = self.getOutputPath('stack_sel.stk')
        print "stackFn: ", stackFn
        writeSetOfImages(protImport.outputParticles, stackFn, selFn)


class TestSpiderWorkflow(TestWorkflow):
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('mda')
        cls.particlesFn = cls.dataset.getFile('particles')
        
    def validateFilesExist(self, files):
        exist = []
        for f in files:
            if os.path.exists(f):
                exist.append(f)
        self.assertEqual(files, exist, "Missing output files")
    
    def test_mdaWorkflow(self):
        """ Run an Import particles protocol. """
        protImport = self.newProtocol(ProtImportParticles, filesPath=self.particlesFn, samplingRate=3.5)
        self.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % self.particlesFn)
        
        protFilter = self.newProtocol(SpiderProtFilter)
        protFilter.inputParticles.set(protImport)
        protFilter.inputParticles.setExtended('outputParticles')
        self.launchProtocol(protFilter)
        self.assertIsNotNone(protFilter.outputParticles,
                             "There was a problem with the SpiderProtFilter outputParticles")
        
        protAPSR = self.newProtocol(SpiderProtAlignAPSR)
        protAPSR.inputParticles.set(protFilter.outputParticles)
        self.launchProtocol(protAPSR)
        self.assertIsNotNone(protAPSR.outputParticles,
                             "There was a problem with the SpiderProtAlignAPSR outputParticles")
        self.assertIsNotNone(protAPSR.outputAverage,
                             "There was a problem with the SpiderProtAlignAPSR outputAverage")
        self.assertTrue(protAPSR.outputParticles.hasAlignment2D(),
                        "outputParticles have no alignment registered")
        
        protPairwise = self.newProtocol(SpiderProtAlignPairwise)
        protPairwise.inputParticles.set(protFilter.outputParticles)
        self.launchProtocol(protPairwise)       
        self.assertIsNotNone(protPairwise.outputParticles,
                             "There was a problem with the SpiderProtAlignPairwise outputParticles")
        self.assertIsNotNone(protPairwise.outputAverage,
                             "There was a problem with the SpiderProtAlignPairwise outputAverage")
        self.assertTrue(protPairwise.outputParticles.hasAlignment2D(),
                        "outputParticles have no alignment registered")
         
        protMask = self.newProtocol(SpiderProtCustomMask)
        protMask.inputImage.set(protAPSR.outputAverage)
        self.launchProtocol(protMask)       
        self.assertIsNotNone(protMask.outputMask,
                             "There was a problem with the SpiderProtCustomMask outputAverage")
              
        protCAPCA = self.newProtocol(SpiderProtCAPCA)
        protCAPCA.maskType.set(1)
        protCAPCA.maskImage.set(protMask.outputMask)
        protCAPCA.inputParticles.set(protAPSR.outputParticles)
        self.launchProtocol(protCAPCA)
        self.assertIsNotNone(protCAPCA.imcFile,
                             "There was a problem with the SpiderProtCAPCA imcFile")
        self.assertIsNotNone(protCAPCA.seqFile,
                             "There was a problem with the SpiderProtCAPCA seqFile")
        
        
        protWard = self.newProtocol(SpiderProtClassifyWard)
        protWard.pcaFile.set(protCAPCA.imcFile)
        protWard.inputParticles.set(protAPSR.outputParticles)
        self.launchProtocol(protWard)
        nativeFiles = []
        dendroFile = protWard._getFileName('dendroDoc')
        averages = protWard._getFileName('averages')
        nativeFiles.append(dendroFile)
        nativeFiles.append(averages)
        self.validateFilesExist(nativeFiles)
        
        protKmeans = self.newProtocol(SpiderProtClassifyKmeans)
        protKmeans.pcaFile.set(protCAPCA.imcFile)
        protKmeans.inputParticles.set(protAPSR.outputParticles)
        protKmeans.numberOfClasses.set(4)
        self.launchProtocol(protKmeans)
        self.assertIsNotNone(protKmeans.outputClasses,
                             "There was a problem with the SpiderProtClassifyKmeans outputClasses")
        
        protDiday = self.newProtocol(SpiderProtClassifyDiday)
        protDiday.pcaFile.set(protCAPCA.imcFile)
        protDiday.inputParticles.set(protAPSR.outputParticles)
        self.launchProtocol(protDiday)
        self.validateFilesExist(nativeFiles)
