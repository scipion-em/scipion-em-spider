# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from pwem.protocols import ProtImportParticles
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from pwem.tests.workflows.test_workflow import TestWorkflow

from ..convert import writeSetOfImages
from ..protocols import (SpiderProtFilter, SpiderProtAlignAPSR,
                         SpiderProtAlignPairwise, SpiderProtCustomMask,
                         SpiderProtCAPCA, SpiderProtClassifyWard,
                         SpiderProtClassifyDiday, SpiderProtClassifyKmeans)


class TestSpiderConvert(TestWorkflow):
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('mda')
        cls.particlesFn = cls.dataset.getFile('particles')
    
    def test_convert(self):
        """ Run an Import particles protocol. """
        print(magentaStr("\n==> Importing data - particles:"))
        protImport = self.newProtocol(ProtImportParticles,
                                      filesPath=self.particlesFn, samplingRate=3.5)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputParticles,
                             "SetOfParticles has not been produced.")

        stackFn = self.getOutputPath('stack.stk')
        selFn = self.getOutputPath('stack_sel.stk')
        print("stackFn: ", stackFn)
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
        print(magentaStr("\n==> Importing data - particles:"))
        protImport = self.newProtocol(ProtImportParticles,
                                      filesPath=self.particlesFn, samplingRate=3.5)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputParticles,
                             "SetOfParticles has not been produced.")

        print(magentaStr("\n==> Testing spider - filter particles:"))
        protFilter = self.newProtocol(SpiderProtFilter)
        protFilter.inputParticles.set(protImport)
        protFilter.inputParticles.setExtended('outputParticles')
        self.launchProtocol(protFilter)
        self.assertIsNotNone(protFilter.outputParticles,
                             "There was a problem with the SpiderProtFilter outputParticles")

        print(magentaStr("\n==> Testing spider - align ap sr:"))
        protAPSR = self.newProtocol(SpiderProtAlignAPSR)
        protAPSR.inputParticles.set(protFilter.outputParticles)
        self.launchProtocol(protAPSR)
        self.assertIsNotNone(protAPSR.outputParticles,
                             "There was a problem with the SpiderProtAlignAPSR outputParticles")
        self.assertIsNotNone(protAPSR.outputAverage,
                             "There was a problem with the SpiderProtAlignAPSR outputAverage")
        self.assertTrue(protAPSR.outputParticles.hasAlignment2D(),
                        "outputParticles have no alignment registered")

        print(magentaStr("\n==> Testing spider - align pairwise:"))
        protPairwise = self.newProtocol(SpiderProtAlignPairwise)
        protPairwise.inputParticles.set(protFilter.outputParticles)
        self.launchProtocol(protPairwise)       
        self.assertIsNotNone(protPairwise.outputParticles,
                             "There was a problem with the SpiderProtAlignPairwise outputParticles")
        self.assertIsNotNone(protPairwise.outputAverage,
                             "There was a problem with the SpiderProtAlignPairwise outputAverage")
        self.assertTrue(protPairwise.outputParticles.hasAlignment2D(),
                        "outputParticles have no alignment registered")

        print(magentaStr("\n==> Testing spider - custom mask 2d:"))
        protMask = self.newProtocol(SpiderProtCustomMask)
        protMask.inputImage.set(protAPSR.outputAverage)
        self.launchProtocol(protMask)       
        self.assertIsNotNone(protMask.outputMask,
                             "There was a problem with the SpiderProtCustomMask outputAverage")

        print(magentaStr("\n==> Testing spider - ca pca:"))
        protCAPCA = self.newProtocol(SpiderProtCAPCA)
        protCAPCA.maskType.set(1)
        protCAPCA.maskImage.set(protMask.outputMask)
        protCAPCA.inputParticles.set(protAPSR.outputParticles)
        self.launchProtocol(protCAPCA)
        self.assertIsNotNone(protCAPCA.imcFile,
                             "There was a problem with the SpiderProtCAPCA imcFile")
        self.assertIsNotNone(protCAPCA.seqFile,
                             "There was a problem with the SpiderProtCAPCA seqFile")

        print(magentaStr("\n==> Testing spider - classify ward:"))
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

        print(magentaStr("\n==> Testing spider - classify k-means:"))
        protKmeans = self.newProtocol(SpiderProtClassifyKmeans)
        protKmeans.pcaFile.set(protCAPCA.imcFile)
        protKmeans.inputParticles.set(protAPSR.outputParticles)
        protKmeans.numberOfClasses.set(4)
        self.launchProtocol(protKmeans)
        self.assertIsNotNone(protKmeans.outputClasses,
                             "There was a problem with the SpiderProtClassifyKmeans outputClasses")

        print(magentaStr("\n==> Testing spider - classify diday:"))
        protDiday = self.newProtocol(SpiderProtClassifyDiday)
        protDiday.pcaFile.set(protCAPCA.imcFile)
        protDiday.inputParticles.set(protAPSR.outputParticles)
        self.launchProtocol(protDiday)
        self.validateFilesExist(nativeFiles)
