# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se)
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

from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from pwem.protocols import ProtImportParticles
from pwem.tests.workflows.test_workflow import TestWorkflow

from ..protocols import (SpiderProtFilter, SpiderProtAlignAPSR,
                         SpiderProtAlignPairwise)

   
class TestSpiderAlign(TestWorkflow):
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('mda')
        cls.particlesFn = cls.dataset.getFile('particles')

    def runAlignment(self, protFilter, AlignmentClass, **kwargs):
        protAlign = self.newProtocol(AlignmentClass, **kwargs)
        protAlign.inputParticles.set(protFilter.outputParticles)
        self.launchProtocol(protAlign)
        className = protAlign.getClassName()
        self.assertIsNotNone(protAlign.outputParticles,
                             "There was a problem with the %s outputParticles" % className)
        self.assertIsNotNone(protAlign.outputAverage,
                             "There was a problem with the %s outputAverage" % className)
        self.assertTrue(protAlign.outputParticles.hasAlignment2D(),
                        "outputParticles have no alignment registered")
        return protAlign

    def test_align(self):
        """ Run an Import particles protocol. """
        print(magentaStr("\n==> Importing data - particles:"))
        protImport = self.newProtocol(ProtImportParticles,
                                      filesPath=self.particlesFn,
                                      samplingRate=3.5)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputParticles,
                             "SetOfParticles has not been produced.")

        print(magentaStr("\n==> Running spider - filter particles:"))
        protFilter = self.newProtocol(SpiderProtFilter)
        protFilter.inputParticles.set(protImport)
        protFilter.inputParticles.setExtended('outputParticles')
        self.launchProtocol(protFilter)
        self.assertIsNotNone(protFilter.outputParticles,
                             "There was a problem with the SpiderProtFilter outputParticles")

        print(magentaStr("\n==> Testing spider - align ap sr:"))
        self.runAlignment(protFilter, SpiderProtAlignAPSR,
                          objLabel='align apsr')
        print(magentaStr("\n==> Testing spider - align pairwise:"))
        self.runAlignment(protFilter, SpiderProtAlignPairwise,
                          objLabel='align pairwise')
        print(magentaStr("\n==> Testing spider - align pairwise with 180 deg. rotation:"))
        self.runAlignment(protFilter, SpiderProtAlignPairwise,
                          objLabel='align pairwise - RT180', cgOption=2)  # RT180
