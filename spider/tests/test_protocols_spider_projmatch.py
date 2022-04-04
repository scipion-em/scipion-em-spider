# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from pwem.protocols import ProtImportParticles, ProtImportVolumes

from ..protocols import SpiderProtRefinement


class TestSpiderBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')
        cls.vol = cls.dsRelion.getFile('import/case2/volume.mrc')

    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         filesPath=pattern,
                                         samplingRate=samplingRate,
                                         checkStack=checkStack)
        cls.launchProtocol(cls.protImport)
        return cls.protImport

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import volumes protocol. """
        cls.protImportVols = cls.newProtocol(ProtImportVolumes,
                                             filesPath=pattern,
                                             samplingRate=samplingRate)
        cls.launchProtocol(cls.protImportVols)
        return cls.protImportVols


class TestSpiderRefinement(TestSpiderBase):
    def test_ProjMatchSpider(self):
        print(magentaStr("\n==> Importing data - particles:"))
        protImportPart = self.newProtocol(ProtImportParticles,
                                          objLabel='from scipion (to-reconstruct)',
                                          importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                          sqliteFile=self.dsRelion.getFile('import/case2/particles.sqlite'),
                                          magnification=10000,
                                          samplingRate=7.08)
        self.launchProtocol(protImportPart)
        self.assertIsNotNone(protImportPart.getFiles(),
                             "There was a problem with the import")

        print(magentaStr("\n==> Importing data - volume:"))
        protImportVol = self.newProtocol(ProtImportVolumes,
                                         filesPath=self.vol,
                                         samplingRate=7.08)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(),
                             "There was a problem with the import")

        print(magentaStr("\n==> Testing spider - refinement:"))
        protRefine = self.newProtocol(SpiderProtRefinement,
                                      numberOfIterations=2,
                                      alignmentShift=2,
                                      radius=26,
                                      smallAngle=True,
                                      angStepSm=1,
                                      thetaRange=2)
        protRefine.inputParticles.set(protImportPart.outputParticles)
        protRefine.input3DReference.set(protImportVol.outputVolume)
        self.launchProtocol(protRefine)
        self.assertIsNotNone(protRefine.outputVolume,
                             "There was a problem with Spider refinement protocol")
