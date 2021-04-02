# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se)
#                Tapu Shaikh            (shaikh@ceitec.muni.cz)
# *              Grigory Sharov         (gsharov@mrc-lmb.cam.ac.uk)
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

from os.path import join
from glob import glob
import re

import pyworkflow.utils as pwutils
from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params
from pwem.protocols import ProtRefine3D
from pwem.emlib.image import ImageHandler
from pwem.constants import ALIGN_PROJ
from pwem.objects import Volume, FSC

from .. import Plugin
from ..utils import (SpiderDocFile, SpiderDocAliFile,
                     writeScript, runScript)
from ..convert import convertEndian, alignmentToRow
from ..constants import *
from .protocol_base import SpiderProtocol


class SpiderProtRefinement(ProtRefine3D, SpiderProtocol):
    """ Reference-based refinement using SPIDER AP SHC and AP REF commands.

    Iterative refinement improves the accuracy in the determination of orientations.
    This improvement is accomplished by successive use of
    more finely-sampled reference projections.
    
    Two different workflows are suggested: with defocus groups or
    without (gold-standard refinement).
    
    For more information, see:
    [[http://spider.wadsworth.org/spider_doc/spider/docs/techs/recon/mr.html][SPIDER documentation on projection-matching]]
    """
    _label = 'refine 3D'
    _devStatus = PROD

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('protType', params.EnumParam,
                      choices=['with defocus groups', 'gold-standard'],
                      default=GOLD_STD,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Choose refinement type',
                      help='In the first method (with defocus groups), '
                           'ensembles of particles with similar defocus '
                           'values are grouped together. A separate 3D '
                           'reconstruction is computed for each "defocus group." '
                           'CTF-correction is performed when these reconstructions '
                           'are merged into a final reconstruction.\n\n'
                           'In the second, gold-standard method, the CTF-correction '
                           'is applied at the level of windowed particle images. '
                           'More information on Spider '
                           '[[https://spider.wadsworth.org/spider_doc/spider/docs/techs/recon1a/Docs/mr1.html][web-site]]')

        form.addParam('inputParticles', params.PointerParam, 
                      pointerClass='SetOfParticles',
                      pointerCondition='hasCTF',
                      label="Input particles", important=True,
                      help='Select the input particles.\n')  
        
        form.addParam('input3DReference', params.PointerParam,
                      pointerClass='Volume', important=True,
                      label='Initial 3D reference volume:',
                      help='Input 3D reference reconstruction.\n')
        
        form.addParam('numberOfIterations', params.IntParam, default=10,
                      label='Number of iterations:',
                      help='Set the number of iterations. \n\n'
                           'Multi-reference alignment is computationally intensive, '
                           'so rather than use the finest-spaced reference projections immediately, '
                           'we start out with a coarse spacing, and '
                           'then search a restricted range of orientations '
                           'a set of gradually more finely-spaced reference projections.')
        
        form.addParam('alignmentShift', params.IntParam, default=7,
                      label='Shift range',
                      help="Alignments are tested in the range +- this value")
        form.addParam('radius', params.IntParam, default=50,
                      label='Particle radius (px)',
                      help="Radius of the structure (px) used in alignment search\n"
                           "(Default is for ribosome. EDIT as needed.)\n"
                           "This value is used to find radius for last alignment radius.\n")
        form.addParam('winFrac', params.FloatParam, default=0.95,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Projection diameter',
                      help="Fraction of window diameter used in projection\n"
                           " (.95= use 95% window size)\n")
        form.addParam('sphDeconAngle', params.IntParam, default=0,
                      expertLevel=params.LEVEL_ADVANCED,
                      condition='protType == 1',
                      label='Spherical deconvolution angle (deg)',
                      help="Spherical deconvolution angle in degreees (0 == Do not deconvolve)")
        form.addParam('bpType', params.EnumParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      condition='protType == 1',
                      choices=['BP CG', 'BP 3F', 'BP RP', 'BP 3N'],
                      default=BP_3F,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Backprojection method',
                      help="Choose backprojection method (BP CG, BP 3F, BP RP or BP 3N). "
                           "More information on Spider "
                           "[[https://spider.wadsworth.org/spider_doc/spider/docs/bp_overview.html][web-site]]")
        
        form.addParam('smallAngle', params.BooleanParam, default=False,
                      label='Use small angle refinement?',
                      help="In regular, non-small-angle refinement, "
                           "a set of reference projections is computed for all experimental images. \n\n"
                           "In small-angle refinement, "
                           "a set of reference projections is computed for each particle on the fly, "
                           "using the SPIDER command "
                           "[[http://spider.wadsworth.org/spider_doc/spider/docs/man/voras.html][VO RAS]]. ")        
        # GLO [ang-steps]  = '3.3,3.,2.,2.,2.,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5'  ; Angular degree steps
        # GLO [ang-limits] = '0.,0.,15.,8.,6.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.'            ; Angular limits
        form.addParam('angSteps', params.StringParam, default='3.3 3 3x2 1.5',
                      condition='not smallAngle',
                      label='Angular increment',
                      help="This parameter determines how finely spaced the reference will be projected. "
                           "Each value in the list corresponds to the value for a given iteration. "
                           "A value *V* can be repeated *N* times by using the notation *NxV*. "
                           "If more iterations are requested than the number of values specified here, "
                           "the last value will be repeated. \n\n"
                           "For example, in the default *3.3 3 3x2 1.5*, "
                           "the increment will be: iter 1 - 3.3 degrees, iter 2 - 3 degrees, iter 3,4,5 - 2 degrees, "
                           "and iterations from the sixth onward will use 1.5 degrees.")
        form.addParam('angLimits', params.StringParam, default='2x0 15 8 6 5',
                      condition='not smallAngle',
                      label='Angular range',
                      help="This parameter determines the range of reference projections that will be searched. "
                           "A value of *0* corresponds to an unrestricted search. "
                           "A small value will result in a faster alignment, but may not find the correct orientation. "
                           "A value *V* can be repeated *N* times by using the notation *NxV*. "
                           "If more iterations are requested than the number of values specified here, "
                           "the last value will be repeated. \n\n"
                           "For example, in the default *2x0 15 8 6 5*, "
                           "the increment will be: iter 1,2 - unrestricted, iter 3 - 15 degrees, iter 4 - 8 degrees, "
                           "iter 5 - 6 degrees and iterations from the sixth onward will use 5 degrees.")
        form.addParam('angStepSm', params.FloatParam, default=0.5,
                      condition='smallAngle',
                      label='Angular increment',
                      help="This parameter determines how finely spaced the reference will be projected.")          
        form.addParam('thetaRange', params.FloatParam, default=2.0,
                      condition='smallAngle',
                      label='Angular range ',
                      help="This parameter determines the range of reference projections that will be searched.")          

        form.addParallelSection(threads=4, mpi=0)
    
    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):        
        # Create new stacks and selfiles per groups
        self._insertFunctionStep('convertInputStep',
                                 self.inputParticles.get().getObjId())

        self._insertFunctionStep('runScriptStep', 'refine.pam')
                
        self._insertFunctionStep('createOutputStep')
    
    # --------------------------- STEPS functions -----------------------------
    
    def convertInputStep(self, particlesId):
        """ Convert all needed inputs before running the refinement script. """
        partSet = self.inputParticles.get()
        protType = self.protType.get()

        self._writeParamsFile(partSet)
        self._writeGroupFiles(partSet, protType)
        
        # Convert the input volume
        volPath = self._getExtraPath('ref_vol.vol')
        ImageHandler().convert(self.input3DReference.get(), volPath)
        pwutils.moveFile(volPath, volPath.replace('.vol', '.stk'))
        
        self._writeRefinementScripts(protType)
                
    def _writeRefinementScripts(self, protType):
        """ Write the needed scripts to run refinement
        and substitute some values.
        """
        refPath = self._getExtraPath('Refinement')
        pwutils.makePath(refPath)
        
        def path(p):
            """ Escape path with '' and add ../ """
            return "'%s'" % join('..', p)
        
        def script(name, paramsDict={}, protType=protType):
            if protType == DEF_GROUPS:
                dirName = 'defocus-groups'
            else:
                dirName = 'no-defocus-groups'

            outputScript = join(refPath, name)
            writeScript(Plugin.getScript('projmatch', 'Refinement', dirName, name),
                        outputScript, paramsDict)
            
        nIter = self.numberOfIterations.get()
        
        def getListStr(valueStr):
            return "'%s'" % ','.join(pwutils.getListFromValues(valueStr, nIter))

        diam = int(self.radius.get() * 2 * self.inputParticles.get().getSamplingRate())
        params = {'[alignsh]': self.alignmentShift.get(),
                  # shrange was renamed to alignsh in new versions
                  '[shrange]': self.alignmentShift.get(),
                  '[iter-end]': self.numberOfIterations.get(),
                  '[diam]': diam,
                  '[win-frac]': self.winFrac.get(),
                  # '[converg]': self.convergence.get(),
                  '[small-ang]': '1' if self.smallAngle else '0',
                  '[ang-steps]': getListStr(self.angSteps.get()),
                  '[ang-limits]': getListStr(self.angLimits.get()),
                  '[ang-step-sm]': "'(%0.2f)'" % self.angStepSm.get(),
                  '[theta-range]': "'(%0.2f)'" % self.thetaRange.get(),
                  
                  '[vol_orig]': path('ref_vol'),
                  '[sel_group_orig]': path('sel_group'),
                  '[sel_particles_orig]': path('group{***[grp]}_selfile'),
                  '[group_align_orig]': path('group{***[grp]}_align'),
                  '[unaligned_images_orig]': path('group{***[grp]}_stack'),
                  '[out_align]': path('stack_alignment')
                  }

        if self.protType == GOLD_STD:
            params.update({'sphdecon': self.sphDeconAngle.get(),
                           'bp-type': self.bpType.get() + 1})

        script('refine_settings.pam', params)
        if protType == DEF_GROUPS:
            scriptList = ['refine', 'prepare', 'grploop', 'mergegroups',
                          'enhance', 'endmerge', 'smangloop', 'endrefine']
        else:
            scriptList = ['refine', 'refine-setrefangles',
                          'refine-prjrefs', 'refine-loop', 'refine-smangloop',
                          'refine-bp', 'merge-fsc-filt', 'sphdecon', 'enhance',
                          'show-r2']

        for s in scriptList:
            script('%s.pam' % s)
        
    def _writeParamsFile(self, partSet):
        acq = partSet.getAcquisition()
        params = {'datetime': 'now',
                  'pixelSize': partSet.getSamplingRate(),
                  'voltage': acq.getVoltage(),
                  'sphericalAberration': acq.getSphericalAberration(),
                  'windowSize': partSet.getDimensions()[0],
                  'nummps': self.numberOfThreads.get()
                  }
        
        paramFile = open(self._getExtraPath('params.stk'), 'w+')
        paramFile.write("""
 ;spi/dat  Generated by Scipion on %(datetime)s  params.stk
    5 1    %(pixelSize)f           ; pixel size (A)
    6 1    %(voltage)f             ; electron energy (kV)
    7 1    %(sphericalAberration)f ; spherical aberration (mm)
   17 1    %(windowSize)f          ; window size (pixels)
   18 1    %(nummps)d              ; number of threads to use
                        """ % params)
        paramFile.close()        
        
    def _writeGroupFiles(self, partSet, protType):
        """Write files that are needed by each group:
        - stack
        - selfile
        - docfile
        """
        ih = ImageHandler()
        # Keep a dict with all groups found in particles
        groupDict = {}
        template = self._getExtraPath('group%03d_%s.stk')

        if protType == DEF_GROUPS:
            for part in partSet:
                defocusGroup = self._getDefocusGroup(part)

                if defocusGroup not in groupDict:
                    groupInfo = DefocusGroupInfo(defocusGroup, template, ih)
                    groupDict[defocusGroup] = groupInfo
                else:
                    groupInfo = groupDict[defocusGroup]

                groupInfo.addParticle(part)
        else:
            numProcs = self.numberOfMpi.get()
            # explicitly set a minimum of 2 groups
            # GS: maybe not necessary?
            numGroups = 2 if numProcs < 3 else numProcs
            d, r = divmod(len(partSet), numGroups)
            numParts = d + 1 if r > 0 else d
            groupId = 1

            for part in partSet:
                if groupId not in groupDict:
                    groupInfo = DefocusGroupInfo(groupId, template, ih)
                    groupDict[groupId] = groupInfo
                else:
                    groupInfo = groupDict[groupId]
                    groupSize = groupDict[groupId].counter
                    if groupSize > numParts:
                        groupId += 1
                        groupInfo = DefocusGroupInfo(groupId, template, ih)
                        groupDict[groupId] = groupInfo

                groupInfo.addParticle(part)

        # Write the docfile with the group information
        # like the number of particles (and the defocus)
        groupsDoc = SpiderDocFile(self._getExtraPath('sel_group.stk'), 'w+')
        for gi in groupDict.values():
            if protType == DEF_GROUPS:
                groupsDoc.writeValues(gi.number, gi.counter, gi.defocus)
            else:
                groupsDoc.writeValues(gi.number, gi.counter)
            # Convert the endianness of the stack
            convertEndian(gi.stackfile, gi.counter)
            # Close each group docfile
            gi.close()

        groupsDoc.close()
        
    def runScriptStep(self, script):
        """ Just run the script that was generated in convertInputStep. """
        refPath = self._getExtraPath('Refinement')
        runScript(script, 'pam/stk', program=Plugin.getProgram(),
                  nummpis=1, cwd=refPath, log=self._log)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        vol = Volume()
        lastIter = self._getLastIterNumber()
        if self.protType == GOLD_STD:
            vol.setFileName(self._getExtraPath('Refinement/final/vol_%02d.stk' % lastIter))
            half1 = self._getExtraPath('Refinement/final/vol_%02d_s1.stk' % lastIter)
            half2 = self._getExtraPath('Refinement/final/vol_%02d_s2.stk' % lastIter)
        else:
            vol.setFileName(self._getExtraPath('Refinement/final/bpr%02d.stk' % lastIter))
            half1 = self._getExtraPath('Refinement/final/bpr%02d_sub1.stk' % lastIter)
            half2 = self._getExtraPath('Refinement/final/bpr%02d_sub2.stk' % lastIter)
        vol.setSamplingRate(imgSet.getSamplingRate())
        vol.setHalfMaps([half1, half2])

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromDoc(outImgSet)

        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputParticles, vol)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles, outImgSet)

        fsc = FSC(objLabel=self.getRunName())
        resolution, fscData = self._getFscData(it=lastIter)
        fsc.setData(resolution, fscData)

        self._defineOutputs(outputFSC=fsc)
        self._defineSourceRelation(vol, fsc)
    
    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        if self.smallAngle and not self.inputParticles.get().hasAlignmentProj():
            errors.append('*Small angle* option can only be used if '
                          'the particles have angular assignment.')
        return errors
        
    def _warnings(self):
        warns = []
        if not self.smallAngle:
            niter = self.numberOfIterations.get()
            if niter < len(pwutils.getListFromValues(self.angSteps.get())):
                warns.append('*Angular steps* have more values than iterations')           
    
        return warns
    
    def _summary(self):
        summary = list()
        summary.append('Number of iterations: *%s*' % self.numberOfIterations)
        
        if self.smallAngle:
            summary.append('Small-angle refinement: ')
            summary.append('    Angular increment: *%s* degrees' % self.angStepSm)
            summary.append('    Angular range: *%s* degrees' % self.thetaRange)
        else:
            summary.append('Angular increments: *%s*' % self.angSteps)
            summary.append('Angular range: *%s*' % self.angLimits)

        diam = int(self.radius.get() * 2 * self.inputParticles.get().getSamplingRate())
        summary.append('Particle diameter: *%d* Angstroms' % diam)
        summary.append('Shift range: *%s* pixels' % self.alignmentShift)
        # summary.append('Projection diameter: *%s* of window size' % self.winFrac)

        return summary

    def _citations(self):
        return ['Penczek1992']
    
    def _methods(self):
        msg = "\nInput particles %s " % self.getObjectTag('inputParticles')
        msg += "were subjected to 3D refinement by projection matching ([Penczek1992]) "
        msg += "for %s iterations " % self.numberOfIterations
        msg += "using %s as an initial reference. " % self.getObjectTag('input3DReference')
        if self.protType == GOLD_STD:
            msg += "\nUsing gold-standard refinement procedure."
        else:
            msg += "\nUsing defocus group-based procedure."
        
        return [msg]
        
    # --------------------------- UTILS functions -----------------------------
    def _getDefocusGroup(self, img):
        return img.getMicId()

    def _getLastIterNumber(self):
        """ Return the list of iteration files, give the iterTemplate. """
        result = None
        if self.protType == DEF_GROUPS:
            template = self._getExtraPath('Refinement/final/bpr??.stk')
        else:
            template = self._getExtraPath('Refinement/final/vol_??.stk')
        files = sorted(glob(template))
        if files:
            f = files[-1]
            s = re.compile('(\d{2}).stk').search(f)
            if s:
                result = int(s.group(1))
        return result

    def _fillDataFromDoc(self, imgSet):
        outDocFn = SpiderDocAliFile(self._getExtraPath('stack_alignment.stk'))
        imgSet.setAlignmentProj()
        initPartSet = self.inputParticles.get()
        partIter = iter(initPartSet.iterItems(orderBy=['id'], direction='ASC'))
        imgSet.copyItems(partIter,
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=iter(outDocFn))

    def _createItemMatrix(self, item, row):
        from ..convert import createItemMatrix
        from pwem.constants import ALIGN_PROJ
        createItemMatrix(item, row, align=ALIGN_PROJ)

    def _getFscData(self, it):
        if self.protType == GOLD_STD:  # gold std
            fn = self._getExtraPath("Refinement/final/fscdoc_m_%02d.stk" % it)
        else:  # def groups
            fn = self._getExtraPath("Refinement/final/ofscdoc_%02d.stk" % it)

        resolution = []
        fscData = []
        fscDoc = SpiderDocFile(fn)
        for values in fscDoc:
            resolution.append(1 / values[1])
            fscData.append(values[2])

        return resolution, fscData


class DefocusGroupInfo:
    """ Helper class to store some information about 
    defocus groups like the number of particles
    or the docfile to be generated.
    """
    def __init__(self, defocusGroup, template, ih):
        self.ih = ih
        self.number = defocusGroup
        self.selfile = template % (defocusGroup, 'selfile')
        self.docfile = template % (defocusGroup, 'align')
        self.stackfile = template % (defocusGroup, 'stack')
        self.counter = 0  # number of particles in this group

        self.sel = SpiderDocFile(self.selfile, 'w+')
        self.doc = SpiderDocFile(self.docfile, 'w+')
        self.doc.writeComment(self.docfile)
        header = ['KEY', 'PSI', 'THE', 'PHI', 'REF#', 'EXP#', 'CUM.{ROT',
                  'SX', 'SY}', 'NPROJ', 'DIFF', 'CCROT', 'ROT', 'SX', 'SY', 'MIR-CC']
        self.doc.writeHeader(header)

    def addParticle(self, img):
        self.counter += 1
        if self.counter == 1:
            ctf = img.getCTF()
            self.defocus = (ctf.getDefocusU() + ctf.getDefocusV()) / 2.
        self.ih.convert(img, (self.counter, self.stackfile))
        self.sel.writeValues(self.counter)
        alignRow = {ANGLE_PSI: 0.,
                    ANGLE_THE: 0.,
                    ANGLE_PHI: 0.,
                    SHIFTX: 0.,
                    SHIFTY: 0.}
        alignment = img.getTransform()
        
        if alignment is not None:
            alignmentToRow(alignment, alignRow, ALIGN_PROJ)
            
        values = [0.00, alignRow[ANGLE_THE], alignRow[ANGLE_PHI], 
                  0.00, self.counter, 
                  alignRow[ANGLE_PSI], alignRow[SHIFTX], alignRow[SHIFTY], 
                  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.0]
        self.doc.writeValues(*values)
        
    def close(self):
        self.sel.close()
        self.doc.close()
