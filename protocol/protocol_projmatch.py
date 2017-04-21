# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
#                Tapu Shaikh            (shaikh@ceitec.muni.cz)
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

from os.path import join
from glob import glob
import re

import pyworkflow.utils as pwutils
import pyworkflow.em as em
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtRefine3D
from pyworkflow.em.constants import ALIGN_PROJ
from pyworkflow.em.data import Volume

from ..spider import SpiderDocFile, writeScript, getScript, runScript, getVersion
from ..Spiderutils import nowisthetime
from ..convert import ANGLE_PHI, ANGLE_PSI, ANGLE_THE, SHIFTX, SHIFTY, convertEndian, alignmentToRow
from protocol_base import SpiderProtocol

# Protocol type
DEF_GROUPS = 0
GOLD_STD = 1

# Backprojection method
BP_CG = 0
BP_3F = 1
BP_RP = 2
BP_3N = 3

class SpiderProtRefinement(ProtRefine3D, SpiderProtocol):
    """ Iterative reference-based refinement of particles orientations, 
    based on the Spider AP SHC and AP REF programs.
    
    Two different workflows are suggested: with defocus groups or
    without (gold-standard refinement).
    
    Iterative refinement improves the accuracy in the determination of orientations.
    This improvement is accomplished by successive use of
    more finely-sampled reference projections.
    
    For more information, see:
    [[http://spider.wadsworth.org/spider_doc/spider/docs/techs/recon/mr.html]
    [SPIDER documentation on projection-matching]]
    
    """
    _label = 'refinement'

    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('protType', params.EnumParam,
                      choices=['with defocus groups','gold-standard'],
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
                           'This method presents certain advantages. First, it '
                           'circumvents one of the approximations when using '
                           'defocus groups, namely that all particles in a defocus '
                           'group follow the same CTF profile. At high resolution, '
                           'where the CTF oscillates more rapidly, this assumption '
                           'may not hold. Second, parallelization can be more '
                           'efficient, since groups can be of identical size, '
                           'independent of the number of particles at each defocus. '
                           'Third, particles from what would be sparsely populated '
                           'defocus groups need not be thrown out.\n\nHowever the '
                           'strategy of using defocus groups may have some advantages. '
                           'First is that it can readily account for the non-uniform '
                           'distribution of signal-to-noise in projection data '
                           '(Penczek, 2012). Second, we find that reconstructions '
                           'using particle-level CTF-correction sometimes show '
                           'artifacts when using iterative backprojection methods, '
                           'such as *BP RP* or *BP CG*, whereas the use of defocus '
                           'groups does not present such limitations.')

        form.addParam('inputParticles', params.PointerParam, 
                      pointerClass='SetOfParticles',
                      pointerCondition='hasCTF',
                      label="Input particles", important=True,
                      help='Select the input particles.\n')  
        
        form.addParam('input3DReference', params.PointerParam,
                      pointerClass='Volume', 
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
        form.addParam('convergence', params.FloatParam, default=0.05,
                      expertLevel=params.LEVEL_ADVANCED,
                      condition='protType == 0 and not _protGoldStdIsSupported',
                      label='Convergence criterion fraction',
                      help="Refinement has converged when this fraction of all images move < 1.5 * stepsize.")
        form.addParam('sphDeconAngle', params.IntParam, default=0,
                      expertLevel=params.LEVEL_ADVANCED,
                      condition='_protGoldStdIsSupported and protType == 1',
                      label='Spherical deconvolution angle (deg)',
                      help="Spherical deconvolution angle in degreees (0 == Do not deconvolve)")
        form.addParam('bpType', params.EnumParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      condition='_protGoldStdIsSupported and protType == 1',
                      choices=['BP CG', 'BP 3F', 'BP RP', 'BP 3N'],
                      default=BP_3F,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Backprojection method',
                      help="Choose backprojection method (BP CG, BP 3F, BP RP or BP 3N). "
                           "More information on Spider "
                           "[https://spider.wadsworth.org/spider_doc/spider/docs/bp_overview.html][web-site]]")
        
        form.addParam('smallAngle', params.BooleanParam, default=False,
                      label='Use small angle refinement?',
                      help="In regular, non-small-angle refinement, "
                           "a set of reference projections is computed for all experimental images. \n\n"
                           "In small-angle refinement, "
                           "a set of reference projections is computed for each particle on the fly, "
                           "using the SPIDER command "
                           "[[http://spider.wadsworth.org/spider_doc/spider/docs/man/voras.html][VO RAS]]. ")        
        #GLO [ang-steps]  = '3.3,3.,2.,2.,2.,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5'  ; Angular degree steps   
        #GLO [ang-limits] = '0.,0.,15.,8.,6.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.'            ; Angular limits
        form.addParam('angSteps', params.StringParam, default='3.3 3 3x2 1.5',
                      condition='not smallAngle',
                      label='Angular increment',
                      help="This parameter determines how finely spaced the reference will be projected. "
                           "Each value in the list corresponds to the value for a given iteration. "
                           "A value *V* can be repeated *N* times by using the notation *NxV*. "
                           "If more iterations are requested than the number of values specified here, "
                           "the last value will be repeated. \n\n"
                           "For example, in the default *3.3 3 3x2 1.5*, "
                           "the value of 2 degrees will be repeated three times, "
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
                           "and iterations from the sixth onward will use 5 degrees, "
                           "an unrestricted search will be performed twice.")
        form.addParam('angStepSm', params.FloatParam, default=0.5,
                      condition='smallAngle',
                      label='Angular increment',
                      help="This parameter determines how finely spaced the reference will be projected.")          
        form.addParam('thetaRange', params.FloatParam, default=2.0,
                      condition='smallAngle',
                      label='Angular range ',
                      help="This parameter determines the range of reference projections that will be searched.")          

        form.addParallelSection(threads=4, mpi=0)
    
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self):        
        # Create new stacks and selfiles per defocus groups
        self._insertFunctionStep('convertInputStep', self.inputParticles.get().getObjId())

        self._insertFunctionStep('runScriptStep', 'refine.pam')
                
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions --------------------------------------------
    
    def convertInputStep(self, particlesId):
        """ Convert all needed inputs before running the refinement script. """
        partSet = self.inputParticles.get()

        self._writeParamsFile(partSet)
        self._getLastIterNumber()
        self._writeGroupFiles(partSet)
        
        # Convert the input volume
        volPath = self._getExtraPath('vol01.vol')
        em.ImageHandler().convert(self.input3DReference.get(), volPath)
        pwutils.moveFile(volPath, volPath.replace('.vol', '.stk'))
        
        self._writeRefinementScripts()
                
    def _writeRefinementScripts(self):
        """ Write the needed scripts to run refinement
        and substitute some values.
        """
        
        refPath = self._getExtraPath('Refinement')
        pwutils.makePath(refPath)
        
        def path(p):
            """ Escape path with '' and add ../ """
            return "'%s'" % join('..', p)
        
        def script(name, paramsDict={}, protType=DEF_GROUPS):
            if protType == DEF_GROUPS:
                dirName = 'defocus-groups'
            else:
                dirName = 'no-defocus-groups'

            outputScript=join(refPath, name)
            writeScript(getScript('projmatch', 'Refinement', dirName, name), outputScript, paramsDict)
            
        nIter = self.numberOfIterations.get()
        
        def getListStr(valueStr):
            return "'%s'" % ','.join(pwutils.getListFromValues(valueStr, nIter))

        diam = int(self.radius.get() * 2 * self.inputParticles.get().getSamplingRate())
        params = {'[alignsh]': self.alignmentShift.get(),  # shrange is renamed to alignsh in new versions
                  '[shrange]': self.alignmentShift.get(),
                  '[iter-end]': self.numberOfIterations.get(),
                  '[diam]': diam,
                  '[win-frac]': self.winFrac.get(),
                  '[converg]': self.convergence.get(),
                  '[small-ang]': '1' if self.smallAngle else '0',
                  '[ang-steps]': getListStr(self.angSteps.get()),
                  '[ang-limits]': getListStr(self.angLimits.get()),
                  '[ang-step-sm]': '(%0.2f)' % self.angStepSm.get(),
                  '[theta-range]': '(%0.2f)' % self.thetaRange.get(),
                  
                  '[vol_orig]': path('vol001'),
                  '[sel_group_orig]': path('sel_group'),
                  '[sel_particles_orig]': path('group{***[grp]}_selfile'),
                  '[group_align_orig]': path('group{***[grp]}_align'),
                  '[unaligned_images_orig]': path('group{***[grp]}_stack')
                  }

        if self._protGoldStdIsSupported() and self.protType == GOLD_STD:
            params.update({'sphdecon': self.sphDeconAngle.get(),
                           'bp-type': self.bpType.get() + 1})

        script('refine_settings.pam', params)
        for s in ['refine', 'prepare', 'grploop', 'mergegroups', 
                  'enhance', 'endmerge', 'smangloop', 'endrefine']:
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
        
    def _writeGroupFiles(self, partSet):
        """Write files that are needed by each defocus group:
        - stack
        - selfile
        - docfile
        """
        ih = em.ImageHandler()
        # Keep a dict with all groups found in particles
        groupDict = {}
        template = self._getExtraPath('group%03d_%s.stk')
        
        for part in partSet:
            defocusGroup = self._getDefocusGroup(part)
            
            if defocusGroup not in groupDict:
                groupInfo = DefocusGroupInfo(defocusGroup, template, ih)
                groupDict[defocusGroup] = groupInfo
            else:
                groupInfo = groupDict[defocusGroup]
                
            groupInfo.addParticle(part)

        # Write the docfile with the defocus groups information
        # like the number of particles and the defocus
        groupsDoc = SpiderDocFile(self._getExtraPath('sel_group.stk'), 'w+')
        for gi in groupDict.values():
            groupsDoc.writeValues(gi.number, gi.counter, gi.defocus)
            # Convert the endianness of the stack
            convertEndian(gi.stackfile, gi.counter)
            # Close each group docfile
            gi.close()

        groupsDoc.close()
        
    def runScriptStep(self, script):
        """ Just run the script that was generated in convertInputStep. """
        refPath = self._getExtraPath('Refinement')
        runScript(script, 'pam/stk', cwd=refPath, log=self._log)
        
    def projectStep(self, volumeId):
        pass
     
    def alignStep(self):
        """Create new stacks and selfiles per defocus groups """
        pass
    
    def reconstructStep(self):
        pass   
    
    def mergeStep(self):
        pass    
    
    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        vol = Volume()
        vol.setFileName(self._getExtraPath('Refinement/final/bpr%02d.stk' % self._getLastIterNumber()))
        vol.setSamplingRate(imgSet.getSamplingRate())

        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputParticles, vol)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        if self.protType == GOLD_STD and getVersion() == '21.03':
            errors.append('Gold-standard refinement is supported only '
                          'for newer Spider versions. Please update your installation.')
        if self.smallAngle:
            if not self.inputParticles.get().hasAlignmentProj():
                errors.append('*Small angle* option can only be used if '
                              'the particles have angular assignment.')
        return errors
        
    def _warnings(self):
        print
        print "protocol_projmatch._warning"
        print
        warns = []
        if not self.smallAngle:
            niter = self.numberOfIterations.get()
            if niter < len(pwutils.getListFromValues(self.angSteps.get())):
                warns.append('*Angular steps* have more values than iterations')           
    
        return warns
    
    def _summary(self):
        summary = []
        summary.append('Number of iterations: *%s*' % self.numberOfIterations)
        
        if self.smallAngle:
            summary.append('Small-angle refinement: ')
            summary.append('    Angular increment: *%s* degrees' % self.angStepSm)
            summary.append('    Angular range: *%s* degrees' % self.thetaRange)
        else:
            summary.append('Angular increments: *%s*' % self.angSteps)
            summary.append('Angular range: *%s*' % self.angLimits)

        diam = int(self.radius.get() * 2 * self.inputParticles.get().getSamplingRate())
        summary.append('Particle diameter: *%s* Angstroms' % diam)
        summary.append('Shift range: *%s* pixels' % self.alignmentShift)
        summary.append('Projection diameter: *%s* of window size' % self.winFrac)

        return summary
        
    def _citations(self):
        return ['Penczek1992']
    
    def _methods(self):
        msg  = "\nInput particles %s " % self.getObjectTag('inputParticles')
        msg += "were subjected to refinement of orientations ([Penczek1992]) "
        msg += "for %s iterations " % self.numberOfIterations
        msg += "using %s as an initial reference. " % self.getObjectTag('input3DReference')
        
        return [msg]
        
    #--------------------------- UTILS functions --------------------------------------------
    def _getDefocusGroup(self, img):
        return img.getMicId()

    def _getLastIterNumber(self):
        """ Return the list of iteration files, give the iterTemplate. """
        result = None
        template = self._getExtraPath('Refinement/final/bpr??.stk')
        files = sorted(glob(template))
        if files:
            f = files[-1]
            s = re.compile('bpr(\d{2})').search(f)
            if s:
                result = int(s.group(1))
        return result

    def _protGoldStdIsSupported(self):
        return True if getVersion() != '21.03' else False

    
class DefocusGroupInfo():
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
        date, time, _ = nowisthetime()
        self.doc.writeComment("spi/dat   Generated by Scipion on %s AT %s" % (date, time))
        self.doc.writeComment("  KEY       PSI,    THE,    PHI,   REF#,    EXP#,  CUM.{ROT,   SX,    SY},  NPROJ,   DIFF,      CCROT,    ROT,     SX,     SY,   MIR-CC")
        
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
