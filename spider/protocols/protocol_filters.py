# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se)
# *              Tapu Shaikh            (shaikh@ceitec.muni.cz)
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

from enum import Enum

from pwem.protocols import ProtFilterParticles
from pwem.objects import SetOfParticles
from pyworkflow.protocol.params import (EnumParam, BooleanParam,
                                        DigFreqParam, FloatParam)
from pyworkflow.utils.path import removeBaseExt
from pyworkflow.constants import PROD

from ..constants import (FILTER_BUTTERWORTH, FILTER_FERMI,
                         FILTER_LOWPASS, FILTER_SPACE_REAL)
from ..utils import SpiderShell
from .protocol_base import SpiderProtocol


class outputs(Enum):
    outputParticles = SetOfParticles

      
class SpiderProtFilter(ProtFilterParticles, SpiderProtocol):
    # the URL doesn't load
    """
    Applies Fourier-space filtering operations to particle images or volumes
    using several classical frequency-domain filters commonly employed in
    cryo-EM image processing. The protocol is designed to enhance structural
    signal, suppress noise, or emphasize specific spatial frequency ranges in
    preparation for downstream reconstruction, classification, alignment, or
    visualization tasks. More info:
    https://spider.wadsworth.org/spider_doc/spider/docs/man/fq.html

    AI Generated:

    Filter Particles (SpiderProtFilter) - User Manual
        Overview

        The Filter Particles protocol applies frequency-domain filters to
        cryo-EM particle images or volumetric data in order to control the
        balance between signal preservation and noise suppression. Filtering is
        one of the most common preprocessing operations in structural biology
        because experimental images typically contain substantial levels of
        high-frequency noise, low-frequency background variation, or unwanted
        frequency components introduced during acquisition and reconstruction.

        In practical cryo-EM workflows, filtering is frequently used before
        particle alignment, classification, refinement, or visualization.
        Depending on the scientific objective, the protocol may either suppress
        noisy high-frequency information, remove low-frequency background
        gradients, or isolate a specific frequency band that contains relevant
        structural detail. Proper filtering often improves algorithmic
        stability and visual interpretability, although excessive filtering may
        remove biologically meaningful information.

        Inputs and General Workflow

        The protocol accepts a set of particles or images that will be
        transformed through Fourier filtering operations. The resulting output
        preserves the original particle organization while generating a new
        filtered dataset suitable for subsequent analysis.

        From a biological perspective, filtering should be understood as a
        signal-conditioning step rather than a reconstruction method. The goal
        is not to create new structural information but to improve the
        interpretability of information already present in the data. Care must
        therefore be taken to avoid introducing misleading visual features or
        suppressing weak but biologically relevant signals.

        Filter Types and Their Biological Meaning

        Several filtering models are available, each emphasizing different
        frequency behaviors and suited to different biological scenarios.

        The Top-hat filter performs an abrupt truncation of frequencies beyond
        a selected cutoff. This approach is computationally simple and useful
        for exploratory analysis or aggressive denoising, but the sharp
        transition may introduce ringing artifacts near strong density
        boundaries. Biological users should therefore apply this filter
        cautiously when interpreting fine structural details.

        The Gaussian filter produces a smooth attenuation of frequencies and is
        one of the most commonly used filtering strategies in cryo-EM. Because
        the transition is gradual, Gaussian filtering tends to preserve overall
        structural continuity while reducing high-frequency noise. This option
        is often preferred for visualization, initial preprocessing, or gentle
        denoising workflows.

        The Fermi filter provides an intermediate behavior between abrupt and
        smooth filtering. By adjusting the temperature parameter, the user can
        control how sharply the transition occurs around the cutoff frequency.
        This flexibility is particularly useful when balancing noise reduction
        against preservation of intermediate-resolution features.

        The Butterworth filter is widely used in signal processing because it
        allows controlled frequency attenuation with adjustable steepness. The
        order parameter determines how rapidly frequencies are suppressed
        beyond the cutoff region. Lower orders provide smoother transitions,
        while higher orders behave more aggressively. In biological practice,
        moderate Butterworth filters are often effective for reducing noise
        while preserving interpretable density boundaries.

        The Raised cosine filter creates a smooth transition between retained
        and suppressed frequencies within a specified frequency interval. This
        approach is especially useful when the user wishes to isolate or
        emphasize a frequency band without introducing abrupt Fourier-space
        discontinuities.

        Low-Pass and High-Pass Filtering

        The protocol supports both low-pass and high-pass filtering modes.
        Understanding the biological implications of these modes is essential
        for correct interpretation.

        Low-pass filtering suppresses high-frequency components and is commonly
        used to reduce noise. This operation smooths the data and highlights
        large-scale structural organization. In cryo-EM workflows, low-pass
        filtering is frequently applied during early refinement stages,
        visualization, or preparation of low-resolution references. Excessive
        low-pass filtering, however, may remove secondary structure features
        or blur biologically meaningful conformational differences.

        High-pass filtering suppresses low-frequency information and enhances
        local contrast or fine structural details. This may improve visibility
        of edges or local features but can also amplify noise if applied too
        aggressively. Biological users should employ high-pass filtering
        carefully, particularly when analyzing weak densities or flexible
        regions.

        Frequency Selection and Resolution Interpretation

        Frequency parameters define which spatial frequencies are retained or
        attenuated. In cryo-EM terms, these frequencies are directly related
        to structural resolution.

        Lower spatial frequencies correspond to broad global shapes and overall
        molecular architecture, while higher frequencies contain finer
        structural information such as secondary structure elements or local
        side-chain detail. Selecting an appropriate cutoff therefore depends on
        the biological question being addressed.

        For example, strong low-pass filtering may be appropriate when
        inspecting overall domain organization or particle orientation, whereas
        milder filtering is preferred when evaluating secondary structure
        quality or map interpretability.

        Padding and Boundary Effects

        The protocol optionally applies padding during filtering operations.
        Padding extends the image boundaries before Fourier transformation,
        reducing edge discontinuities and minimizing artifacts introduced by
        periodic boundary assumptions in Fourier processing.

        From a practical perspective, padding is generally recommended because
        it improves boundary behavior and reduces artificial ringing near the
        particle edges. This becomes particularly important for particles that
        occupy a large fraction of the image box or contain strong density
        gradients near the borders.

        Disabling padding may reduce computational overhead slightly, but it
        can increase the risk of edge artifacts that interfere with alignment,
        classification, or interpretation.

        Outputs and Their Interpretation

        The protocol produces a filtered particle dataset that preserves the
        original metadata and organizational structure while replacing the
        image content with filtered versions. The resulting particles can be
        used directly in downstream cryo-EM workflows including alignment,
        classification, reconstruction, or visualization.

        Biologically, filtered particles should always be interpreted in the
        context of the chosen frequency parameters. Features removed by
        filtering are not necessarily absent from the specimen itself but may
        simply have been suppressed computationally. For this reason,
        publication-quality interpretation should ideally compare filtered and
        unfiltered representations whenever possible.

        Practical Recommendations

        In routine cryo-EM processing, Gaussian or moderate Butterworth
        low-pass filtering provides a good starting point for reducing noise
        while preserving structural continuity. Padding should generally remain
        enabled unless there is a strong computational reason to disable it.

        Aggressive filtering strategies may improve visual appearance but can
        distort biological interpretation if overused. When evaluating flexible
        complexes, weak ligand densities, or heterogeneous conformations,
        conservative filtering is usually preferable.

        High-pass filtering is best reserved for specialized applications such
        as contrast enhancement or feature detection rather than routine
        structural interpretation.

        Final Perspective

        Fourier filtering is a foundational operation in cryo-EM image
        analysis because it directly shapes the balance between signal and
        noise. Although mathematically straightforward, filtering decisions can
        strongly influence biological interpretation, alignment stability, and
        downstream reconstruction quality. Careful selection of filter type,
        cutoff frequency, and transition behavior is therefore essential for
        producing reliable and biologically meaningful results.
    """
    _label = 'filter particles'
    _devStatus = PROD
    _possibleOutputs = outputs
    
    def __init__(self, **kwargs):
        ProtFilterParticles.__init__(self, **kwargs)
        SpiderProtocol.__init__(self, **kwargs)
        # To avoid showing MPI box due to duplicated init
        self.allowMpi = False

        self._op = "FQ"
        self._params = {'ext': 'stk', 
                        'particles': 'particles_filtered',
                        'particlesSel': 'particles_filtered_sel'}

    # --------------------------- DEFINE param functions ----------------------
    def _defineProcessParams(self, form):
        form.addParam('filterType', EnumParam,
                      choices=['Top-hat', 'Gaussian', 'Fermi',
                               'Butterworth', 'Raised cosine'],
                      label="Filter type", default=FILTER_BUTTERWORTH,
                      help="""Select what type of filter do you want to apply.
                      
*Top-hat*: Filter is a "top-hat" function 
that truncates the Fourier transform at spatial frequency.
            
*Gaussian*: Filter is the Gaussian function: EXP(-F**2 / (2 * SPF**2)), 
where F is the frequency.
              
*Fermi*: Filter is: 1 / (1 + EXP((F - SPF) / T)) which negotiates 
between "Top-hat" and Gaussian characteristics, depending on 
the value of the temperature T.

*Butterworth* Filter is: 1 / (SQRT(1 + F / RAD)**(2 * ORDER)) 
The ORDER determines the filter fall off and RAD corresponds 
to the cut-off radius. 

*Raised cosine* Filter is: 0.5 * (COS(PI * (F - Flow) / (Flow - Fup)) + 1) 
if Flow < F < Fup, 1 if F < Flow, and 0 if F > Fup

See detailed description of the filter at [[https://spider.wadsworth.org/spider_doc/spider/docs/man/fq.html][SPIDER's FQ online manual]]
                           """)
        form.addParam('filterMode', EnumParam, choices=['low-pass', 'high-pass'],
                      label='Filter mode', default=FILTER_LOWPASS)
        form.addParam('usePadding', BooleanParam, default=True, 
                      label='Use padding?',
                      help='If set to *Yes*, to improve boundary quality\n'
                           'the image is padded with the average value to\n'
                           'twice the original size during filtration.\n\n'
                           'If *No* padding is applied, this may lead to\n'
                           'artifacts near boundary of image.')
        form.addParam('filterRadius', DigFreqParam, default=0.12, 
                      label='Filter radius (0 < f < 0.5)',
                      condition='filterType <= %d or filterType == %d' %
                                (FILTER_SPACE_REAL, FILTER_FERMI),
                      help='Low frequency cutoff to apply the filter.\n')  
        
        line = form.addLine('Frequency', 
                            condition='filterType > %d' % FILTER_FERMI,
                            help='Range to apply the filter. Expected values between 0 and 0.5.')
        line.addParam('lowFreq', DigFreqParam, default=0.1, label='Lowest')
        line.addParam('highFreq', DigFreqParam, default=0.2, label='Highest')
         
        form.addParam('temperature', FloatParam, default=0.3, 
                      label='Temperature T:',
                      condition='filterType == %d' % FILTER_FERMI,
                      help='Enter a temperature parameter T The filter falls off roughly within \n'
                           'this reciprocal distance (in terms of frequency units).')     
        
    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        # Define some names
        self.particlesStk = self._getPath('%(particles)s.%(ext)s' % self._params)
        # Insert processing steps
        self._insertFunctionStep('convertInput', 'inputParticles', 
                                 self._getFileName('particles'),
                                 self._getFileName('particlesSel'),
                                 needsGPU=False)
        self._insertFunctionStep('filterStep', self.filterType.get(),
                                 needsGPU=False)
        self._insertFunctionStep('createOutputStep', needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def filterStep(self, filterType):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        particles = self.inputParticles.get()
        n = particles.getSize()
        OP = self._op
        args = []

        if not self.usePadding:
            OP += ' NP'
        
        if filterType <= FILTER_FERMI:
            args.append(self.filterRadius.get())
        else:
            args.append('%f %f' % (self.lowFreq, self.highFreq))
            
        if filterType == FILTER_FERMI:
            args.append(self.temperature.get())
        
        # Map to expected filter number in Spider for operation FQ
        filterNumber = filterType * 2 + 1
        # Consider low-pass or high-pass
        filterNumber += self.filterMode.get()

        self._enterWorkingDir()  # Do operations inside the run working dir
        
        spi = SpiderShell(ext=self.getExt())  # Create the Spider process to send commands
        particlesStk = removeBaseExt(self.particlesStk)
        
        # Run a loop for filtering
        locStr = particlesStk + '@******[part]'
        cmds = ['do lb5 [part] = 1,%d' % n,
                OP, locStr, locStr, filterNumber] + args + ['lb5']
        
        for c in cmds:
            spi.runCmd(c)
            
        spi.close()
            
        self._leaveWorkingDir()  # Go back to project dir

    def createOutputStep(self):
        particles = self.inputParticles.get()
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(particles)

        updateItem = lambda p, i: p.setLocation(i, self.particlesStk)
        imgSet.copyItems(particles,
                         updateItemCallback=updateItem,
                         itemDataIterator=iter(range(1, particles.getSize()+1)))

        self._defineOutputs(**{outputs.outputParticles.name: imgSet})
        self._defineTransformRelation(particles, imgSet)
        
# --------------------------- INFO functions ----------------------------------
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        pixelSize = self.inputParticles.get().getSamplingRate()
        
        summary = list()
        summary.append('Used filter: *%s %s*' %
                       (self.getEnumText('filterType'),
                        self.getEnumText('filterMode')))
 
        if self.filterType <= FILTER_FERMI:
            summary.append('Filter radius: *%s px^-1*' % self.filterRadius)
            radiusAngstroms = pixelSize / self.filterRadius.get()
            summary.append('Filter radius: *%s Angstroms*' % radiusAngstroms)
        else:
            summary.append('Frequency range: *%s - %s*' % (self.lowFreq,
                                                           self.highFreq))

        if self.filterType == FILTER_FERMI:
            summary.append('Temperature factor: *%s*' % self.temperature)

        summary.append('Padding set to: *%s*' % self.usePadding)
        return summary
    
    def _methods(self):
        methods = []
        msg = '\nInput particles %s were %s filtered using a %s filter' %\
              (self.getObjectTag('inputParticles'),
               self.getEnumText('filterMode'),
               self.getEnumText('filterType'))

        if self.filterType <= FILTER_FERMI:
            msg += ', using a radius of %s px^-1' % self.filterRadius
        else:
            msg += ', using a frequency range of %s to %s px^-1' %\
                   (self.lowFreq, self.highFreq)

        if self.filterType == FILTER_FERMI: 
            msg += ' and a temperature factor of of %s px^-1' % self.temperature
            
        if self.usePadding:
            msg += ', padding the images by a factor of two.'
        else:
            msg += ' with no padding.'

        methods.append(msg)
        methods.append('Output particles: %s' %
                       self.getObjectTag('outputParticles'))

        return methods
