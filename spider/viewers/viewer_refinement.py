# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se)
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
from glob import glob

import pyworkflow.protocol.params as params
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO
from pwem.viewers import (EmPlotter, ChimeraView,
                          EmProtocolViewer, ChimeraAngDist)

from ..constants import *
from ..protocols import SpiderProtRefinement
from ..utils import SpiderDocFile


class SpiderViewerRefinement(EmProtocolViewer):
    """ Visualization of Spider refinement results. """

    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [SpiderProtRefinement]
    _label = 'viewer refinement'

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('viewIter', params.EnumParam,
                      choices=['last', 'selection'], default=ITER_LAST,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label="Iteration to visualize", important=True,
                      help="""
*last*: only the last iteration will be visualized.
*selection*: you may specify a range of iterations.
Examples:
"1,5-8,10" -> [1,5,6,7,8,10]
"2,6,9-11" -> [2,6,9,10,11]
"2 5, 6-8" -> [2,5,6,7,8]                      
                           """)
        form.addParam('iterSelection', params.NumericRangeParam,
                      condition='viewIter==%d' % ITER_SELECTION,
                      label="Iterations list",
                      help="Write the iteration list to visualize.")

        group = form.addGroup('Angular assignment')
        group.addParam('displayAngDist', params.EnumParam,
                       choices=['2D plot', 'chimera'],
                       default=ANGDIST_2DPLOT,
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Display angular distribution',
                       help='*2D plot*: display angular distribution as interative 2D in matplotlib.\n'
                            '*chimera*: display angular distribution using Chimera with red spheres.')
        group.addParam('spheresScale', params.IntParam, default=-1,
                       condition="displayAngDist == %d" % ANGDIST_CHIMERA,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Spheres size',
                       help='')

        group = form.addGroup('Volumes')
        group.addParam('displayVol', params.EnumParam,
                       choices=['slices', 'chimera'],
                       default=VOLUME_SLICES, display=params.EnumParam.DISPLAY_HLIST,
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')
        group.addParam('showVolumes', params.EnumParam, default=VOL,
                       choices=self.volList(),
                       label='Volume to visualize',
                       help='Select the volume to visualize')

        group = form.addGroup('Resolution')

        group.addParam('showFSC', params.LabelParam, default=True,
                       important=True,
                       label='Display resolution plots (FSC)')
        group.addParam('groupFSC', params.EnumParam, default=0,
                       choices=['iterations', 'defocus groups'],
                       condition='not isGoldStdProt and viewIter==0',
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Group FSC plots by',
                       help='Select which FSC curve you want to '
                            'show together in the same plot.')
        group.addParam('groupSelection', params.NumericRangeParam,
                       condition='not isGoldStdProt and groupFSC==1 and viewIter==0',
                       label="Groups list",
                       help="Write the group list to visualize. See examples in iteration list")
        group.addParam('resolutionThresholdFSC', params.FloatParam, default=0.143,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Threshold in resolution plots',
                       help='Use 0.5 for refinement with defocus groups, otherwise '
                       'for gold-standard refinement use 0.143')

    def isGoldStdProt(self):
        # True = gold std refinement
        return self.protocol.protType == 1

    def volList(self):
        if self.isGoldStdProt():
            return ['reconstructed', 'half1', 'half2', 'filtered',
                    'filtered & centered']
        else:
            return ['reconstructed', 'half1', 'half2']

    def _getVisualizeDict(self):
        # self._load()
        return {
            'showVolumes': self._showVolumes,
            'showFSC': self._showFSC,
            'displayAngDist': self._displayAngDist
        }

    def _validate(self):
        return []

    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999.
        if value:
            inv = 1 / value
        return "1/%0.2f" % inv

    def _getIterations(self):
        if self.viewIter == ITER_LAST:
            return [self.protocol._getLastIterNumber()]
        else:
            return self._getListFromRangeString(self.iterSelection.get(''))

    def _getGroups(self):
        return self._getListFromRangeString(self.groupSelection.get(''))

    def _getFinalPath(self, *paths):
        return self.protocol._getExtraPath('Refinement', 'final', *paths)

# =========================================================================
# ShowVolumes
# =========================================================================
    def _createVolumesSqlite(self):
        """ Write an sqlite with all volumes selected for visualization. """

        volSqlite = self.protocol._getExtraPath('viewer_volumes.sqlite')
        samplingRate = self.protocol.inputParticles.get().getSamplingRate()
        self.createVolumesSqlite(self.getVolumeNames(),
                                 volSqlite, samplingRate)

        return [self.objectView(volSqlite)]

    def getVolumeNames(self, it=None):
        """ If it is not none, return the volume of this iteration only. """
        if it is None:
            iterations = self._getIterations()
        else:
            iterations = [it]

        if self.isGoldStdProt():
            volTemplate = VOLNAMES_GOLDSTD[self.showVolumes.get()]
        else:
            volTemplate = VOLNAMES_DEFGROUPS[self.showVolumes.get()]

        volumes = [self._getFinalPath(volTemplate % i) + '.stk'
                   for i in iterations]

        return volumes

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self.getVolumeNames()
        cmdFile = self._getFinalPath('chimera_volumes.cxc')
        with open(cmdFile, 'w+') as f:
            for vol in volumes:
                localVol = os.path.relpath(vol,
                                           self.protocol._getFinalPath())
                if os.path.exists(vol):
                    f.write("open %s format spider\n" % localVol)
            f.write('tile\n')
        view = ChimeraView(cmdFile)
        return [view]

    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()

        elif self.displayVol == VOLUME_SLICES:
            return self._createVolumesSqlite()

# ===============================================================================
# plotFSC
# ===============================================================================

    def _plotFSC(self, a, fscFile):
        resolution = []
        fsc = []

        fscDoc = SpiderDocFile(fscFile)
        for values in fscDoc:
            resolution.append(1 / values[1])
            fsc.append(values[2])

        self.maxfsc = max(fsc)
        self.minInv = min(resolution)
        self.maxInv = max(resolution)
        a.plot(resolution, fsc)
        from matplotlib.ticker import FuncFormatter
        a.xaxis.set_major_formatter(FuncFormatter(self._formatFreq))
        a.set_ylim([-0.1, 1.1])
        fscDoc.close()

    def _showFSC(self, paramName=None):
        threshold = self.resolutionThresholdFSC.get()
        iterations = self._getIterations()
        groups = self._getGroups()

        if self.isGoldStdProt():
            template = 'fscdoc_m_%02d.stk'
            title = 'Masked FSC'
        else:
            template = 'fscdoc_%02d.stk'
            title = 'FSC'

        if self.groupFSC == 0:  # group by iterations
            files = [(it, self._getFinalPath(template % it)) for it in iterations]
            legendPrefix = 'iter'
        else:
            it = iterations[-1]  # show only last iteration
            legendPrefix = 'group'

            def group(f):  # retrieve the group number
                return int(f.split('_')[-1].split('.')[0])

            groupFiles = glob(self._getFinalPath('ofscdoc_%02d_???.stk' % it))
            groupFiles.sort()
            files = [(group(f), f) for f in groupFiles if group(f) in groups]
            if not files:  # empty files
                return [self.errorMessage("Please select valid groups to display",
                                          title="Wrong groups selection")]

        plotter = EmPlotter(windowTitle='Resolution FSC')
        a = plotter.createSubPlot(title, 'Angstroms^-1', 'FSC')
        legends = []
        for it, fscFile in files:
            if os.path.exists(fscFile):
                self._plotFSC(a, fscFile)
                legends.append('%s %d' % (legendPrefix, it))
            else:
                print("Missing file: ", fscFile)

        # plot final FSC curve (from BP)
        if self.groupFSC == 0 and not self.isGoldStdProt():
            lastIter = self.protocol._getLastIterNumber()
            if lastIter in iterations:
                fscFinalFile = self._getFinalPath('ofscdoc_%02d.stk' % lastIter)
                if os.path.exists(fscFinalFile):
                    self._plotFSC(a, fscFinalFile)
                    legends.append('final')

        if threshold < self.maxfsc:
            a.plot([self.minInv, self.maxInv], [threshold, threshold],
                   color='black', linestyle='--')

        plotter.showLegend(legends)
        a.grid(True)

        return [plotter]

    def _iterAngles(self, it):
        """ Iterate over the angular distribution for a given iteration. """
        # Get the alignment files of each group for this iteration
        if self.isGoldStdProt():
            template = 'align_%02d_???_s?.stk'
        else:
            template = 'align_%02d_???.stk'

        files = glob(self._getFinalPath(template % it))
        for anglesFile in files:
            alignDoc = SpiderDocFile(anglesFile)
            for values in alignDoc:
                theta = values[1]
                phi = values[2]

                if theta > 90:
                    theta = abs(180. - theta)
                    phi += 180
                yield phi, theta
            alignDoc.close()

    def _displayAngDist(self, *args):
        iterations = self._getIterations()
        nparts = self.protocol.inputParticles.get().getSize()
        views = []

        if self.displayAngDist == ANGDIST_2DPLOT:
            for it in iterations:
                if it == 1:
                    print("Orientations for the first iteration cannot be plotted. "
                          "Skipping..")
                    continue
                anglesSqlite = self._getFinalPath('angular_dist_%03d.sqlite' % it)
                title = 'Angular distribution iter %03d' % it
                plotter = EmPlotter(windowTitle=title)
                self.createAngDistributionSqlite(anglesSqlite, nparts,
                                                 itemDataIterator=self._iterAngles(it))
                plotter.plotAngularDistributionFromMd(anglesSqlite, title)
                views.append(plotter)
        else:
            it = iterations[-1]
            print("Using last iteration: ", it)
            anglesSqlite = self._getFinalPath('angular_dist_%03d.sqlite' % it)
            self.createAngDistributionSqlite(anglesSqlite, nparts,
                                             itemDataIterator=self._iterAngles(it))
            volumes = self.getVolumeNames(it)
            vol = self.protocol.outputVolume
            volOrigin = vol.getOrigin(force=True).getShifts()
            samplingRate = vol.getSamplingRate()
            
            views.append(ChimeraAngDist(volumes[0],
                                        self.protocol._getTmpPath(),
                                        voxelSize=samplingRate,
                                        volOrigin=volOrigin,
                                        angularDistFile=anglesSqlite,
                                        format="spider"))

        return views
