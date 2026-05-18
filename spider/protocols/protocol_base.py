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

from pwem.protocols import EMProtocol

from .. import Plugin
from ..utils import runTemplate
from ..convert import writeSetOfImages


class SpiderProtocol(EMProtocol):
    """
    Provides a general execution framework for SPIDER-based image processing
    workflows within cryo-EM pipelines. The protocol acts as a bridge between
    Scipion image datasets and SPIDER processing utilities, enabling image
    conversion, management of intermediate files, execution of SPIDER scripts,
    and validation of processing results in a consistent and reproducible way.

    AI Generated:

    SPIDER Base Protocol (SpiderProtocol) — User Manual
        Overview

        The SPIDER Base Protocol provides the infrastructure required to run
        SPIDER image-processing procedures inside a Scipion workflow. Rather
        than representing a specific biological analysis itself, this protocol
        serves as the operational foundation for many higher-level SPIDER
        classification, dimensionality reduction, and image analysis tools.

        In practical cryo-EM workflows, users frequently combine multiple
        software ecosystems. SPIDER remains historically important for
        multivariate statistical analysis, classification, and several
        particle-processing operations. This protocol allows those operations
        to be integrated into modern Scipion projects while maintaining
        compatibility with SPIDER conventions and file formats.

        General Workflow

        The protocol prepares cryo-EM image datasets so they can be processed
        by SPIDER procedures. Input particle sets are converted into SPIDER-
        compatible stacks together with associated selection files, ensuring
        that downstream SPIDER operations receive the data in the expected
        structure and organization.

        Once the data have been prepared, the protocol launches predefined
        SPIDER workflows using configurable execution templates. These
        templates define the scientific operation to be performed, while the
        protocol itself manages execution logistics, working directories,
        parameter passing, and output organization.

        Biological Role in Cryo-EM Pipelines

        Although this protocol does not perform biological interpretation on
        its own, it enables many biologically meaningful analyses that rely on
        SPIDER infrastructure. Examples include particle classification,
        multivariate statistical analysis, correspondence analysis, principal
        component analysis, and clustering workflows.

        In cryo-EM studies, these analyses are essential for separating
        heterogeneous particle populations, identifying conformational states,
        improving class averages, and reducing structural variability prior to
        three-dimensional reconstruction. The protocol therefore functions as
        a key interoperability layer between data management and advanced image
        analysis.

        File and Dataset Management

        A major responsibility of the protocol is maintaining a consistent
        organization of intermediate and final processing files. SPIDER
        workflows often generate large collections of stacks, documents, and
        auxiliary metadata. Proper organization is important both for
        reproducibility and for later inspection of results.

        The protocol automatically manages working directories and processing
        paths so that users can focus on the biological workflow rather than
        on manual file handling. This is especially valuable in large cryo-EM
        projects where multiple classification or dimensionality-reduction
        analyses may be executed simultaneously.

        Parallel Execution and Computational Considerations

        SPIDER workflows can vary substantially in computational demand.
        Depending on the downstream analysis, processing may involve large
        particle stacks and computationally intensive statistical operations.
        The protocol supports execution environments ranging from simple
        workstation runs to larger parallelized infrastructures.

        In routine practice, users should consider the size of the particle
        dataset and the complexity of the intended analysis before launching
        large SPIDER jobs. Proper allocation of computational resources can
        significantly reduce runtime and improve workflow stability.

        Error Detection and Reliability

        The protocol includes mechanisms to monitor execution status and detect
        failures during SPIDER processing. This is particularly important
        because numerical image-analysis workflows may fail due to invalid
        inputs, incompatible dimensions, insufficient disk space, or unstable
        parameter combinations.

        From a practical perspective, users should always inspect processing
        logs and intermediate outputs when unexpected results occur. Early
        detection of problematic datasets or unstable analyses can prevent the
        propagation of artifacts into downstream biological interpretation.

        Integration Within Scipion

        One of the major strengths of this protocol is its integration into
        the broader Scipion ecosystem. SPIDER-based analyses can be connected
        directly to particle extraction, preprocessing, classification,
        refinement, and visualization workflows without requiring manual file
        conversion outside the platform.

        This interoperability allows biological users to combine historical
        SPIDER methodologies with modern cryo-EM processing pipelines,
        facilitating both reproducibility and methodological flexibility.

        Practical Recommendations

        In routine biological workflows, users should ensure that particle
        stacks are properly curated before launching SPIDER analyses. Poorly
        centered particles, inconsistent preprocessing, or highly noisy data
        may reduce the quality of downstream statistical analyses and
        classifications.

        It is also advisable to maintain clear documentation of processing
        parameters and intermediate results, particularly in projects
        involving heterogeneous conformational states or exploratory
        classification strategies.

        Final Perspective

        The SPIDER Base Protocol is best understood as an enabling framework
        for advanced cryo-EM image analysis rather than a standalone analysis
        method. By managing data conversion, workflow execution, and
        integration with SPIDER utilities, it allows classical statistical
        and classification methodologies to remain fully accessible within
        modern Scipion environments. Its value lies in providing reliable and
        reproducible access to SPIDER capabilities while simplifying the
        operational complexity traditionally associated with these workflows.
    """
    _label = None
    _params = None
            
    def convertInput(self, attrName, stackFn, selFn):
        """ Convert from an input pointer of SetOfImages to Spider.
        Params:
            attrName: the attribute name of the input pointer
            stackFn: the name of the stack for converted images
            selFn: the name of the selection file.
        """
        imgSetPointer = getattr(self, attrName)
        writeSetOfImages(imgSetPointer.get(), stackFn, selFn)
        
    def _getFileName(self, key, *args):
        """ Give a key, append the extension
        and prefix the protocol working dir. 
        """
        template = '%(' + key + ')s.%(ext)s'
        
        return self._getPath(template % self._params)
    
    def getExt(self):
        """ Return the extension used in the script,
        stored in a dictionary called self._params. 
        """
        return self._params['ext']
    
    def getScript(self):
        return getattr(self, '_script', None)
    
    def runTemplate(self, inputScript, ext, paramsDict, nummpis=1):
        """ This function will create a valid Spider script
        by copying the template and replacing the values in dictionary.
        After the new file is read, the Spider interpreter is invoked.
        """
        self._enterWorkingDir()

        log = getattr(self, '_log', None)
        mpiFlag = True if nummpis > 1 else False
        program = Plugin.getProgram(mpiFlag)
        runTemplate(inputScript, ext, paramsDict, nummpis=nummpis,
                    program=program, log=log)
        self._leaveWorkingDir()
    
        f = open(self.getLogPaths()[0], 'r')
        for line in f.readlines():
            if 'FATAL ERROR ENCOUNTERED IN BATCH MODE' in line:
                raise RuntimeError('Spider script error!')
        f.close()
