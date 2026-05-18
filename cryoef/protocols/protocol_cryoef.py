# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pwem.protocols import ProtAnalysis3D
from pwem.objects import Volume

from cryoef import Plugin
from ..convert import writeAnglesFn, parseOutput


class ProtCryoEF(ProtAnalysis3D):
    """
    Analyses the orientation distribution of single-particle cryo-EM data
    to evaluate sampling quality, directional coverage, and expected
    reconstruction performance in both real and Fourier space.

    AI Generated:

    Orientation Analysis (ProtCryoEF) - User Manual
        Overview

        The Orientation Analysis protocol evaluates how particle
        orientations are distributed within a single-particle cryo-EM
        dataset. Its main purpose is to determine whether the angular
        sampling of the experimental images is sufficiently uniform to
        support a reliable three-dimensional reconstruction. Uneven
        orientation distributions are a common limitation in cryo-EM and
        can lead to directional resolution anisotropy, missing structural
        information, or reconstruction artifacts.

        In practical biological workflows, this analysis is especially
        useful after refinement or reconstruction steps where particle
        orientations have already been assigned. By examining the
        orientation coverage, users can better understand the strengths
        and weaknesses of their dataset before investing additional effort
        into high-resolution interpretation, atomic modeling, or
        classification strategies.

        Inputs and Experimental Context

        The protocol requires a set of particles containing projection
        alignment information. These angular assignments are interpreted
        as the experimental sampling of projection directions across the
        molecule. Because the analysis depends entirely on the quality of
        these orientations, the protocol is most meaningful when applied
        to particles obtained from a well-converged refinement.

        The symmetry definition is biologically important because it
        determines how orientations are interpreted during the analysis.
        Symmetric particles may compensate for incomplete angular
        sampling through redundancy introduced by symmetry operations,
        whereas asymmetric particles require much broader experimental
        coverage. Selecting the correct symmetry group is therefore
        essential for obtaining biologically meaningful conclusions.

        The approximate particle diameter is also required because the
        protocol estimates the effective sampling behavior relative to the
        physical size of the reconstructed object. Accurate diameter
        values improve the realism of the resulting point spread function
        estimations and resolution predictions.

        Angular Accuracy and Resolution Estimation

        The protocol allows users to specify the expected angular
        accuracy of the dataset. This parameter reflects the uncertainty
        associated with particle orientation assignments during
        refinement. Smaller angular uncertainties generally indicate
        better refinement quality and more reliable directional sampling,
        while larger uncertainties may reduce the effective resolution
        achievable in the reconstruction.

        An optional B-factor can also be provided to estimate the decay
        of high-resolution signal within the dataset. In biological cryo-EM
        studies, the B-factor often reflects cumulative effects from beam
        damage, conformational heterogeneity, particle flexibility, and
        alignment inaccuracies. Providing a realistic estimate improves
        the interpretation of the predicted resolution distribution.

        Users may additionally define an FSC resolution value when an
        independent resolution estimate is already available from a
        previous reconstruction workflow. This allows the protocol to
        incorporate experimentally derived resolution information into the
        orientation analysis. When no value is supplied, the protocol
        estimates the effective resolution behavior automatically.

        Tilt Constraints and Preferred Orientation

        Preferred orientation is one of the most common challenges in
        cryo-EM specimen preparation. Many biological complexes adopt
        limited orientations on the grid, producing anisotropic sampling
        in Fourier space. The protocol includes a maximum tilt parameter
        that helps model the expected angular accessibility during data
        collection and reconstruction.

        From a biological perspective, datasets suffering from strong
        preferred orientation may exhibit good resolution in some
        directions but poor resolution in others. This effect can limit
        the interpretability of flexible regions, membrane domains, or
        elongated assemblies. Evaluating orientation coverage early in
        the workflow can guide decisions such as collecting tilted data,
        modifying grid preparation conditions, or optimizing biochemical
        stabilization strategies.

        Outputs and Their Interpretation

        The protocol generates representations of the point spread
        function in both real space and Fourier space. These outputs
        provide complementary views of how angular sampling influences
        the expected reconstruction quality.

        The real-space representation helps visualize how directional
        information propagates through the reconstructed volume. The
        Fourier-space representation reveals regions of strong or weak
        sampling and is especially useful for identifying anisotropy or
        missing angular coverage.

        Quantitative metrics are also produced to summarize the
        efficiency and uniformity of the orientation distribution. These
        include estimates of mean resolution behavior, variability across
        directions, and the best and worst expected reconstruction
        conditions. Together, these measurements provide a practical
        assessment of whether the dataset supports isotropic
        high-resolution interpretation.

        Practical Recommendations

        In routine cryo-EM workflows, orientation analysis is most useful
        immediately after obtaining a stable three-dimensional refinement.
        Datasets with broad angular coverage and low directional bias are
        generally more suitable for downstream atomic interpretation and
        flexible-region analysis.

        When strong anisotropy is detected, biological users should
        consider strategies to improve orientation diversity. These may
        include collecting tilted datasets, changing support films,
        adjusting buffer composition, modifying detergent conditions for
        membrane proteins, or exploring alternative sample preparation
        protocols.

        Symmetry should always be verified carefully before analysis.
        Incorrect symmetry assumptions can artificially improve or worsen
        the apparent orientation distribution and may lead to misleading
        conclusions regarding reconstruction quality.

        Final Perspective

        Orientation distribution analysis is not simply a geometric
        diagnostic but an essential biological quality-control step in
        cryo-EM reconstruction. Understanding how particles sample
        different viewing directions provides critical insight into the
        reliability, isotropy, and interpretability of the final map.
        Careful interpretation of angular coverage can help guide both
        experimental optimization and downstream structural analysis.
    """
    _label = 'orientation analysis'
    _devStatus = PROD
    _possibleOutputs = {
        'outputVolume1': Volume,
        'outputVolume2': Volume
    }

    def __init__(self, **kwargs):
        ProtAnalysis3D.__init__(self, **kwargs)

    def _initialize(self):
        """ This function is mean to be called after the
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
                  'anglesFn': self._getExtraPath('input_angles.dat'),
                  'projections': self._getExtraPath('input_projections.sqlite'),
                  'output_log': self._getExtraPath('input_angles.log'),
                  'real space PSF': self._getExtraPath('input_angles_R.mrc'),
                  'fourier space PSF': self._getExtraPath('input_angles_K.mrc'),
                  'output_hist': self._getExtraPath('input_angles_PSFres.dat')
                  }

        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles", important=True,
                      help='Provide input particles with angular information.')
        form.addParam('symmetryGroup', params.StringParam, default='c1',
                      label="Symmetry",
                      help='If the molecule is asymmetric, set Symmetry group '
                           'to C1. Look at the XMIPP Wiki for more details:'
                           ' https://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/'
                           'WebHome?topic=Symmetry')
        form.addParam('diam', params.IntParam, default=200,
                      label='Particle diameter (A)',
                      help='Approximate particle diameter, in Angstroms.')
        form.addParam('angAcc', params.IntParam, default=1,
                      label='Angular accuracy (deg)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Angular accuracy in degrees.')
        form.addParam('Bfact', params.IntParam, default=160,
                      label='B-factor (A^2)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='B-factor estimate for your data, if one was '
                           'estimated for the 3D reconstruction.')
        form.addParam('FSCres', params.FloatParam, default=-1,
                      label='FSC resolution (A)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='FSC resolution using 0.143 criterion. '
                           'Default (-1) value means that resolution will be '
                           'automatically estimated from B-factor.')
        form.addParam('maxTilt', params.IntParam, default=45,
                      label='Max tilt angle (deg)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Maximum tilt angle allowed for prediction '
                           'algorithm, in degrees.')

    # --------------------------- INSERT steps functions ----------------------
    
    def _insertAllSteps(self):
        # Insert processing steps
        self._initialize()
        self._insertFunctionStep('convertInputStep', needsGPU=False)
        self._insertFunctionStep('runCryoEFStep', needsGPU=False)
        self._insertFunctionStep('createOutputStep', needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    
    def convertInputStep(self):
        """ Convert input angles as expected by cryoEF."""
        partSet = self._getInputParticles()
        anglesFn = self._getFileName('anglesFn')
        with open(anglesFn, 'a') as f:
            for part in partSet:
                writeAnglesFn(part, f)

    def runCryoEFStep(self):
        """ Call cryoEF with the appropriate parameters. """
        args = self._getArgs()
        param = ' '.join(['%s %s' % (k, str(v)) for k, v in args.items()])
        program = Plugin.getProgram()

        self.runJob(program, param, env=Plugin.getEnviron())

    def createOutputStep(self):
        partSet = self._getInputParticles()

        vol = Volume()
        vol.setSamplingRate(partSet.getSamplingRate())
        vol.setObjLabel('real space PSF')
        vol.setFileName(self._getFileName('real space PSF'))

        vol2 = Volume()
        vol2.setSamplingRate(partSet.getSamplingRate())
        vol2.setObjLabel('fourier space PSF')
        vol2.setFileName(self._getFileName('fourier space PSF'))

        outputs = {'outputVolume1': vol,
                   'outputVolume2': vol2}
        self._defineOutputs(**outputs)
        self._defineSourceRelation(self.inputParticles, vol)
        self._defineSourceRelation(self.inputParticles, vol2)

    # --------------------------- INFO functions ------------------------------
    
    def _summary(self):
        summary = []

        if hasattr(self, 'outputVolume1'):
            results = list(parseOutput(self._getExtraPath('input_angles.log')))
            eff, meanRes, stdev, worstRes, bestRes = results
            summary.append('Efficiency of the orientation distribution: *%0.2f*' % eff)
            summary.append('Mean PSF resolution: *%0.2f A*' % meanRes)
            summary.append('Standard deviation: *%0.2f A*' % stdev)
            summary.append('Worst PSF resolution: *%0.2f A*' % worstRes)
            summary.append('Best PSF resolution: *%0.2f A*' % bestRes)
        else:
            summary.append("Output is not ready yet.")

        return summary
    
    def _validate(self):
        errors = []

        return errors
    
    # --------------------------- UTILS functions -----------------------------
 
    def _getArgs(self):
        """ Prepare the args dictionary."""
        args = {'-f': self._getFileName('anglesFn'),
                '-b': self._getInputParticles().getFirstItem().getXDim(),
                '-a': self.angAcc.get(),
                '-B': self.Bfact.get(),
                '-D': self.diam.get(),
                '-g': self.symmetryGroup.get() or 'c1',
                '-m': self.maxTilt.get()
                }
        if self.FSCres.get() != -1:
            args['-r'] = self.FSCres.get()

        return args

    def _getInputParticles(self):
        return self.inputParticles.get()
