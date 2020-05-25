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
"""
Protocol wrapper around the cryoEF tool for analysing the orientation
distribution of single-particle EM data
"""

import pyworkflow.protocol.params as params
from pwem.protocols import ProtAnalysis3D
from pwem.objects import Volume
from cryoef import Plugin
from cryoef.convert import writeAnglesFn, parseOutput


class ProtCryoEF(ProtAnalysis3D):
    """
    Protocol for analysing the orientation distribution of single-particle EM data.
     
    Find more information at http://www.mrc-lmb.cam.ac.uk/crusso/cryoEF/
    """
    _label = 'orientation analysis'

    def __init__(self, **kwargs):
        ProtAnalysis3D.__init__(self, **kwargs)

    def _initialize(self):
        """ This function is mean to be called after the
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
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
                           ' http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/'
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
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runCryoEFStep')
        self._insertFunctionStep('createOutputStep')

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
