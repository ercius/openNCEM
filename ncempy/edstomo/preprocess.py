import sys, os, shutil
import re
from collections import OrderedDict
import json
import numpy as np
from skimage.external import tifffile
from scipy.ndimage.interpolation import shift
from ncempy.edstomo.CharacteristicEmission import GetFluorescenceLineEnergy
from ncempy.io.emd import fileEMD

def BinEDSSpatialDimensions(Arr, Binning=4):
    if len(Arr.shape) != 3:
        print('Input stack needs to have a shape like (tilt, x, y).')
        return Arr

    if ((Arr.shape[1] % Binning) != 0) or ((Arr.shape[2] % Binning) != 0):
        print('Input stack needs to have dimensions that are evenly divided by the binning factor.')
        return Arr

    # Create new axes with a length of Binning and fold the array into them.
    ArrNew = Arr.reshape((  Arr.shape[0],                                   # Tilt dimension untouched
                            Binning,int(Arr.shape[1]/Binning),              # Spatial dimension binned
                            Binning,int(Arr.shape[2]/Binning)), order='F')   # Spatial dimension binned
    # Sum along the binning directions.
    ArrNew = np.sum(np.sum(ArrNew, axis=3), axis=1)
    # print('Old shape: ' + str(Arr.shape))
    # print('New shape: ' + str(ArrNew.shape))
    return ArrNew

def ExtractSignalsFromEMD(InputEMD=None, SignalNames=['HAADF', 'Mg_K', 'Fe_Ka'], Binning=4):
    ''' Read in a set of Bruker bcf files containing EDS acquisitions.

    Parameters:
        InputEMD (str): Name of EMD file containing the tilt stacks.

        SignalNames (list of str): A list of signals that should be extracted from the Bruker files.
            Valid values include:
            HAADF: For the HAADF signal
            Element_Line: e.g. Fe_Ka for a weighted sum of Fe_Ka1 and Ka2, or Fe_Ka1 for just that line, or Mg_K for a weighted sum of all Mg-K lines.  There is no brehmsstrahlung removal.
            eV1-eV2: Start and stop energies in eV.  The signal will be the sum of all energies over the range [eV1, eV2).

        Binning (int): How much binning to do on EDS signals to reduce noise.  1 means no binning.  2 means each output voxel is 2x2x2 input voxels.  HAADF signals are not rebinned as they are usually not noisy.

    Returns:
        SignalDict (OrderedDict of np.ndarray): Dictionary with names matching SignalNames and dimensions of (tilt, x, y).
        Tilts (list of floats): A list of tilt angles.

    '''

    # for s in SignalNames:
    #     if '_' in s:
    #         El, Line = s.split('_')
    #         print(GetFluorescenceLineEnergy(El, Series=Line[0], Line=Line))
    
    # Make an ordered dictionary with one entry for each signal.
    SignalDict = OrderedDict()
    for n in SignalNames:
        SignalDict[n] = []

    # Open the emd file that has our data.
    EMD = fileEMD(InputEMD, readonly=True)

    # We will compute the FWHM of peaks using the Mn-Ka reported in the EMD.
    EnergyResolutionMnKa = EMD.microscope.attrs['MnKaResolution[eV]']
    Mn_Ka_Energy = GetFluorescenceLineEnergy('Mn', Series='K', Line='Ka')
    K = EnergyResolutionMnKa / np.sqrt(Mn_Ka_Energy)
    print('Energy resolution of Mn-Ka is: ' + str(EnergyResolutionMnKa) + ' eV.')
    print('Assumed FWHM of peaks will be (%g*sqrt(E))/2., hence at Mn-Ka: %g eV.' % (K, K*np.sqrt(Mn_Ka_Energy)/2))

    # Get links to the data we'll need from the EMD.
    HAADF, HAADF_dims = EMD.get_emdgroup(EMD.data['HAADF_TiltStack'])
    EDS, EDS_dims  = EMD.get_emdgroup(EMD.data['EDS_TiltStack'])
    Tilts = HAADF_dims[0][0]
    print('HAADF dimensions are (%d, %d).'%(len(HAADF_dims[1][0]), len(HAADF_dims[1][0])))

    # Calculate the rebinning size for the EDS data.
    rebinsize_m = int(len(EDS_dims[1][0])/Binning)
    rebinsize_n = int(len(EDS_dims[2][0])/Binning)
    print('Binning is %d so rebinned EDS cubes will have spatial dimension (%d, %d).' % (Binning, rebinsize_m, rebinsize_n))

    for sig in SignalDict.keys():
        print(sig, end='')
        if sig == 'HAADF':
            SignalDict['HAADF'] = HAADF[:].astype("float32") # The HAADF doesn't get rebinned.
            print('')

        if '_' in sig:
            # This is a fluorescence line.  
            El, Line = sig.split('_')
            CenterEnergy = GetFluorescenceLineEnergy(El, Series=Line[0], Line=Line)
            if CenterEnergy is None:
                print('Unrecognized fluorescence line: %s-%s, ignoring this signal.'%(El,Line))
                continue
            LowEnergy = CenterEnergy - K*np.sqrt(CenterEnergy)/2
            HighEnergy = CenterEnergy + K*np.sqrt(CenterEnergy)/2
            LowEnergyIndex = np.argmin((EDS_dims[3][0] - LowEnergy)**2)
            HighEnergyIndex = np.argmin((EDS_dims[3][0] - HighEnergy)**2)
            print(', %g-%g eV window, energy bins: %d-%d.'% (LowEnergy, HighEnergy, LowEnergyIndex, HighEnergyIndex))
            Cube = np.sum(EDS[:, :, :, LowEnergyIndex:HighEnergyIndex+1], axis=-1)
            SignalDict[sig] = BinEDSSpatialDimensions(Cube, Binning)
            # fluor[0].rebin((rebinsize_m,rebinsize_n)).data.copy().astype("float32"))
        
        if '-' in sig:
            print('Energy range signals not implemented yet.')
    
    print('Signals Extracted.')

    # Turn those signal readouts into 3D numpy arrays (tilt, x, y).
    for k, v in SignalDict.items():
        npStack = np.array(v)
        SignalDict[k] = npStack

    return SignalDict, Tilts
    
def NormalizeSignals(SignalDict, Tilts, NormalizationSignalName=None, NormalizationImageFraction=1):
    ''' 
    Parameters:
        SignalDict (OrderedDict of np.ndarray): Dictionary with signal stacks.  Each entry is one tilt sequence with dimension (tilt, x, y).  The key is the name of the signal, e.g. 'Fe_Ka' or 'HAADF'.

        Tilts (list of floats): A list of tilt angles.

        NormalizationSignalName (str): Signal name to use for normalizing signal.  Preferably a higher energy but still high signal line.  Note, using HAADF for normalization produced bad results since the HAADF signal is not shadowed by the tomography holder.  None means don't normalize, 'Theory' means to normalize using a theoretical curve.

        NormalizationImageFraction (float): Fraction of frame (from center) to use to produce a normalization curve.  1 means to use the entire frame.  0.5 means to use the center half of the frame (0.25 - 0.75 of the frame width and height).

    Returns:
        SignalDict (OrderedDict of np.ndarray): Same as passed in but normalized.
        NormCurve (np.ndarray): Values to multiply by each tilt image in order to normalize the signal.
    '''

    # Normalize the spectra.
    NormCurve = np.ones(len(Tilts)) # Default to no normalization -- we're just multiplying by identity.
    if NormalizationSignalName is not None:
        if NormalizationSignalName == 'Theory':
            print('Theoretical normalization curve not implemented yet.')
        elif NormalizationSignalName == 'Independent':
            print('Producing normalization curve for each signal independently.')
            for k, v in SignalDict.items():
                NormCurve = GetNormalizationCurve(SignalDict[k], Tilts, NormalizationImageFraction)
                v = v * NormCurve[:,None,None]
                SignalDict[k] = v
        else:
            print('Producing normalization curve for all signals based on '+NormalizationSignalName+'.')
            NormCurve = GetNormalizationCurve(SignalDict[NormalizationSignalName], Tilts, NormalizationImageFraction)
            for k, v in SignalDict.items():
                v = v * NormCurve[:,None,None]
                SignalDict[k] = v

    return SignalDict, NormCurve

def GetNormalizationCurve(Signal, Tilts, ImageFraction=1):
    ''' Given a signal, finds the normalization curve to correct for the fact that sample thickness changes,
    and EDS detector shadowing changes as a function of tilt angle.

    Parameters:
        Signal (np.ndarray): Shape of (Tilts, x, y).  Signal from which to derive a normalization curve.

        Tilts: (np.ndarray): Stage tilts for the tomographic acquisition.

        ImageFraction (float): Fraction of the image to use to produce a norm curve.  Often higher tilt images see additional material coming into view which negates the values of the normalization curve.  1 means use the whole frame and would be valid for an isolated object entirely in the FOV.  0.5 means use only the range in x and y from 0.25-0.75 of the width and height of the frame.  Must be in the range 0 < ImageFraction <= 1.

    Returns:
        NormCurve (np.ndarray): Values to multiply by each tilt image in order to normalize the signal.
    '''

    # Convert to float
    Signal = Signal.astype(float)

    # Figure out the portion out of the center of the image to use to determine the curve.
    xmin = int(np.floor(Signal.shape[1] * ((1/2)-ImageFraction/2)))
    xmax = int(np.ceil(Signal.shape[1] * ((1/2)+ImageFraction/2)))
    ymin = int(np.floor(Signal.shape[2] * ((1/2)-ImageFraction/2)))
    ymax = int(np.ceil(Signal.shape[2] * ((1/2)+ImageFraction/2)))
    # Crop the signal to only include the center portion.
    Signal = Signal[:, xmin:xmax, ymin:ymax]

    # Get the index of the stack with the least tilt.
    MinTilt = np.argmin(np.array(Tilts)**2)
    # Now sum up the minimum tilt frame (usually zero tilt.)
    CenterIntensity = np.sum(Signal[MinTilt,:,:])
 
    # The normalization curve will have one value for each tilt.
    NormCurve = np.zeros(len(Tilts))

    # Loop through each frame and integrate it to get intensity.
    for i in range(Signal.shape[0]):
        NormCurve[i] = np.sum(Signal[i,:,:])
 
    # Set this so the center intensity is 1.  Therefore, after normalization, the least tilt map is unchanged.
    NormCurve /= CenterIntensity

    # Right now higher intensities mean we counted more counts due to a decrease in detector shadowing.
    # We want to multiply each stack by the value of NormCurve, so we need inverse.
    NormCurve = 1/NormCurve

    return NormCurve

def ReadImageJTranslations(FileName, Tilts):
    ''' The ImageJ multistackreg plugin is one way to align stacks.  This method reads in the translations file written by ImageJ.

    Parameters:
        FileName (str): Name of the file that was created by ImageJ.
        Tilts: (np.ndarray): Stage tilts for the tomographic acquisition.

    Returns:
        np.ndarray with shape (Tilts, 2) containing an x and y translation for each tilt.
    '''

    with open(FileName, 'r') as f:
        s = f.read()
        
    translations = re.findall(pattern='([\d\.]*)\t([\d\.]*)', string=s)

    # There are 6 lines in each translation from the output file.  
    # The zeroth is the new position, the third is the original position.
    # The other lines are for rotations which we are not using.
    # Example:
    # [('236.81374221505683', '288.85170101351565'),
    #  ('0.0', '0.0'),
    #  ('0.0', '0.0'),
    #  ('256.0', '256.0'),
    #  ('0.0', '0.0'),
    #  ('0.0', '0.0'),
    #  ('246.23227371186195', '229.14170151652533'),
    #  ('0.0', '0.0'),
    #  ('0.0', '0.0'),
    #  ('256.0', '256.0'),
    #  ('0.0', '0.0'),
    #  ('0.0', '0.0'),
    #  ('264.2877416986045', '281.32142526192746'),
    #  ('0.0', '0.0'),
    #  ('0.0', '0.0'),
    #  ('256.0', '256.0'),
    #  ('0.0', '0.0'),
    #  ('0.0', '0.0'),
    #  ...
     

    # There should be as many translations as there are images-1 because we are registering to the first image.
    assert(len(translations)/6 == len(Tilts)-1)

    offsets = []
    offsets.append(np.array((0,0))) # MultiStackReg will ignore the first image since others are registered relative to it.
    # Loop through all the translation info, image by image (6 lines at a time).
    for i in range(0, len(translations), 6):
        NewPos = translations[i] # Image shifted to here
        OrigPos = translations[i+3] # Image started here
        # How much to shift:
        offset = -np.array((float(NewPos[1]) - float(OrigPos[1]), float(NewPos[0]) - float(OrigPos[0])))
        if len(offsets) > 0:
            offset += offsets[-1]
        offsets.append(offset)

    offsets = np.array(offsets)
    print(offsets)
    return offsets

def ReadTomVizTranslations(FileName, Tilts):
    ''' TomViz has a capability for aligning stacks.  The manual method works very well because EDS tomographs usually only have a small number of frames.  This function reads in alignements generated by the TomViz manual method and saved to a json file.

    Parameters:
        FileName (str): Name of the file that was created by TomViz.
        Tilts: (np.ndarray): Stage tilts for the tomographic acquisition.

    Returns:
        np.ndarray with shape (Tilts, 2) containing an x and y translation for each tilt.
    '''

    with open(FileName, 'r') as f:
        s = json.load(f)
    s = np.array(s)
    # The coordinates in tomviz are different.  X and Y are swapped, and the sign of one axis is different.
    s = np.roll(s, 1, axis=1)
    s[:,0] = -s[:,0]
    return s

def ApplyTranslations(SignalDict, Translations, AlignSignalName):
    ''' Apply translations to all the images in all the signals.

    Parameters:
        SignalDict (OrderedDict of np.ndarray): Dictionary with signal stacks.  Each entry is one tilt sequence with dimension (tilt, x, y).  The key is the name of the signal, e.g. 'Fe_Ka' or 'HAADF'.

        Translations (np.ndarray with shape (Tilts, 2)): An x and y translation for each tilt.
        
        AlignSignalName (str): Name of the signal used to align the stack.  e.g. 'Fe_Ka' or 'HAADF'

    Returns:
        AlignedSignalDict: SignalDict after applying translations to each image in each signal.

    '''

    AlignedSignalDict = OrderedDict()

    # Alignment is done with one signal.  Usually this is HAADF but sometimes it is an EDS channel.
    # Since the binning is not the same, we have to scale the alignment to the other channels.
    AlignSignalShape = SignalDict[AlignSignalName].shape

    for SignalName, Sig in SignalDict.items():
        print('Aligning', SignalName)
        SigAlign = []

        # Calculate how much to translate this image based on it's shape relative to the signal used for alignment.
        TranslationScale = np.ones(2)
        TranslationScale[0] = Sig.shape[1]/AlignSignalShape[1]
        TranslationScale[1] = Sig.shape[2]/AlignSignalShape[2]
        print('TranslationScale: '+str(TranslationScale))

        # Apply the shifts to every image in this signal stack.
        for i in range(len(Translations)):
            ThisShift = np.round(Translations[i]*TranslationScale).astype(int)
            SigAlign.append(shift(Sig[i], ThisShift))

        # Turn the individual shifted images back into a 3D numpy array with dimension (tilts,x,y).
        SigAlign = np.array(SigAlign)

        AlignedSignalDict[SignalName] = SigAlign
    return AlignedSignalDict

def WriteSignalsToTIFFs(OutputDirectory, SignalDict):
    ''' Write the signals we have to tiff stacks so that we can do alignment, GENFIRE, etc. on them.  Also write auxiliary files that are needed by these processes.

    Parameters:
        OutputDirectory (str): Name of directory to place the tif files.

        SignalDict (OrderedDict of np.ndarray): Dictionary with signal stacks.  Each entry is one tilt sequence with dimension (tilt, x, y).  The key is the name of the signal, e.g. 'Fe_Ka' or 'HAADF'.

    Returns:
        Files written.
    '''

    # Ensure the output directory exists and is ready for us to write our output.
    if not os.path.exists(OutputDirectory):
        print('Creating directory: '+OutputDirectory)
        os.makedirs(OutputDirectory)

    # Write out all the stacks as tiff files.
    for k, v in SignalDict.items():
        # Reformat the stacks to have the correct shape and bit depth to write to tifs.
        npStack = v.astype(float)

        # Things need to fit in 16 bits.
        if np.max(npStack) > 2**16-1:
            npStack /= np.max(npStack)
            npStack *= 2**16-1
        
        # Now convert to uint16 for 16 bit unsigned tifs.
        npStack = npStack.round().astype('uint16')

        # Finally we get to save the actual tif.
        print('Writing '+str(k)+'.tif')
        tifffile.imsave(os.path.join(OutputDirectory,k+'.tif'), npStack)

def ReadSignalsFromTIFFs(InputDirectory, SignalNames=['HAADF', 'Mg_Ka', 'Fe_Ka']):
    ''' Read from TIFF files each of the requested signals.
    Parameters:
        InputDirectory (str): Name of directory containing the bcf files.

        SignalNames (list of str): A list of signals that should be extracted from the Bruker files.
            Valid values include:
            HAADF: For the HAADF signal
            Element_edge: e.g. Fe_Ka for a specific line.  There is no brehmsstrahlung removal.
            eV1-eV2: Start and stop energies in eV.  The signal will be the sum of all energies over the range [eV1, eV2).

    Returns:
        SignalDict (OrderedDict of np.ndarray): Dictionary with signal stacks.  Each entry is one tilt sequence with dimension (tilt, x, y).  The key is the name of the signal, e.g. 'Fe_Ka' or 'HAADF'.


    '''
    SignalDict = OrderedDict()
    for SigName in SignalNames:
        npStack = tifffile.imread(os.path.join(InputDirectory, SigName+'.tif')).astype(float)
        SignalDict[SigName] = npStack
    return SignalDict

def WriteMetaDataFiles(OutputDirectory, Tilts, NormCurve=None, NormalizationSignalName=None):
    ''' Writes out meta data obtained while processing the inputs which may be useful.
    Parameters:
        OutputDirectory (str): Name of directory to place the tif files.

        Tilts (list of floats): A list of tilt angles.

        NormCurve (np.ndarray): Values to multiply by each tilt image in order to normalize the signal.

    Returns:
        Files written.
    '''
    # Make the tilts file
    with open(os.path.join(OutputDirectory,'tilts.txt'), 'w') as f:
        print('Writing tilts file.')
        for n in Tilts:
            f.write('0.0 %0.1f 0.0\n'%(n))
    # Write out the normalization curve.
    if NormCurve is not None:
        np.savetxt(os.path.join(OutputDirectory, 'Normalization_Curve.txt'), NormCurve.T, header=NormalizationSignalName)

def ReadMetaDataFiles(InputDirectory):
    ''' Writes out meta data obtained while processing the inputs which may be useful.
    Parameters:
        InputDirectory (str): Name of directory that has the metadata files.

    Returns:
        Tilts (list of floats): A list of tilt angles.
        NormCurve (np.ndarray): Values to multiply by each tilt image in order to normalize the signal.
        NormalizationSignalName (str): Name of the signal that was used to normalize.
    '''
    # Read the tilts file
    print('Reading tilts file.')
    TiltsTemp = np.genfromtxt(os.path.join(InputDirectory,'tilts.txt'))
    Tilts = TiltsTemp[:,1]

    # Read the normalization curve.
    print('Reading normalization curve.')
    NormCurve = np.genfromtxt(os.path.join(InputDirectory, 'Normalization_Curve.txt'), names=True)
    NormalizationSignalName = NormCurve.dtype.names[0]

    return Tilts, NormCurve, NormalizationSignalName

def WriteSignalsToGENFIRE(OutputDirectory, SignalDict, Tilts, GENFIRETemplateDir=''):
    ''' Write npy files that GENFIRE can read and ensure all the metadata it needs are present.

    Parameters:
        OutputDirectory (str): Name of directory to place the tif files.

        SignalDict (OrderedDict of np.ndarray): Dictionary with signal stacks.  Each entry is one tilt sequence with dimension (tilt, x, y).  The key is the name of the signal, e.g. 'Fe_Ka' or 'HAADF'.

    Returns:
        A npy file for each signal.
        A script file for each signal to run GENFIRE using DoGenfire.py.
        DoGenfire.py is present in the output directory.
        A script file to run all the script files.
    '''

    # For now, rotation axis is hard coded.  We actually need to address the full rotation issue later.
    # Right now, I just assume the user did a perfect job of aligning the EDS stack around the stage alpha axis.  >;-)
    RotationAxisVertical = True # If false, then the rotation axis is going horizontally across the image.

    # Right now, we only support unix-like systems.  This currently runs on OS-X and Scientific Linux 7.  But it could easily be extended to other OS options. 

    # Read in the template for the script file to run GENFIRE.
    with open(os.path.join(os.path.dirname(__file__), 'slurm_template.sh'), 'r') as f:
        slurmscripttemplate = f.read()

    # Make sure the output directory has the DoGenfire.py file.
    shutil.copyfile(os.path.join(os.path.dirname(__file__), 'DoGenfire.py'), os.path.join(OutputDirectory, 'DoGenfire.py'))

    # Start a script file that will run all the script files.
    runall = open(os.path.join(OutputDirectory,'runall.sh'), 'w')

    # Write a npy file for each signal.
    for SigName, Sig in SignalDict.items():

        # Swap the numpy axes to be how GENFire expects it and save npy files.
        Sig = np.swapaxes(np.swapaxes(Sig, 0,1), 1,2)
        if(RotationAxisVertical==True):
            Sig = np.swapaxes(Sig, 0,1)
        print('Writing '+str(SigName)+'_aligned.npy, ', end='')
        np.save(os.path.join(OutputDirectory,SigName+'_aligned.npy'), np.array(Sig))

        # Write out a script file which will run GENFIRE for this MRC.
        # HAADF images are usually too big to do the full genfire.  We have avoid using oversampling so they fit in memory.
        if SigName == 'HAADF':
            osetting = '-o 1'
        else:
            osetting = ''

        slurmscript = slurmscripttemplate.replace('$JOBNAME', 'GENFIRE_'+SigName)
        slurmscript = slurmscript.replace('$GENFIRESTRING', SigName+' '+osetting+'\n')
        with open(os.path.join(OutputDirectory,SigName+'_slurm.sh'), 'w') as f:
            print('Writing '+SigName+'_slurm.sh')
            f.write(slurmscript)
        
        runall.write('sbatch '+SigName+'_slurm.sh\n')
        runall.write('sleep 3\n')
    print('Writing runall.sh')
    runall.close()

if __name__ == '__main__':
    print('------------------------------ Begin preprocess.py test ------------------------------')


    # InputDirectory = os.path.join('..', 'data', 'L2083-K-4-1', 'Input')
    # OutputDirectory = os.path.join('..', 'data', 'L2083-K-4-1', 'TestOutput', 'Unaligned')
    # AlignedDirectory = os.path.join('..', 'data', 'L2083-K-4-1', 'TestOutput', 'Aligned')
    OutputDirectory = os.path.join('..', 'data', 'L2083-K-4-1', 'TestOutput')
    SignalNames = ['HAADF', 'Al_Ka', 'C_Ka', 'Ca_Ka', 'Cr_Ka', 'Fe_Ka', 'Ga_Ka', 'Mg_Ka', 'Na_Ka', 'Ni_Ka', 'O_Ka', 'P_Ka', 'Pt_La', 'S_Ka', 'Si_Ka']
    EMDFileName = os.path.join('..', 'data', 'L2083-K-4-1', 'ConvertBrukerToEMD', 'Output', 'TomoData.emd')

    Signals, Tilts = ExtractSignalsFromEMD(EMDFileName, SignalNames=['HAADF', 'Fe_Ka', 'Si_K', 'Mg_K'])
    Signals, NormCurve = NormalizeSignals(Signals, Tilts, NormalizationSignalName='Si_K', NormalizationImageFraction=0.5)
    # WriteSignalsToTIFFs(OutputDirectory, Signals)

    # # Read in the data and process it up to the point where we need to do stack alignment in an external program.
    # if False:
    #     Signals, Tilts = ExtractRawSignalsFromBrukerSequence(InputDirectory=InputDirectory, SignalNames=SignalNames)
    #     print('Signal Names: '+str(Signals.keys()))
    #     print('HAADF signal shape: '+str(Signals["HAADF"].shape))
    #     print('Fe_Ka signal shape: '+str(Signals["Fe_Ka"].shape))
    #     print('Tilts: '+str(Tilts))
    #     Signals, NormCurve = NormalizeSignals(Signals, Tilts, NormalizationSignalName='Fe_Ka', NormalizationImageFraction=0.5)
    #     print(NormCurve)
    #     WriteSignalsToTIFFs(OutputDirectory, Signals)
    #     WriteMetaDataFiles(OutputDirectory, Tilts, NormCurve, NormalizationSignalName='Fe_Ka')
    #     input('Please do the stack alignment in ImageJ with the Multistackreg plugin now.  Press any key to continue...')

    # # At this point you should align the stacks.

    # # Apply the alignment file to generate aligned stacks and produce input into GENFIRE.
    # if True:
    #     SignalDict = ReadSignalsFromTIFFs(OutputDirectory, SignalNames=SignalNames)
    #     Tilts, NormCurve, NormalizationSignalName = ReadMetaDataFiles(OutputDirectory)
    #     Translations = ReadTomVizTranslations(os.path.join(OutputDirectory, 'TomVizAlignments.json'), Tilts)
    #     # Translations = ReadImageJTranslations(os.path.join(OutputDirectory, 'TransformationMatrices.txt'), Tilts)
    #     AlignedSignals = ApplyTranslations(SignalDict, Translations, 'HAADF')
    #     WriteSignalsToTIFFs(AlignedDirectory, AlignedSignals)
    #     WriteMetaDataFiles(AlignedDirectory, Tilts, NormCurve, NormalizationSignalName='Fe_Ka')
    #     WriteSignalsToGENFIRE(AlignedDirectory, AlignedSignals, Tilts)
    #     print('Now run GENFIRE using the (SignalName)_aligned.sh scripts.')

    # # At this point run GENFIRE.

    print('------------------------------ End preprocess.py test ------------------------------')
