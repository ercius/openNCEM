import sys, os, shutil
from collections import OrderedDict
import glob2 as glob
import numpy as np
import hyperspy.api as hs
from ncempy.io.emd import fileEMD

def GetTiltsFromBrukerSequence(Directory=None):
    ''' Find the set of Bruker bcf files and figure out their tilts based on filenames.

    There should be a set of bcf files in a directory with names: ... -10.bcf, -5.bcf, 0.bcf, 5.bcf ...
    The number represents the stage tilt of that stack.
    The angles do not matter but the names must be coded as above.  They will be automatically read and sorted.

    Parameters:
        Directory (str): Name of directory containing the bcf files.

    Returns:
        (list of floats): A list of tilt angles.
    '''

    Bcfs = glob.glob(os.path.join(Directory, '*.bcf'))
    NewBcfs = []
    for i in Bcfs:
        _ , FileName = os.path.split(i)
        NewBcfs.append(FileName)
    Tilts = list(map(lambda s: float(s[:-4]), NewBcfs))
    Tilts.sort()

    return Tilts

def GetSpatialDimension(Dim):
    ''' Return a spatial value in meters. '''
    unit = Dim.units
    if unit == 'nm':
        scalefac = 1e-9
    elif unit == 'µm':
        scalefac = 1e-6
    elif unit == 'mm':
        scalefac = 1e-3
    elif unit == 'cm':
        scalefac = 1e-2
    else:
        print('Unrecognized spatial dimension unit in Bruker file: ' + unit)
        scalefac = 1
    Axis = Dim.index2value(range(Dim.size))
    return Axis * scalefac

def GetEnergyDimension(Dim):
    ''' Return an energy value in eV.  '''
    unit = Dim.units
    if str(unit) == 'eV':
        scalefac = 1
    elif str(unit) == str('keV'):
        scalefac = 1e3
    elif unit == 'MeV':
        scalefac = 1e6
    else:
        print('Unrecognized energy dimension unit in Bruker file: ' + unit)
        scalefac = 1
    Axis = Dim.index2value(range(Dim.size))
    return Axis * scalefac

def ExtractRawSignalsFromBrukerSequence(InputDirectory=None, OutputEMD=None):
    ''' Read in a set of Bruker bcf files containing EDS acquisitions and write an EMD file with the data.

    Parameters:
        InputDirectory (str): Name of directory containing the bcf files.

        OutputEMD (str): Name of the output file (may include a path).

    Returns:
        None.
    '''

    # Find all the bcf files first.
    Tilts = GetTiltsFromBrukerSequence(Directory=InputDirectory)
    Tilts = np.array(Tilts)

    # We need some filename if none was given to us.
    if OutputEMD is None:
        OutputEMD = 'test.emd'
        
    print('Extracting Signals:')
    for i, t in enumerate(Tilts):
        # Load the bruker file for this tilt.
        fname = os.path.join(InputDirectory, str(int(t))+'.bcf')
        x = hs.load(fname)

        # Only the first time, we need to calculate the sizes of all the arrays we are going to make.
        if 'HAADFsize_m' not in locals():
            # HAADF's have x,y dimensions.
            HAADFsize_m = GetSpatialDimension(x[0].axes_manager['width'])
            HAADFsize_n = GetSpatialDimension(x[0].axes_manager['height'])
            HAADFDim = (len(Tilts), len(HAADFsize_m), len(HAADFsize_n))
            print('HAADF has dimensions ' + str(HAADFDim))

            # EDS cubes have x,y, energy dimensions.
            EDSsize_m = GetSpatialDimension(x[1].axes_manager['width'])
            EDSsize_n = GetSpatialDimension(x[1].axes_manager['height'])
            EDSsize_e = GetEnergyDimension(x[1].axes_manager['Energy'])
            EDSDim = (len(Tilts), len(EDSsize_m), len(EDSsize_n), len(EDSsize_e))
            print('EDS has dimensions ' + str(EDSDim))

            HAADF = np.zeros(HAADFDim)
            EDS = np.zeros(EDSDim, dtype='float32')

            BeamEnergy = x[0].metadata['Acquisition_instrument']['TEM']['beam_energy']*1000 # eV
            EnergyResolutionMnKa = x[1].metadata['Acquisition_instrument']['TEM']['Detector']['EDS']['energy_resolution_MnKa'] # eV
            DetectorTiltAngle = x[1].metadata['Acquisition_instrument']['TEM']['Detector']['EDS']['elevation_angle'] # degrees
            RealTime = x[1].metadata['Acquisition_instrument']['TEM']['Detector']['EDS']['real_time'] # seconds

        print(str(fname))

        HAADF[i,:,:] = x[0].data.astype("float32")
        EDS[i,:,:,:] = x[1].data.astype("float32") 

    # open nonexisting file for writing
    if os.path.isfile(OutputEMD):
        os.remove(OutputEMD)
    EMD = fileEMD(OutputEMD)

    data = HAADF
    dims = ( (Tilts, 'angle', '[deg]'),
             (HAADFsize_m, 'x', '[m]'),
             (HAADFsize_n, 'y', '[m]'))
    print('Writing HAADF tilt stack.')
    EMD.put_emdgroup('HAADF_TiltStack', data, dims)

    data = EDS
    dims = ( (Tilts, 'angle', '[deg]'),
             (EDSsize_m, 'x', '[m]'),
             (EDSsize_n, 'y', '[m]'),
             (EDSsize_e, 'E', '[eV]')) 
    print('Writing EDS tilt stack.')
    EMD.put_emdgroup('EDS_TiltStack', data, dims)

    EMD.microscope.attrs['BeamVoltage[eV]'] = BeamEnergy
    EMD.microscope.attrs['MnKaResolution[eV]'] = EnergyResolutionMnKa
    EMD.microscope.attrs['DetectorTiltAngle[deg]'] = DetectorTiltAngle
    EMD.microscope.attrs['RealTime[s]'] = RealTime * len(Tilts)

    EMD.put_comment('File created.')
    del EMD
    print('Created file ' + OutputEMD)

