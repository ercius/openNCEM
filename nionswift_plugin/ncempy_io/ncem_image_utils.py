import numpy as np
import logging

from ncempy.io.ser import serReader
from ncempy.io.mrc import mrcReader
from ncempy.io.emd import fileEMD
from ncempy.io.emdVelox import fileEMDVelox

from nion.data.Calibration import Calibration
from nion.data import DataAndMetadata


def loadSER(fpath):
	logging.debug('Reading ser file with openNCEM')

	ser = serReader(fpath)

	data = ser['data']

	dimensional_calibrations = [Calibration(ser['pixelOrigin'][i],ser['pixelSize'][i],ser['pixelUnit'][i]) for i in range(len(ser['pixelSize']))]
	
	# currently only supporting 2D images and 3D sequences
	if len(data.shape) == 2:
		descriptor = DataAndMetadata.DataDescriptor(False,0,2)
	if len(data.shape) == 3:
		descriptor = DataAndMetadata.DataDescriptor(False,1,2)
	else:
		descriptor = DataAndMetadata.DataDescriptor(False,0,len(data.shape))

	return DataAndMetadata.new_data_and_metadata(data,dimensional_calibrations=dimensional_calibrations,data_descriptor=descriptor)

def loadMRC(fpath):
	pass





def loadEMD(fpath):
	# first try to open as a properly formatted EMD file (Berkeley specification)
	emd = fileEMD(fpath,readonly=True)

	# check if there are any valid EMD datasets:
	if len(emd.list_emds) > 0:
		# the file is a Berkeley EMD. Currently supporting only first dataset
		data = np.array(emd.list_emds[0]['data'])

		if len(data.shape) == 2:
			descriptor = DataAndMetadata.DataDescriptor(False,0,2)
		if len(data.shape) == 3:
			descriptor = DataAndMetadata.DataDescriptor(False,1,2)
		if len(data.shape) == 4:
			descriptor = DataAndMetadata.DataDescriptor(False,2,2)
		if len(data.shape) == 5:
			descriptor = DataAndMetadata.DataDescriptor(True,2,2)
			data = np.moveaxis(data,4,0)
			data = np.swapaxes(data,1,3)
			data = np.swapaxes(data,2,4)
		else:
			descriptor = DataAndMetadata.DataDescriptor(False,0,len(data.shape))

		origin = []
		pixsize = []
		units = []
		for j in range(len(data.shape)):
			dim = emd.list_emds[0]['dim'+str(j+1)]
			if isinstance(dim[0],bytes):
				# this axis is a complex datatype, so use 0 as the origin, 1 as the step, and imag as units
				origin.append(0)
				pixsize.append(1)
				units.append('imag')
			else:
				origin.append(dim[0])
				pixsize.append( dim[1]-dim[0])
				units.append('pixels')


		dimensional_calibrations = [Calibration(origin[i],pixsize[i],units[i]) for i in range(len(origin))]

		return DataAndMetadata.new_data_and_metadata(data,dimensional_calibrations=dimensional_calibrations,data_descriptor=descriptor)		

