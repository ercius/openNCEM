import numpy as np
import logging

from ncempy.io.ser import serReader

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
	pass