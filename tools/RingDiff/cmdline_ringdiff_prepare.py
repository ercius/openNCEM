'''
Command line tool to prepare an emd evaluation file for diffraction ring analysis.
'''

import argparse
import os
import ncempy.io.emd
import ncempy.eval.ring_diff

# parse the commandline arguments
parser = argparse.ArgumentParser(description='This tool prepares an emd evaluation file for diffraction ring analysis. The file is created with a blank parental settings group. Evaluation groups for all emdgroups in the given emdfiles are created. The emd evaluation file is meant to be edited with an external hdf5 viewer before the evaluation is executed.')
parser.add_argument('-i','--input', nargs='*', help='one or more emdfiles with emdgroups to evaluate')
parser.add_argument('-o','--output', nargs=1, help='path to the output EMD file (will be overwritten)')
args = parser.parse_args()

print('{}'.format(args.output))
print('{}'.format(args.input))

# open output emdfile
if os.path.isfile(args.output[0]):
    os.remove(args.output[0])
femd_out = ncempy.io.emd.fileEMD(args.output[0])
grp_eva = femd_out.file_hdl.create_group('evaluation')

print('Evaluation file {} created.'.format(args.output[0]))

# insert dummy settings
grp_set = ncempy.eval.ring_diff.put_settings(grp_eva, ncempy.eval.ring_diff.dummy_settings)

print('Dummy settings written to {}.'.format(grp_set.name))
print('.. to be edited with external hdf5 viewer. Note that you can copy this settings group and all attributes to the evaluation subgroups, if you want to use customized settings for single evaluations. A number of settings is optional and the corresponding attributes/datasets can be deleted to use default values during evaluation: plt_imgminmax, rad_rmax, rad_dr, rad_sigma, mask, fit_maxfev.')

# gather all emdtype groups from input files
if not args.input is None:

    for fname in args.input:
    
        print('Collecting emdtype groups from {}.'.format(fname))
        
        hdl = grp_eva.create_group(os.path.basename(fname))
        
        femd_in = ncempy.io.emd.fileEMD(fname, readonly=True)

        for i in range(len(femd_in.list_emds)):
        
            ncempy.eval.ring_diff.put_sglgroup(hdl, '{}'.format(femd_in.list_emds[i].name.split('/')[-1]), femd_in.list_emds[i])
            print('.. {}'.format(femd_in.list_emds[i].name))

print('Remember to delete evaluation groups for emdtype groups which are not intended to be evaluated.')

del femd_out, femd_in

