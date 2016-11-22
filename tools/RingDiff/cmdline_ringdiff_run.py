'''
Command line tool to run diffraction ring analysis on an emd evaluation file.
'''

import argparse
import os
import ncempy.fio.emd
import ncempy.eval.ring_diff

# parse the commandline arguments
parser = argparse.ArgumentParser(description='This tool runs the ring diffraction analysis for a given emd evaluation file.')
parser.add_argument('input', help='path to the emd evaluation file')
parser.add_argument('-f', action='store_true', help='overwrite existing results')
parser.add_argument('-v', action='store_true', help='verbose mode')
parser.add_argument('-p', action='store_true', help='show the plots')

args = parser.parse_args()


# open input file
femd = ncempy.fio.emd.fileEMD(args.input)

# execute
ncempy.eval.ring_diff.run_all(femd.file_hdl, femd, args.f, args.v, args.p)

# close output    
del femd

# wait for plots
if args.p:
    import matplotlib.pyplot as plt
    plt.show()

