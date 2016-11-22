'''
Command line tool to convert a SER file to an EMD file.
'''

import argparse
import ncempy.fio.ser

# parse the commandline arguments
parser = argparse.ArgumentParser(description='Tool to convert a SER file to an EMD file.')
parser.add_argument('input', help='path to the input SER file')
parser.add_argument('--emi', help='path to corresponding EMI file')
parser.add_argument('output', help='path to the output EMD file')
args = parser.parse_args()

fser = ncempy.fio.ser.fileSER(args.input, emifile=args.emi)
fser.writeEMD(args.output)
