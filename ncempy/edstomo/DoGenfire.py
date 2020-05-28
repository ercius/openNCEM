# import genfire as gf
from genfire import reconstruct, fileio
import sys, os, shutil
import argparse

# Parse command line arguments.
parser = argparse.ArgumentParser()
parser.add_argument("fileroot", help="There should be an aligned tomo stack named fileroot_aligned.npy in the directory and a the reconstruction will be named fileroot_reconstruction.mrc.")
parser.add_argument("-n", "--numIterations", help="Default number of GENFIRE iterations.  See GENFIRE documentation.  Default=100.", type=int, default=100, action="store")
parser.add_argument("-o", "--oversamplingRatio", help="How many times bigger to make the fourier space.  Bigger is more accurate, but you may run out of memory.  See GENFIRE documentation.  Default=4.", type=int, default=4, action="store")
args = parser.parse_args()

# Make a reconstruction using paramaters from default or from the command line.
GF = reconstruct.GenfireReconstructor(
        projections = args.fileroot+"_aligned.npy",
        eulerAngles = "tilts.txt",
        resultsFilename = args.fileroot+"_reconstruction.mrc",
        numIterations = args.numIterations,
        interpolationCutoffDistance=0.7,
        oversamplingRatio=args.oversamplingRatio,
        resolutionExtensionSuppressionState=2,
        )

# Number crunch.
results = GF.reconstruct()

# Save results.
fileio.saveResults(results, GF.params.resultsFilename)
