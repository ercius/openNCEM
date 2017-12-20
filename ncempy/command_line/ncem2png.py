""" This scripts converts SER, DM3, and DM4 files into PNG """

import argparse
from ncempy.io.dm import fileDM
from ncempy.io.ser import fileSER
import ntpath
import os

from matplotlib import cm
from matplotlib.image import imsave

def _discover_emi(file_route):
    file_name = ntpath.basename(file_route)
    folder = os.path.dirname(file_route)
    parts = file_name.split("_")
    if len(parts)==1:
        # No "_" in the filename.
        return None
    emi_file_name = "{}.emi".format("_".join(parts[:-1]))
    emi_file_route = os.path.join(folder, emi_file_name)
    if not os.path.isfile(emi_file_route):
        # file does not exist
        return None
    
    return emi_file_route

def dm_to_png(source_file, dest_file):
    f = fileDM(source_file, on_memory=True)
    f.parseHeader()
    ds = f.getDataset(0)
    img = ds['data']
    dimension=len(img.shape)
    print(img.shape)
    if dimension == 3:
        img = img[:,:,int(img.shape[2]/2),]
    elif dimension == 4:
        img = img[:,:,
                  int(img.shape[2]/2),
                  int(img.shape[3]/2)]
    else:
        raise ValueError("This scripts cannot extract PNGs from DM files with"
                         " more than four dimensions.")
    imsave(dest_file, img, format="png", cmap=cm.gray)
    return f

def ser_to_png(source_file, dest_file):
    emi_file = _discover_emi(source_file)
    f = fileSER(source_file, emi_file)
    ds = f.getDataset(0)
    img = ds[0]
    imsave(dest_file, img, format="png")
    return f
    
def main():
    parser = argparse.ArgumentParser(description='Extracts a preview png from'
                                     ' a SER, DM3, or DM4 file.')
    
    parser.add_argument('source_file', metavar='source_file', type=str, nargs=1,
                    help='Source file, must have ser, dm3, o dm4 extension.')
    
    parser.add_argument('dest_file', metavar='dest_file', type=str, nargs='?',
                        help='Filename to write the preview png file. If not set'
                        ' it defaults to the name of the source file appending'
                        ' png.',
                        default=None)
    
    args = parser.parse_args()
    
    source_file = args.source_file[0]
    dest_file = args.dest_file
    if dest_file is None:
        dest_file="{}.png".format(source_file)
    
    extension = source_file.split(".")[-1].lower()
    if  not extension in ["dm3", "dm4", "ser"]:
        raise ValueError("Extension/filetype {} not supported!".format(
                                                            extension))
    
    if extension in ["dm3","dm4"]:
        dm_to_png(source_file, dest_file)
    
    if extension in ["ser"]:
        ser_to_png(source_file, dest_file)

if __name__ =="__main__":
    main()