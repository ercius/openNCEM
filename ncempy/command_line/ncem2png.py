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

def extract_dimension(img, fixed_dimensions=None):
    out_img=img

    dimension=len(img.shape)
    if fixed_dimensions is None:
        if dimension == 3:
            fixed_dimensions=[int(img.shape[0]/2),
                              '','']
        elif dimension == 4:
            fixed_dimensions=['',
                              int(img.shape[2]/2),
                              int(img.shape[3]/2)]
        elif dimension > 4:
            raise ValueError("This scripts cannot extract PNGs from DM files"
                            " with"
                             " more than four dimensions without explicit"
                             " use of --fixed_dimensions.")

    if len(fixed_dimensions)!=dimension:
        raise ValueError("Number of values to index the image ({}) do not"
                         " match the image dimension ({}).".format(
                             len(fixed_dimensions), dimension))
    d_tuple=()
    selected_dimensions=0
    print("{}D Selecting frame ({})".format(dimension,
                                            fixed_dimensions))
    for (d, s) in zip(fixed_dimensions, img.shape):
        if not d:
            d_tuple+=(slice(None, None, None),)
        elif d=="m":
            d_tuple+=(int(s/2),)
            selected_dimensions+=1
        elif int(d)<s:
            selected_dimensions+=1
            d_tuple+=(int(d),)
    if dimension-selected_dimensions!=2:
        raise ValueError("Dimension extractor do not fix enough components"
                         " to produce a 2D image. Needs to fix: {}"
                         "".format(dimension-2))
    out_img = img[d_tuple] 
    return out_img

def dm_to_png(source_file, dest_file, fixed_dimensions=None):
    """ Saves the DM3 or DM4 source_file as PNG dest_file. If the data has three
    of four dimensions. The image taken is from the middle image in those
    dimensions."""
    f = fileDM(source_file, on_memory=True)
    f.parseHeader()
    ds = f.getDataset(0)
    img = ds['data']
    img = extract_dimension(img, fixed_dimensions)
    imsave(dest_file, img, format="png", cmap=cm.gray)
    return f

def ser_to_png(source_file, dest_file):
    """ Saves the SER source_file as PNG dest_file."""
    emi_file = _discover_emi(source_file)
    f = fileSER(source_file, emi_file)
    ds = f.getDataset(0)
    img = ds[0]
    imsave(dest_file, img, format="png")
    return f
    
def main():
    parser = argparse.ArgumentParser(description='Extracts a preview png from'
                                     ' a SER, DM3, or DM4 file.')
    
    parser.add_argument('source_files', metavar='source_files', type=str,
                        nargs="+",
                    help='Source files, must have ser, dm3, o dm4 extension.')
    
    parser.add_argument('--out_file', dest='dest_file', action='store', 
                        type=str, nargs=1,
                        help='Filename to write the preview png file. If not set'
                        ' it defaults to the name of the source file appending'
                        ' png. Only valid when a single input file is '
                        ' processed.',
                        default=None)
    
    parser.add_argument('--fixed_dimensions', dest='fixed_dimensions',
                        action='store', 
                        type=str, nargs=1,
                        help="List of numbers, 'm'(for middle), or nothing to"
                        " extract"
                        " a particular slice of data. e.g., '2,,' will extract"
                        " an png with the values y,z of x=2. ',m,,2' will "
                        " extract all the values x,z for y=1/2shapeY, and w=2.",
                        default=None)
    
    args = parser.parse_args()
    
    
    if args.dest_file is not None and len(args.source_files)>1:
        raise ValueError("--out_file only can be used when a single input file"
                         " is processed.")
    
    fixed_dimensions = args.fixed_dimensions
    if fixed_dimensions is not None:
        fixed_dimensions=fixed_dimensions[0].split(',')
    
    for source_file in args.source_files:
        
        if args.dest_file is None:
            dest_file="{}.png".format(source_file)
        else:
            dest_file=args.dest_file[0]
        extension = source_file.split(".")[-1].lower()
        if  not extension in ["dm3", "dm4", "ser"]:
            raise ValueError("Extension/filetype {} not supported!".format(
                                                                extension))
        print("Extracting from {}, saving image as {}".format(source_file,
                                                      dest_file ))
        if extension in ["dm3","dm4"]:
            dm_to_png(source_file, dest_file, fixed_dimensions=fixed_dimensions)
        
        if extension in ["ser"]:
            ser_to_png(source_file, dest_file)

if __name__ =="__main__":
    main()