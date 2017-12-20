import argparse
from ncempy.io.dm import fileDM
from matplotlib import cm
from matplotlib.image import imsave

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
    raise ValueError("Extension/filetype {} not supported!".format(extension))

if extension in ["dm3","dm4"]:
    f = fileDM(source_file, on_memory=True)
    f.parseHeader()
    ds = f.getDataset(0)
    img = ds['data']
    if len(img.shape)>2:
        img = img[int(img.shape[0]/2),:,:]
    imsave(dest_file, img, format="png", cmap=cm.gray)