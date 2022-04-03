import argparse
import numpy as np
# import matplotlib.pyplot as plt
# import skimage.io
from cellpose import models, io
from PIL import Image

parser = argparse.ArgumentParser(description='Segments DAPI with cellpose as priors')
parser.add_argument('-i', '--input', nargs='+', help='input (list) of DAPI tiff files')
parser.add_argument('-o', '--output', nargs='+', help='output (list) of tiff files of segmentation masks')
parser.add_argument('--diameter', type=float, default=65)
parser.add_argument('-mask_threshold', type=float, default=0)
parser.add_argument('--crop', action='store_true', help='only run on a cropped image for testing')

def segment_image(model, dapi_file, mask_file, args):
    
    channel = [0,0]

    print('Loading %s' % filename)
    if args.crop:
        img = io.imread(filename)[4000:5000, 4000:5000]
    else:
        img = io.imread(filename)

    print('Segmenting %s' % filename)
    masks, flows, _, diams = model.eval(img, 
        diameter=args.diameter, mask_threshold=args.mask_threshold, channels=channel)

    # save results so you can load in gui
    print('Saving %s to npy file' % filename)
    io.masks_flows_to_seg(img, masks, flows, diams, filename, channel)

    print('Saving mask file %s to tiff for baysor')
    im = Image.fromarray(masks)
    im.save(mask_files[i])

    print('Done.\n\n')


if __name__ == "__main__":

    args = parser.parse_args()

    model = models.Cellpose(model_type='nuclei')

    dapi_files = args.input
    mask_files = args.output
    assert len(dapi_files) == len(mask_files)  

    for i,filename in enumerate(dapi_files):
        segment_image(model, filename, mask_file=mask_files[i], args=args)