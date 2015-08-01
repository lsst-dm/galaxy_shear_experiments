from great3sims import constants, run
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("base", help="name of the output directory")
parser.add_argument("-p", "--galaxy_stamp_size", help="Size of galaxy images to be cutout", type = int, default=64)
parser.add_argument("-d", "--ndims", help="each image file is ndims x ndims galaxies", type = int, default=100)
parser.add_argument("-n", "--n_subfields", help="number of subfields to process", type = int, default=1)
parser.add_argument("-e", "--n_epochs", help="number of epochs to process", type = int, default=1)
parser.add_argument("-s", "--shear_value", help="shear value, if fixed", type = float, default=None)
parser.add_argument("-a", "--shear_angle", help="shear angle, if fixed", type = float, default=None)
parser.add_argument("-f", "--fov", help="field of view, in degrees", type = float, default=1.50)
parser.add_argument("-g", "--get_arg", help="return argument value", type = str, default=None)
parser.add_argument("-t", "--exp_type", help="type of great3sims experiment", type = str, default="control")
parser.add_argument("-r", "--gal_dir", help="location of galaxy data", type = str, default="control")

args = parser.parse_args()
if not args.get_arg == None:
    print eval(args.get_arg, args.__dict__)
    sys.exit(0)

constants.image_size_deg = args.fov
constants.nrows = args.ndims
constants.ncols = args.ndims
constants.n_subfields =  args.n_subfields
constants.xsize["ground"][True] = args.galaxy_stamp_size
constants.xsize["ground"][False] = args.galaxy_stamp_size
constants.n_subfields_per_field["constant"][True] = 1
constants.subfield_grid_subsampling = 1
constants.n_deep_subfields = 0
constants.deep_frac = 0.0
subfield_max = constants.n_subfields + constants.n_deep_subfields - 1
run(args.base, gal_dir=args.gal_dir, steps=["metaparameters", "catalogs", "config",], experiments=[args.exp_type], obs_type="ground", shear_type=["constant"], draw_psf_src = '%s/control/ground/constant/psfs/psfs.index'%args.base, subfield_max=subfield_max, shear_value=args.shear_value, shear_angle=args.shear_angle)
