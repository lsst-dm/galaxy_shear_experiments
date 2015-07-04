#!/usr/bin/env python
# encoding: utf-8
"""
generate_star_grid_instance_catalog.py

Make a PhoSim instance catalog with a grid of identical stars.

Many things copied from the DESC ImageSimulationRecipes simGridStars.py

Created by Michael Schneider on 2014-10-03
"""

import argparse
import itertools
import sys
import math
import random
import os.path
import subprocess, time
import numpy as np
import logging


# Print log messages to screen:
# logging.basicConfig(level=logging.DEBUG,
#                     format='%(asctime)s - %(levelname)s - %(message)s')
# Print log messages to file:
logging.basicConfig(filename='log_generate_star_grid_instance_catalog.log',
                    level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


###########################################################
# Utilities to be exported to python libraries
###########################################################
def gen_star_list(outfilename, type, fov, ra_cen, dec_cen, nstars_per_dim, seeing, filter, seed, histid):
	sed_filename = "../sky/sed_flat.txt"
        nstars = nstars_per_dim * nstars_per_dim
	### Pointing from simGridStars in ImageSimulationRecipes

	### Coordinate ranges (degrees)
	###  Space stars 1 arcminute apart

        raCorr = math.cos(dec_cen*math.pi/180.0)
        if type == "random":
	    ra = []
            dec = []
            dec_min = dec_cen - fov/2.0
            dec_max = dec_cen + fov/2.0
            ra_min = ra_cen - fov/2.0/raCorr
            ra_max = ra_cen + fov/2.0/raCorr

            for i in range(nstars):
                ra.append(random.uniform(ra_min, ra_max))
                dec.append(random.uniform(dec_min, dec_max))

        elif type == "grid":
	    grid_radius = np.min(((nstars_per_dim / 60.) / 2., fov/2.))
	    ra_min = ra_cen - (fov/2.0/raCorr)
	    ra_max = ra_cen + (fov/2.0/raCorr)
	    dec_min = dec_cen - fov/2.0
	    dec_max = dec_cen + fov/2.0
            ra, dec = np.meshgrid(np.linspace(ra_min, ra_max, num=nstars_per_dim),
						  np.linspace(dec_min, dec_max, num=nstars_per_dim))
	    ra = ra.ravel()
	    dec = dec.ravel()
        else:
            raise StandardException("Illegal type requested: %s"%type)

        #   Do a grid distribution
        steps = 10
        hist = np.zeros((steps, steps), dtype=int)
        rastep = (ra_max-ra_min)/steps
        decstep = (dec_max-dec_min)/steps
        for i in range(len(ra)):
            ra_bin = int((ra[i]-ra_min)/rastep)
            dec_bin = int((dec[i]-dec_min)/decstep)
            if ra_bin == steps and ra[i] == ra_max:
                ra_bin = steps - 1
            if dec_bin == steps and dec[i] == dec_max:
                dec_bin = steps - 1
            hist[ra_bin, dec_bin] = hist[ra_bin, dec_bin]+1
        print hist
        print hist.sum()
        fout = open(outfilename, "w")
	# stars = pd.DataFrame({"object": "object",
	# 	"id": np.arange(nstars_per_dim * nstars_per_dim, dtype=float),
	### Write header of instance catalog file with observing and simulation parameters
	print >> fout, 'Unrefracted_RA', ra_cen
	print >> fout, 'Unrefracted_Dec', dec_cen
	print >> fout, 'Opsim_obshistid %d'%histid
	print >> fout, 'Unrefracted_RA 81.1760879'
	print >> fout, 'Unrefracted_Dec -12.2103036'
	print >> fout, 'Opsim_moonra 313.898938'
	print >> fout, 'Opsim_moondec -12.7605726'
	print >> fout, 'Opsim_rotskypos 184.222551'
	print >> fout, 'Opsim_rottelpos 0'
	print >> fout, 'Opsim_filter %d'%filter
	print >> fout, 'Opsim_rawseeing %f'%seeing
	print >> fout, 'Opsim_sunalt -18.1759974'
	print >> fout, 'Opsim_moonalt -25.2607921'
	print >> fout, 'Opsim_dist2moon 122.048036'
	print >> fout, 'Opsim_moonphase 0.394149005'
	print >> fout, 'Opsim_expmjd 50486.0469'
	print >> fout, 'Opsim_altitude 72.4911951'
	print >> fout, 'Opsim_azimuth 355.24931'
	print >> fout, 'SIM_SEED %d'%seed
	print >> fout, 'SIM_EXPTIME 30.0'
	print >> fout, 'SIM_VISTIME 30.0'
	print >> fout, 'SIM_NSNAP 1'

	### Write sources to instance catalog file
	### Columns in the instance catalog:
	###   ["object", "id", "ra", "dec", "mag_norm", "sed_name",
	###	   "redshift", "gamma1", "gamma2", "kappa", "delta_ra", "delta_dec",
	###	   "source_type", "dust_rest_name", "dust_lab_name"]
	for i in xrange(len(ra)):
		print >> fout, 'object', str(i), str(ra[i]), str(dec[i]), '17.5', sed_filename, '0 0 0 0 0 0 star none none'
	fout.close()

	# stars = pd.DataFrame({"object": "object",
	# 	"id": np.arange(nstars_per_dim * nstars_per_dim, dtype=float),
	# 	"ra": ra.ravel(),
	# 	"dec": dec.ravel(),
	# 	"mag_norm": 30.,
	# 	"sed_name": sed_filename,
	# 	"redshift": 0.0,
	# 	"gamma1": 0,
	# 	"gamma2": 0,
	# 	"kappa": 0,
	# 	"delta_ra": 0,
	# 	"delta_dec": 0,
	# 	"source_type": "star",
	# 	"dust_rest_name": "none",
	# 	"dust_lab_name": "none"},
	# 	columns=["object", "id", "ra", "dec", "mag_norm", "sed_name",
	# 			 "redshift", "gamma1", "gamma2", "kappa", "delta_ra", "delta_dec",
	# 			 "source_type", "dust_rest_name", "dust_lab_name"])
	# stars.to_csv(outfile, sep=" ", header=False, index=False)
	return None


# def run_phosim():
# 	# And run phosim over this new file. Note that we've specified
# 	# no sky background, to speed things up.
# 	subcmd = ['python', phosim_path+"phosim.py", catfilename,  "-c", phosim_path+"examples/nobackground", "-s", sensor, "-w", workdir, "-o", outdir]

# 	subprocess.call(subcmd)


###########################################################
# Test functions for this script
###########################################################


def main():
    parser = argparse.ArgumentParser(
        description='My program description.')
    parser.add_argument('--outfile', default="star_grid.pars", help='Output directory')
    parser.add_argument('--nstars_per_dim', type=int, default=100,
    					help='Number of stars/dimension')
    parser.add_argument('--raw_seeing', type=float, default=.7, help='raw seeing')
    parser.add_argument('--filter', type=int, default=2, help='filter number u:0 g:1 r:2 i:3 z:4 y:5')
    parser.add_argument('--random_seed', type=long, default=12345678, help='random number generator seed')
    parser.add_argument('--histid', type=int, default=12345, help='histid')
    parser.add_argument('--fov', type=float, default=3.5, help='field of view')
    parser.add_argument('--ra_cen', type=float, default=81.176, help='ra at focal plane center')
    parser.add_argument('--dec_cen', type=float, default=-12.210, help='dec at focal plane center')
    parser.add_argument('--type', default="grid", help='specify whether output is grid or random')

    args = parser.parse_args()
    logging.debug('Started generation of star grid catalog')
    outfile = os.path.abspath(args.outfile)
    logging.info('Saving star instance catalog to %s' % outfile)
    gen_star_list(outfile, args.type, args.fov, args.ra_cen, args.dec_cen, args.nstars_per_dim, args.raw_seeing, args.filter, args.random_seed, args.histid)
    logging.debug('Finished generation of star grid catalog')   
    return 0


if __name__ == "__main__":
    sys.exit(main())
