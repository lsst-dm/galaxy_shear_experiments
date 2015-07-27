#!/usr/bin/env python
# encoding: utf-8
"""
generate_star_grid_instance_catalog.py
Make a PhoSim instance catalog with a grid of identical stars.
Many things copied from the DESC ImageSimulationRecipes simGridStars.py
Created by Michael Schneider on 2014-10-03
"""

import argparse
import sys
import os.path
import subprocess, time
import numpy as np
# import pandas as pd
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
def gen_star_list(outfilename, nstars_per_dim=10):
	sed_filename = "../sky/sed_flat.txt"

	### Pointing from simGridStars in ImageSimulationRecipes
	ra_cen = 81.176
	dec_cen = -12.210

	### Coordinate ranges (degrees)
	###  Space stars 1 arcminute apart
	fov = 3.5  
	grid_radius = np.min(((nstars_per_dim / 60.) / 2., fov/2.))
	ra_min = ra_cen - grid_radius
	ra_max = ra_cen + grid_radius
	dec_min = dec_cen - grid_radius
	dec_max = dec_cen + grid_radius
	# ra_min = ra_cen - 0.5 * fov
	# ra_max = ra_cen + 0.5 * fov
	# dec_min = dec_cen - 0.5 * fov
	# dec_max = dec_cen + 0.5 * fov

	ra, dec = np.meshgrid(np.linspace(ra_min, ra_max, num=nstars_per_dim),
						  np.linspace(dec_min, dec_max, num=nstars_per_dim))
	ra = ra.ravel()
	dec = dec.ravel()

	### Write header of instance catalog file with observing and simulation parameters
	outfile = open(outfilename, 'w')
	print >> outfile, 'Unrefracted_RA', ra_cen
	print >> outfile, 'Unrefracted_Dec', dec_cen
	print >> outfile, 'Opsim_obshistid 111111'
	print >> outfile, 'SIM_SEED 97895167'
	print >> outfile, 'SIM_VISTIME 15.0'
	print >> outfile, 'Unrefracted_RA 81.1760879'
	print >> outfile, 'Unrefracted_Dec -12.2103036'
	print >> outfile, 'Opsim_moonra 313.898938'
	print >> outfile, 'Opsim_moondec -12.7605726'
	print >> outfile, 'Opsim_rotskypos 184.222551'
	print >> outfile, 'Opsim_rottelpos 0'
	print >> outfile, 'Opsim_filter 2'
	print >> outfile, 'Opsim_rawseeing 0.6'
	print >> outfile, 'Opsim_sunalt -18.1759974'
	print >> outfile, 'Opsim_moonalt -25.2607921'
	print >> outfile, 'Opsim_dist2moon 122.048036'
	print >> outfile, 'Opsim_moonphase 0.394149005'
	print >> outfile, 'Opsim_expmjd 50486.0469'
	print >> outfile, 'Opsim_altitude 72.4911951'
	print >> outfile, 'Opsim_azimuth 355.24931'

	### Write sources to instance catalog file
	### Columns in the instance catalog:
	###   ["object", "id", "ra", "dec", "mag_norm", "sed_name",
	###	   "redshift", "gamma1", "gamma2", "kappa", "delta_ra", "delta_dec",
	###	   "source_type", "dust_rest_name", "dust_lab_name"]
	for i in xrange(len(ra)):
		print >> outfile, 'object', str(i), str(ra[i]), str(dec[i]), '17.5', sed_filename, '0 0 0 0 0 0 star none none'
	outfile.close()

	# stars = pd.DataFrame({"object": "object",
	#	"id": np.arange(nstars_per_dim * nstars_per_dim, dtype=float),
	#	"ra": ra.ravel(),
	#	"dec": dec.ravel(),
	#	"mag_norm": 30.,
	#	"sed_name": sed_filename,
	#	"redshift": 0.0,
	#	"gamma1": 0,
	#	"gamma2": 0,
	#	"kappa": 0,
	#	"delta_ra": 0,
	#	"delta_dec": 0,
	#	"source_type": "star",
	#	"dust_rest_name": "none",
	#	"dust_lab_name": "none"},
	#	columns=["object", "id", "ra", "dec", "mag_norm", "sed_name",
	#			 "redshift", "gamma1", "gamma2", "kappa", "delta_ra", "delta_dec",
	#			 "source_type", "dust_rest_name", "dust_lab_name"])
	# stars.to_csv(outfile, sep=" ", header=False, index=False)
	return None


# def run_phosim():
#	# And run phosim over this new file. Note that we've specified 
#	# no sky background, to speed things up. 
#	subcmd = ['python', phosim_path+"phosim.py", catfilename,  "-c", phosim_path+"examples/nobackground", "-s", sensor, "-w", workdir, "-o", outdir]

#	subprocess.call(subcmd)


###########################################################
# Test functions for this script
###########################################################


def main():
    parser = argparse.ArgumentParser(
        description='My program description.')
    parser.add_argument('--outfile', default="star_grid.pars", help='Output directory')
    parser.add_argument('--nstars_per_dim', type=int, default=10,
					help='Number of stars per dimension in output grid')

    args = parser.parse_args()
    logging.debug('Started generation of star grid catalog')
    outfile = os.path.abspath(args.outfile)
    logging.info('Saving star instance catalog to %s' % outfile)
    gen_star_list(outfile, nstars_per_dim=args.nstars_per_dim)
    logging.debug('Finished generation of star grid catalog')    
    return 0


if __name__ == "__main__":
    sys.exit(main())
