#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2015 AURA/LSST
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import argparse
import pyfits
import sys
import os.path
import numpy
import lsst.daf.persistence as dafPersist
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.coord as afwCoord
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.meas.base as measBase
import lsst.meas.algorithms as measAlg

# Script to combine a directory of many psf cutouts into a single
# multiple extension fits file

def main(args):
    draw_psf_src = args.input_directory
    listing = os.listdir(draw_psf_src)
    psf_numbers = []
    hdus = None
    
    # Search for appropriately names psfs_nnnn.fits files in the psf directory
    for file in listing:
        index1 = file.find("psfs_")
        index2 = file.find(".fits")
        if index1 != 0 or index2 < 0: continue
        number = file[index1+5:index2]
        if number.isdigit():
            psf_numbers.append(int(number))
    if len(psf_numbers) < 2:
        print "Fewer than 2 psfs_[0-9]*.fits files found. Not creating package"
        sys.exit(1)
    hdus = pyfits.open("%s/psfs_%d.fits" % (draw_psf_src, psf_numbers[0]))
    for psf_number in psf_numbers[1:]:
       hdus.append(pyfits.open("%s/psfs_%d.fits"%(draw_psf_src, psf_number))[0])
    
    # Add the original number to the header of each hdu
    for i in range(len(hdus)):
        hdus[i].header.append(("PSF_NO",psf_numbers[i]))
    
    # Write the multiple extension fits file to the output file
    hdus.writeto(args.output_file, clobber=True)
    print "%s file created with %d files"%(args.output_file, len(psf_numbers))
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_directory", help="Name of the input directory")
    parser.add_argument("output_file", help="Name of multi-psf fits file")
    args = parser.parse_args()
    main(args)    
