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

parser = argparse.ArgumentParser()
parser.add_argument("output_file", help="Name of multi-psf fits file")
parser.add_argument("-d", "--delete", help="delete individual psf.fits files", action="store_true")
args = parser.parse_args()

draw_psf_src = "./"
listing = os.listdir(draw_psf_src)
psf_numbers = []
hdus = None
for file in listing:
    index1 = file.find("psfs_")
    index2 = file.find(".fits")
    if index1 != 0 or index2 < 0: continue
    number = file[index1+5:index2]
    if number.isdigit():
        psf_numbers.append(int(number))
print psf_numbers.sort()
if len(psf_numbers) < 2:
    print "Fewer than 2 psfs_[0-9]*.fits files found"
    sys.exit(1)
hdus = pyfits.open("psfs_%d.fits"%psf_numbers[0])
for psf_number in psf_numbers[1:]:
   hdus.append(pyfits.open("psfs_%d.fits"%psf_number)[0])
hdus.writeto(args.output_file, clobber=True)
print "%s file created with %d files"%(args.output_file, len(psf_numbers))

# delete the individual psf files if requested to
if args.delete:
    for psf_number in psf_numbers:
        os.unlink("psfs_%d.fits"%psf_number)

