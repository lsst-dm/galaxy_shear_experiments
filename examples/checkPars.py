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
import math
import sys
import os.path
import numpy
import pyfits

import lsst.daf.persistence as dafPersist
import lsst.afw.coord as afwCoord
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg
import lsst.meas.base as measBase

"""
checkPar.py is a python program which can be called as a command line script.
It is used to test the distribution of a library over the entire
LSST focal plane, divided into groups (currently by raft)
"""

#   This routine takes a list of images and a catalog name (presumably covering the
#   same area as the images.  Any entry in the catalog whose center is contained in
#   the bounding rectangle of a given image is counted as belonging both to the
#   raft and the sensor of the image
def getRaftStatistics(images, catname, raft):
    imageinfo = []
    #   Go through all the images and see how many sources each contains
    for i, image in enumerate(images):
        exp = afwImage.ExposureF(image)
        width, height = exp.getMaskedImage().getImage().getArray().shape
        llcoord = exp.getWcs().pixelToSky(afwGeom.Point2D(0,0))
        urcoord = exp.getWcs().pixelToSky(afwGeom.Point2D(width-1, height-1))
        ramin = llcoord.toFk5().getRa().asDegrees()
        ramax = urcoord.toFk5().getRa().asDegrees()
        if ramin > ramax:
            temp = ramin
            ramin = ramax
            ramax = temp
        decmin = llcoord.toFk5().getDec().asDegrees()
        decmax = urcoord.toFk5().getDec().asDegrees()
        if decmin > decmax:
            temp = decmin
            decmin = decmax
            decmax = temp
        ind = image.find("_R")
        raft = image[ind+1:ind+5]
        ind = image.find("_S")
        sensor = image[ind+1:ind+5]
        imageinfo.append((raft, sensor, ramin, ramax, decmin, decmax))
    #   hist is a counting array, one entry per image
    hist = numpy.zeros(len(images), dtype = int)
    # read through the OpSim catalog and count which raft, sensor it is in
    # inside the image we are processing
    cat = open(catname, "r")
    for line in cat:
        # catalog lines with objects have x,y position at index 2,3.  Use that for 1st guess
        attrs = line.split()
        if attrs[0] != "object": continue
        ra = float(attrs[2])
        dec = float(attrs[3])
        for i, info in enumerate(imageinfo):
            ramin = info[2]
            ramax = info[3]
            decmin = info[4]
            decmax = info[5]
            if ra <= ramax and ra >= ramin and dec <= decmax and dec >= decmin:
                hist[i] = hist[i] + 1
    for i, info in enumerate(imageinfo):
        print "    ", hist[i], imageinfo[i]
    print hist[0], hist.sum()
    return hist.sum()

#   This main program does counts from images created by PhoSim, over rafts and sensors
#   Each raft is comprised of 9 files named basename_Rnn_Smm.fits in the input directory
if __name__ == "__main__":
    """
    Script to create a library of Psfs from an lsst focal plane.
    Command line parameters are:

    catalog          A pars file from PhoSim which lists all the stars
    input_dir        Directory containing the fits.gz files from PhoSim
    basename         prefix for the fits.gz files, excluding Rnn_Smm.fits.gz
    output_dir       Directory in which to place the output libraries.

    The output index is placed in output_dir/psfs.index
    The library is written to output_dir/basename/psfs/psf_library_nn.fits
    where nn is a two digit form of the raft number
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("catalog", help="Name of the star catalog fed to PhoSim")
    parser.add_argument("input_dir", help="Name of the input dir containing the images")

    parser.add_argument("basename", help="The opthistid provided to PhoSim")
    parser.add_argument("-f", "--filter", help="Use only images containing this string", default=2)
    args = parser.parse_args()
    #   Collect a list of available rafts
    #   The current convention is the output from phosim, which sends us files which
    #   are tagged with Rnn_Smm and are .fits.gz files.  We may need to change this
    #   for HSC files
    rafts = []
    input_dir = os.path.join(args.input_dir, args.basename)
    for file in os.listdir(input_dir):
        if file.find(args.basename) >= 0 and file.find(".fits.gz") >= 0:
            raftIndex = file.find("_R")
            if raftIndex > 0:
                entry = file[raftIndex+1 : raftIndex+4]
                if not entry in rafts:
                    rafts.append(entry)
    total = 0
    #   Now process the images raft-by-raft
    for raft in rafts:
        images = []
        #   From phosim, the files are typically names basename_Rnn_Smm.fits.gz
        #   where base is an identifier name we supply from our PhoSim runs
        for file in os.listdir(input_dir):
            if file.find(args.basename + "_" + "f" + str(args.filter)
                         + "_" + raft) >= 0 and file.find(".fits.gz") > 0:
                imagename = os.path.join(input_dir, file)
                if os.path.isfile(imagename) and not imagename in images:
                    images.append(imagename)
        #   The catalog is either in the input.basename directory, or the current directory
        catname = os.path.join(input_dir, args.catalog)
        if not os.path.isfile(catname):
            catname = os.path.join(".", args.catalog)
        if not os.path.isfile(catname):
            print catname, " also not found.  exiting..."
            sys.exit(1)
        if len(images) > 0:
            total = total + getRaftStatistics(images, catname, raft)
    print "Total = ", total
