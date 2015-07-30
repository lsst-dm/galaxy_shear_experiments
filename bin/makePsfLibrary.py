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
makePsfLibrary is a python program which can be called as a command line script.
It is used to create libraries of PhoSim focal planes which cover the whole
LSST focal plane, divided into groups (currently by raft)
"""

#   This routine takes a list of images and a catalog name (presumably covering the
#   same area as the images.  Any entry in the catalog which is fully contained in an
#   image is cutout, transformed to the center of the pupil, and then placed in a libray
#   (multi-extension fits file).  The entry is also logged to an open file handle, fout
def makeSubLibrary(images, catname, output_file, fout, psfSize=33, satValue=None):
    # Create an XYTransform from TanPixels to Pixels.  This is the transform
    # needed to remove the optical distortion from the Psf. Use R22_S11 as center
    butler = dafPersist.Butler(root=os.path.join(os.path.dirname(__file__), "."))
    camera = butler.get("camera", immediate=True)
    det = camera["R:2,2 S:1,1"]
    tm = det.getTransformMap()
    transform1 = tm.get(det.makeCameraSys(cameraGeom.PIXELS)).invert()
    transform2 = tm.get(det.makeCameraSys(cameraGeom.TAN_PIXELS))
    xytransform = afwGeom.MultiXYTransform([transform1, transform2])

    # A centroid algorithm has to be run here, because the catalog from OpSim
    # is not very accurate
    schema = afwTable.SourceTable.makeMinimalSchema()
    control = measBase.GaussianCentroidControl()
    alg = measBase.GaussianCentroidAlgorithm
    centroidPlugin = alg(control, "centroid", schema)
    alg = measBase.SdssShapeAlgorithm
    control = measBase.SdssShapeControl()
    shapePlugin = alg(control, "shape", schema)

    #   Go through all the images and write out the galaxies psfNo is the line in the input
    #   psfNo is the number of the object in the input file (0 based)
    #   psfIndex is the index of the psf's hdu in its group (0 based)
    psfIndex = 0
    outputHdus = pyfits.hdu.hdulist.HDUList()
    starList = []
    for image in images:
        identifier = image[image.find("_R")+1:image.find("_R")+8]
        exp = afwImage.ExposureF(image)
        cat = open(catname, "r")
        psfNo = -1
        # number the psfs the same way they are listed in the catalog, starting with line 1
        # read through the OpSim catalog and create a psf for any object which is completely
        # inside the image we are processing
        for line in cat:
            # catalog lines with objects have x,y position at index 2,3.  Use that for 1st guess
            attrs = line.split()
            if attrs[0] != "object": continue
            psfNo = psfNo + 1
            position = afwGeom.Point2D(float(attrs[2]), float(attrs[3]))
            coord = afwCoord.Coord(position)
            pixels = exp.getWcs().skyToPixel(coord)
            x = int(pixels.getX())
            y = int(pixels.getY())
            starInfo = {"image":os.path.basename(image), "psfNo":psfNo, "psfIndex":-1, "pixels":pixels}
            starBBox = afwGeom.Box2I(afwGeom.Point2I(x, y), afwGeom.Extent2I(1,1))
            starBBox.grow(psfSize)

            # test to be sure that the star is completely contained
            if exp.getBBox().contains(starBBox):
                subImage = afwImage.ImageF(exp.getMaskedImage().getImage(), starBBox)
                #   high SNR ratio objects with over 50,000 counts.  Crude detection
                #   required allows a large detection threshold of 100
                ds = afwDetection.FootprintSet(subImage, afwDetection.Threshold(100))
                #   if the detection is clean, there should be one footprint with
                #   at least one peak.
                footprintCnt = len(ds.getFootprints())
                starInfo["footprintCnt"] = footprintCnt
                if footprintCnt == 0 or footprintCnt > 1:
                    starInfo["footprints"] = footprintCnt
                if len(ds.getFootprints()[0].getPeaks()) == 0:
                    starInfo["peaks"] = True
                #   Check the peak value if satValue is set
                peakValue = ds.getFootprints()[0].getPeaks()[0].getPeakValue()
                if not satValue is None and peakValue >= satValue:
                    starInfo["peakValue"] = peakValue
                sources = afwTable.SourceCatalog(schema)
                sources.defineCentroid("centroid")
                sources.defineShape("shape")
                record = sources.addNew()
                record.setFootprint(ds.getFootprints()[0])
                try:
                    #   call centroid and shape to test the ellipticity
                    #   this is a crude test at .2, just intended to cull outliers
                    centroidPlugin.measure(record, exp)
                    shapePlugin.measure(record, exp)
                    shape = record.getShape()
                    dims = afwGeom.ellipses.Ellipse(shape).getCore().computeDimensions()
                    ellipticity = abs((dims[0]-dims[1])/dims[0])
                    if ellipticity > .2:
                        starInfo["ellipticity"] = ellipticity
                except:
                    starInfo["ellipticity"] = "NaN"

                # shift the bounding box to the correct centroid
                newCoord = exp.getWcs().pixelToSky(afwGeom.Point2D(record.get("centroid_x"), record.get("centroid_y")))
                dx = int(round(record.get("centroid_x"))) - x
                dy = int(round(record.get("centroid_y"))) - y
                starInfo["centroid"] = record.getCentroid()
                starInfo["shift_x"] = dx
                starInfo["shift_y"] = dy
                subImage = afwImage.ImageF(exp.getMaskedImage().getImage(), starBBox)
                starBBox.shift(afwGeom.Extent2I(dx, dy))
                # if the exposure still contains the cutout, use it.
                if exp.getBBox().contains(starBBox):
                    subImage = afwImage.ImageF(exp.getMaskedImage().getImage(), starBBox)
                    # convert the image to 64 bits, and turn it into a Psf
                    # then warp it using the xytransform from above
                    data = subImage.getArray().astype(numpy.float64)
                    kernel = afwMath.FixedKernel(afwImage.ImageD(data))
                    psf = measAlg.KernelPsf(kernel)
                    warpedPsf = measAlg.WarpedPsf(psf, xytransform, "lanczos3");
                    if len(outputHdus) == 0:
                        hdu = pyfits.hdu.image.PrimaryHDU()
                    else:
                        hdu = pyfits.hdu.image.ImageHDU()
                    try:
                        hdu.data = warpedPsf.computeImage().getArray()
                    except:
                        starInfo["computeImage"] = True
                        starList.append(starInfo)
                        continue
                    #  pixel scale of .2 arcsecs/pixel added for GalSim
                    hdu.header.append("GS_SCALE")
                    hdu.header.set("GS_SCALE", .2)
                    hdu.header.append("PSF_NO")
                    hdu.header.set("PSF_NO", psfNo)
                    #  psf number in input file added in case it is needed downstream
                    outputHdus.append(hdu)
                    starInfo["psfIndex"] = psfIndex
                    starInfo["bbox"] = starBBox
                    psfIndex = psfIndex + 1
                #   otherwise, mark as an edge case
                else:
                    starInfo["edge"] = True
                starList.append(starInfo)
        cat.close()
        del hdu
        del exp
    #   Now go through all the starInfos and delete any the might intersect
    for info in starList:
        for info2 in starList:
                if not info["psfNo"] == info2["psfNo"] and info["image"] == info2["image"]:
                    d = math.sqrt(info["pixels"].distanceSquared(info2["pixels"]))
                    if d < 50:
                        if not "overlap" in info.keys() or info["overlap"] > d:
                            info["overlap"] = str(d)

    outputHdus.writeto(output_file, clobber=True)
    print len(outputHdus), " psf images were cutout from the original images and written to", output_file
    return starList

#   This main program creates a library for each raft in the input_directory.  A raft is assumed
#   to be comprised of 9 files named basename_Rnn_Smm.fits in the input directory
#   The cutouts from each set of 9 files is placed in the output_dir in basename_Rnn.fits
#   An index of all the libraries created in this way is created in basename.index
if __name__ == "__main__":
    """
    Script to create a library of Psfs from an lsst focal plane.
    Command line parameters are:

    catalog          A par file from PhoSim which lists all the stars
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
    parser.add_argument("-o", "--output_dir", help="Name of the output directory", default=None)
    parser.add_argument("-p", "--psf_size", help="Size of psf images to be cutout", type = int, default=33)
    parser.add_argument("-s", "--sat_value", help="saturation value for input images", type = int, default=None)
    parser.add_argument("-i", "--include", help="Use only images containing this string", default=None)
    args = parser.parse_args()
    #   Collect a list of available rafts
    #   The current convention is the output from phosim, which sends us files which
    #   are tagged with Rnn_Smm and are .fits.gz files.  We may need to change this
    #   for HSC files
    rafts = []
    for file in os.listdir(args.input_dir):
        if file.find(args.basename) >= 0 and file.find(".fits.gz") >= 0:
            raftIndex = file.find("_R")
            if raftIndex > 0:
                entry = file[raftIndex+1 : raftIndex+4]
                if not entry in rafts:
                    rafts.append(entry)
    print len(rafts), "rafts found"
    print rafts
    #   Now process the Psfs raft-by-raft
    baseDir = args.output_dir
    if baseDir is None:
        baseDir = os.path.join(args.input_dir, "psfs")

    if not os.path.isdir(baseDir):
        os.mkdir(baseDir)
    fout = open(os.path.join(baseDir, "psfs.index"), "w")
    ferr = open(os.path.join(baseDir, "psfs.errors"), "w")
    for raft in rafts:
        images = []
        #   From phosim, the files are typically names base_Rnn_Smm.fits.gz
        #   where base is an identifier name we supply from our PhoSim runs
        for file in os.listdir(args.input_dir):
            if file.find(args.basename + "_" + raft) >= 0:
                if args.include == None or len(args.include) == 0 or file.find(args.include) >= 0:
                    imagename = os.path.join(args.input_dir, file)
                    if os.path.isfile(imagename) and not imagename in images:
                        images.append(imagename)
        catname = os.path.join(args.input_dir, args.catalog)
        output_file = os.path.join(baseDir, "psf_library_" + raft[1:3] + ".fits")
        #   If there are any images for this raft, cutout the Psfs and place in a mini-library
        if len(images) > 0:
            starInfo = makeSubLibrary(images, catname, output_file, fout, psfSize=args.psf_size, satValue=args.sat_value)
            #   if none of the errors are indicated, it is Ok to write this one to the output file
            #   dump all the errors to the error file
            for info in starInfo:
                nErrors = 0
                for flag in ("overlap", "footprints", "peakValue", "peaks", "edge", "ellipticity", "computeImage"):
                    if flag in info.keys():
                        if nErrors == 0:
                            ferr.write("psfNo: %d"%info["psfNo"])
                        ferr.write(" ")
                        ferr.write(" %s=%s"%(flag, str(info[flag])))
                        nErrors = nErrors + 1
                if nErrors > 0:
                    ferr.write("\n")
                    ferr.flush()
                else:
                    print "library: ", raft, info["psfNo"]
                    fout.write("%s %s %d %d (%d,%d)"%(os.path.basename(output_file), raft,
                               info["psfNo"], info["psfIndex"],
                                info["pixels"].getX(), info["pixels"].getY()))
                    fout.write("\n")
        else:
            print raft, " no images coincided" 
    fout.close()
    ferr.close()
