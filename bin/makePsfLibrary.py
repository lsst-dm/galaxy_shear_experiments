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
import pyfits
import argparse
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

psfSize = 33
rafts = ["R01","R02","R03","R10","R11","R12","R13","R14","R20","R21","R22","R23","R24","R30","R31","R32","R33","R34","R41","R42","R43"]

"""
    Script to create a library of Psfs from an lsst focal plane.

    @param[in]  catalog          A par file from PhoSim which lists all the stars
    @param[in]  input_dir        Directory containing the fits.gz files from PhoSim
    @param[in]  basename         prefix for the fits.gz files, excluding Rnn_Smm.fits.gz
    @param[in]  output_dir       Directory in which to place the output libraries.

    The output index is placed in output_dir/basename.index
    The library is collected in fits files named output_dir/basename_Rnn.psfs.fits
"""

#   This routine takes a list of images and a catalog name (presumably covering the
#   same area as the images.  Any entry in the catalog which is fully contained in an
#   image is cutout, transformed to the center of the pupil, and then placed in a libray
#   (multi-extension fits file).  The entry is also logged to an open file handle, fout
def makeSubLibrary(images, catname, output_file, fout):
    # Create an XYTransform from TanPixels to Pixels.  This is the transform
    # needed to remove the optical distortion from the Psf. Use R22_S11 as center
    butler = dafPersist.Butler(root=os.path.join(os.path.dirname(__file__), "."))
    camera = butler.get("camera", immediate=True)
    det = camera["R:2,2 S:1,1"]
    tm = det.getTransformMap()
    for sys in tm.getCoordSysList():
        if sys.getSysName() == "Pixels":
            transform1 = tm.get(sys).invert()
        if sys.getSysName() == "TanPixels":
            transform2 = tm.get(sys)
    xytransform = afwGeom.MultiXYTransform([transform1, transform2])
    
    # A centroid algorithm has to be run here, because the catalog from OpSim
    # is not very accurate
    schema = afwTable.SourceTable.makeMinimalSchema()
    control = measBase.GaussianCentroidControl()
    alg = measBase.GaussianCentroidAlgorithm
    plugin = alg(control, "", schema)
    
    # count of psfs actually created
    psfIndex = 1
    # Go through all the images and write out the galaxies
    psfCount = 0
    outputHdus = pyfits.hdu.hdulist.HDUList()
    for image in images:
        exp = afwImage.ExposureF(image)
        cat = open(catname, "r")
        # number the psfs the same way they are listed in the catalog, starting with line 1
        # read through the OpSim catalog and create a psf for any object which is completely
        # inside the image we are processing
        while True:
            line = cat.readline()
            if len(line) == 0:
                break
            # catalog lines with objects have x,y position at index 2,3.  Use that for 1st guess
            attrs = line.split()
            if attrs[0] != "object": continue
            coord = afwCoord.Coord(afwGeom.Point2D(float(attrs[2]), float(attrs[3])))
            pixels = exp.getWcs().skyToPixel(coord)
            x = int(pixels.getX())
            y = int(pixels.getY())
            starBBox = afwGeom.Box2I(afwGeom.Point2I(x, y), afwGeom.Extent2I(1,1))
            starBBox.grow(psfSize)
        
            # test to be sure that the star is completely contained
            if exp.getBBox().contains(starBBox):
                subImage = afwImage.ImageF(exp.getMaskedImage().getImage(), starBBox)
                ds = afwDetection.FootprintSet(subImage, afwDetection.Threshold(1000))
                sources = afwTable.SourceCatalog(schema)
                record = sources.addNew()
                record.setFootprint(ds.getFootprints()[0])
                plugin.measure(record, exp)
        
                # shift the bounding box to the correct centroid
                newCoord = exp.getWcs().pixelToSky(afwGeom.Point2D(record.get("_x"), record.get("_y")))
                dx = int(round(record.get("_x"))) - x
                dy = int(round(record.get("_y"))) - y
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
                    hdu.data = warpedPsf.computeImage().getArray()
                    hdu.header.append("GS_SCALE")
                    hdu.header.set("GS_SCALE", .2)
                    hdu.header.append("PSF_NO")
                    hdu.header.set("PSF_NO",psfIndex)
                    
                    outputHdus.append(hdu)
                    fout.write("%s %d %d\n"%(os.path.basename(output_file), psfCount, psfIndex))
                    psfCount = psfCount + 1
            psfIndex = psfIndex + 1
        cat.close()
        del hdu
        del exp
    outputHdus.writeto(output_file, clobber=True)
    print psfCount, " psf images were cutout from the original image", len(outputHdus)

#   This main program creates a library for each raft in the input_directory.  A raft is assumed
#   to be comprised of 9 files named basename_Rnn_Smm.fits in the input directory
#   The cutouts from each set of 9 files is placed in the output_dir in basename_Rnn.fits
#   An index of all the libraries created in this way is created in basename.index
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("catalog", help="Name of the star catalog fed to OpSim")
    parser.add_argument("input_dir", help="Name of the star grid produced by OpSim")
    parser.add_argument("basename", help="Name of the star grid produced by OpSim")
    parser.add_argument("output_dir", help="Name of the output directory")
    args = parser.parse_args()
    # List of available rafts for each focal plane
    fout = open(args.output_dir + "/" + args.basename + ".index", "w")
    for raft in rafts:
        images = []
        for file in os.listdir(args.input_dir):
            if file.find(args.basename + "_" + raft + "_S") >= 0:
                imagename = args.input_dir + "/" + file
                if os.path.isfile(imagename):
                    images.append(imagename)
        catname = args.input_dir + "/" + args.catalog
        output_file = args.output_dir + "/" + args.basename + "." + raft + ".psfs.fits"
        #   If there are any images for this raft, cutout the Psfs and place in a mini-library
        if len(images) > 0:
            makeSubLibrary(images, args.input_dir + "/" + args.catalog, output_file, fout)
        del images 
    fout.close()
