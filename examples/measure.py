#!/usr/bin/env python
import great3sims
import argparse
import numpy
import pyfits
import sys
import os
import time
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.meas.algorithms
import lsst.meas.base as measBase
import lsst.meas.extensions.shapeHSM as shape
from lsst.meas.algorithms import SourceDetectionTask

"""
    Simple test script to use the results from the GalSim Galaxy creation
    and the PhoSim Psf making scripts, and show that we can run measurements
    on the results.  Script was also used to report timings.
"""

def getIndex(name, cat):
    index = 0
    for item in cat.columns:
        if name == item.name:
            return index
        index = index + 1
    return -1

# The psfs may be stored in a directory with names psfs_nnnn.fits
# or in a multiple hdu fits files.  If it is a multiple hdu fits file,
# the header of each psf should contain a "PSF_NO" to indicate the original
# number of the psf before packaging.  If "PSF_NO" is not in the header,
# the Psfs are assumed to be numbered 1, N where N is the number of psfs
def main(args):
    hdus = pyfits.open("%s/ground/constant/epoch_catalog-000-0.fits"%args.type)
    cat = hdus[1]
    exp = afwImage.ExposureF("%s/ground/constant/image-000-0.fits"%args.type)
    psfDict = None
    if not os.path.isdir(args.psfs):
        if not os.path.exists(args.psfs):
            print "Psf location: %s does not exist."%args.psfs
            sys.exit(1)
        psfHdfs = pyfits.open(args.psfs)
        psfDict = {}
        for i in range(len(psfHdfs)):
            psf_number = psfHdfs[i].header.get("PSF_NO")
            if psf_number == None:
                psfDict[i+1] = i
            else:
                psfDict[psf_number] = i
    else:
        psfFormat = args.psfs + "/psf_%d.fits"
    
    i_psf_number = getIndex("psf_number", cat)
    i_index = getIndex("index", cat)
    i_xmax = getIndex("xmax", cat)
    i_ymax = getIndex("ymax", cat)
    i_xmin = getIndex("xmin", cat)
    i_ymin = getIndex("ymin", cat)
    
    schema = afwTable.SourceTable.makeMinimalSchema()
    schema.addField("centroid_x", type=float)
    schema.addField("centroid_y", type=float)
    control = shape.HsmShapeBjControl()
    alg = shape.HsmShapeBjAlgorithm
    #control = measBase.SdssShapeControl()
    #alg = measBase.SdssShapeAlgorithm
    plugin = alg(control, "", schema)
    success = 0
    failed = 0
    count = 0
    startTime = time.time()
    for item in cat.data:
        bbox = afwGeom.Box2I(
                        afwGeom.Point2I(item[i_xmin], item[i_ymin]),
                        afwGeom.Point2I(item[i_xmax], item[i_ymax])
        )
        subexp = afwImage.ExposureF(exp, bbox)
        #subarray = subexp.getMaskedImage().getImage().getArray()
        #subarray = numpy.repeat(numpy.repeat(subarray,2, axis=0), 2, axis=1) 
        #subexp = afwImage.ExposureF(afwImage.MaskedImageF(afwImage.ImageF(subarray)))
        if psfDict:
            psf_number = item[i_psf_number]
            i = psfDict[psf_number]
            data = psfHdfs[i].data.astype(numpy.float64) 
        else:
            psfFile = psfFormat % item[i_psf_number]
            psfImage = afwImage.ImageF(psfFile)
            data = psfImage.getArray().astype(numpy.float64)
        kernel = afwMath.FixedKernel(afwImage.ImageD(data))
        psf = lsst.meas.algorithms.KernelPsf(kernel)
        subexp.setPsf(psf)
        subexp.setXY0(afwGeom.Point2I(0,0))
        detection = SourceDetectionTask()
        table = afwTable.SourceTable.make(schema)
        detections = detection.run(table, subexp)
        detections.sources.defineCentroid("centroid")
        count = count+1
        if len(detections.sources) == 0:
            print "ITEM %d, no sources, "%count, bbox
            continue
        record = detections.sources[0]
        
        record.setFootprint(detections.fpSets.positive.getFootprints()[0])
        cen = record.getFootprint().getPeaks()[0].getCentroid()
        record.set('centroid_x', cen.getX())
        record.set('centroid_y', cen.getY())
        try:
            plugin.measure(record, subexp)
            print item[i_index], record.get("_e1"), record.get("_e2")
            success = success+1
        except lsst.pex.exceptions.RuntimeError as e:
            print item[i_index], e.message
            failed = failed + 1
        except Exception as e:
            print "FAILED: ", item[i_index], e.message
            failed = failed + 1
    print "elapsed time = ", time.time()-startTime
    print "success = ", success, ", failed = ", failed

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("type", help="Name of galaxy type (control or real_galaxy")
    parser.add_argument("psfs", help="location of psfs file or psfs dir")
    args = parser.parse_args()
    main(args) 
