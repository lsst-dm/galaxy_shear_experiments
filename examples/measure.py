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

def getIndex(name, cat):
    index = 0
    for item in cat.columns:
        if name == item.name:
            return index
        index = index + 1
    return -1

parser = argparse.ArgumentParser()
parser.add_argument("type", help="Name of galaxy type (control or real_galaxy")
parser.add_argument("psfs", help="location of psfs file or psfs dir")
args = parser.parse_args()
hdus = pyfits.open("%s/ground/constant/epoch_catalog-000-0.fits"%args.type)
cat = hdus[1]
exp = afwImage.ExposureF("%s/ground/constant/image-000-0.fits"%args.type)

# The psfs may be stored in a directory with names psfs_nnnn.fits
# or in a multiple hdu fits files.  If it is a multiple hdu fits file,
# the header of each psf should contain a "PSF_NO" to indicate the original
# number of the psf before packaging.  If "PSF_NO" is not in the header,
# the Psfs are assumed to be numbered 1, N where N is the number of psfs
psf_dict = None
if not os.path.isdir(args.psfs):
    if not os.path.exists(args.psfs):
        print "Psf location: %s does not exist."%args.psfs
        sys.exit(1)
    psf_hdus = pyfits.open(args.psfs)
    psf_dict = {}
    for i in range(len(psf_hdus)):
        psf_number = psf_hdus[i].header.get("PSF_NO")
        if psf_number == None:
            psf_dict[i+1] = i
        else:
            psf_dict[psf_number] = i
else:
    psfFormat = args.psfs + "/psfs_%d.fits"

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
    bbox = afwGeom.Box2I(afwGeom.Point2I(item[i_xmin], item[i_ymin]), afwGeom.Point2I(item[i_xmax], item[i_ymax]))
    subexp = afwImage.ExposureF(exp, bbox)
    #subarray = subexp.getMaskedImage().getImage().getArray()
    #subarray = numpy.repeat(numpy.repeat(subarray,2, axis=0), 2, axis=1) 
    #subexp = afwImage.ExposureF(afwImage.MaskedImageF(afwImage.ImageF(subarray)))
    if psf_dict:
        psf_number = item[i_psf_number]
        i = psf_dict[psf_number]
        data = psf_hdus[i].data.astype(numpy.float64) 
    else:
        psf_file = psfFormat % item[i_psf_number]
        psf_image = afwImage.ImageF(psf_file)
        data = psf_image.getArray().astype(numpy.float64)
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
    except lsst.pex.exceptions.wrappers.RuntimeError as e:
        if e.message.find("PSF") > 0:
            subexp.writeFits("psf_exp_%d.fits"%item[i_index])
            subexp.getPsf().computeImage(cen).writeFits("psf_psf_%d.fits"%item[i_index])
        if e.message.find("iteration") > 0:
            subexp.writeFits("iter_exp_%d.fits"%item[i_index])
            subexp.getPsf().computeImage(cen).writeFits("iter_psf_%d.fits"%item[i_index])
        if e.message.find("all 0") > 0:
            import pdb
            pdb.set_trace()
            subexp.writeFits("zero_exp_%d.fits"%item[i_index])
            subexp.getPsf().computeImage(cen).writeFits("zero_psf_%d.fits"%item[i_index])
 
        print item[i_index], e.message
        failed = failed + 1
    except BaseException as e:
        print "FAILED: ", item[i_index], e.message
        failed = failed + 1
print "elapsed time = ", time.time()-startTime
print "success = ", success, ", failed = ", failed
