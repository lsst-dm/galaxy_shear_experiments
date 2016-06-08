#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import glob
import unittest
import math
import numpy
import os
import pyfits
import time
import argparse

import lsst.utils.tests
import lsst.shapelet
import lsst.afw.geom
import lsst.afw.table
import lsst.afw.image
import lsst.afw.detection
import lsst.meas.modelfit
import lsst.meas.base
import lsst.meas.algorithms
import lsst.meas.extensions.ngmix

numpy.random.seed(500)


class ShapeCompare():
    #   vals is a list of numerical values.  Return the average and stdev
    #   if thresh is set, clip at 'thresh' sigma for 'iters' iteration
    def getStats(self, vals, thresh=None, iters=1):
        clip = None
        clipped = []
        for i in range(iters):
            if i == 0: clip = None
            else: clip = valAvg + thresh * valStdev
            valSum = 0.0
            valSS = 0.0
            valCount= 0
            for val in vals:
                if clip is None or val < clip:
                    valSum = valSum + val
                    valSS = valSS + val * val
                    valCount = valCount + 1
                else:
                    if i == (iters-1):
                        clipped.append(val)
            valAvg = valSum/valCount
            valStdev = math.sqrt((valCount * valSS - valSum * valSum)/ (valCount * (valCount - 1)))
            if thresh == None:
                 return valAvg, valStdev
        return valCount, valAvg, valStdev, clip, clipped

    #   Compare a psf with an array by using its computeImage
    #   Calculate the average difference and stdev for pixel differences.
    def comparePsfResults(self, psf, msf, weight=False):
        image = psf.computeImage()
        shape = image.getArray().shape
        newi = lsst.afw.image.ImageD(shape[1], shape[0])
        for comp in msf.getComponents():
            eval = comp.evaluate()
            #eval.addToImage(newi.getArray(), lsst.afw.geom.Point2I(-(shape[0]-1)/2, -(shape[1]-1)/2))
            eval.addToImage(newi.getArray(), lsst.afw.geom.Point2I(0, 0))
            core = eval.computeMoments().getCore()
        diffs = newi.getArray() - image.getArray()
        if weight:
            diffs *= image.getArray()
            diffs /= image.getArray().sum()
        return diffs.mean(), diffs.std()

    #   make a square gaussian array of size x size and width sigma
    def makeGaussianArray(self, size, sigma):
        xc = size/2
        yc = size/2
        image = lsst.afw.image.ImageD(lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(0,0),
                    lsst.afw.geom.Point2I(size,size)))
        array = numpy.ndarray(shape = (size, size), dtype = numpy.float64)
        for yi, yv in enumerate(xrange(0, size)):
            for xi, xv in enumerate(xrange(0, size)):
                array[yi, xi] = numpy.exp(-0.5*((xv - xc)**2 + (yv - yc)**2)/sigma**2)
        array /= array.sum()
        return array

    #   Create a config for a named psf model.  nGauss1,2,3 represent the EMPsfApprox.
    #   SingleGaussian, DoubleGaussian represent the ShapeletPsfApprox
    def getConfig(self, model):
        if model.find("nGauss") >= 0:
            return self.getEMConfig(model)
        config = lsst.meas.base.SingleFrameMeasurementTask.ConfigClass()
        config.slots.centroid = "base_GaussianCentroid"
        config.slots.shape = "base_SdssShape"
        config.slots.psfFlux = "base_PsfFlux"
        config.slots.apFlux = None
        config.slots.instFlux = None
        config.slots.modelFlux = None
        config.slots.calibFlux = None
        config.doReplaceWithNoise = False
        psfName="modelfit_GeneralShapeletPsfApprox"
        config.plugins.names = ["base_GaussianCentroid", "base_PsfFlux", "base_SdssShape",
           psfName, "modelfit_CModel"]
        config.plugins["modelfit_CModel"].psfName = psfName + "_%s"%(model,)
        config.plugins[psfName].sequence = [model]
        return config

    def getEMConfig(self, model):
        nGauss = int(model[6:])
        config = lsst.meas.base.SingleFrameMeasurementTask.ConfigClass()
        config.slots.centroid = "base_GaussianCentroid"
        config.slots.shape = None
        config.slots.psfFlux = None
        config.slots.apFlux = None
        config.slots.instFlux = None
        config.slots.modelFlux = None
        config.slots.calibFlux = None
        config.doReplaceWithNoise = False
        psfName="meas_extensions_ngmix_EmPsfApprox"
        config.plugins[psfName].nGauss = nGauss
        config.plugins[psfName].maxIter = 10000
        config.plugins[psfName].tolerance = 1e-6
        config.plugins.names = ["base_GaussianCentroid", psfName, "modelfit_CModel"]
        config.plugins["modelfit_CModel"].psfName = psfName
        return config

    def makePsf(self, params, size):
        if len(params) % 2 != 0:
            raise Exception("Must specify amplitude and sigma for each Gaussian")
        for i in range(len(params)/2):
            amp = float(params[2*i])
            sigma = float(params[2*i + 1])
            temp = self.makeGaussianArray(size, sigma)
            temp *= amp
            print "component %d: amp = %f, sigma = %f"%(i, amp, sigma)
            if i == 0:
                array = temp
            else:
                array += temp
        image = lsst.afw.image.ImageD(array)
        image.writeFits("psf.fits")
        kernel = lsst.afw.math.FixedKernel(image)
        return lsst.meas.algorithms.KernelPsf(kernel)


    #   Run a test of the specified model, collecting measurement information in
    #   a source table and returning it.  Also print information about the time
    #   and the quality of the PsfApprox
    def runShapeCompare(self, template, config, model, count=None, psf=None, verbose=True):
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(config=config, schema=schema)
        measCat = lsst.afw.table.SourceCatalog(schema)
        files = glob.glob("%s*"%template)
        if len(files) == 0:
            task.log.fatal("No files found matching '%s'"%template)
        timesCM = []
        timesSPA = []
        psfDiffs = []
        psfStds = []
        if len(files) > count:
            files = files[0:count]
        for file in files:
            exposure = lsst.afw.image.ExposureF(file)
            if not psf is None:
                exposure.setPsf(psf)
            dettask = lsst.meas.algorithms.SourceDetectionTask()
            dettask.log.setThreshold(dettask.log.WARN)
            task.log.setThreshold(task.log.WARN)
            footprints = dettask.detectFootprints(exposure, sigma=4.0).positive.getFootprints()
            measRecord = measCat.addNew()
            measRecord.setFootprint(footprints[0])
            psfName = config.plugins["modelfit_CModel"].psfName
            msfKey = lsst.shapelet.MultiShapeletFunctionKey(schema[psfName],
                                                             lsst.shapelet.HERMITE)
            quadKey = lsst.afw.table.QuadrupoleKey(schema["modelfit_CModel_exp_ellipse"])
            if verbose:
                print "Using model %s: on %s"%(model, file[file.find(template):])
            try:
                for plugin in task.plugins.keys():
                    startTime = time.time()
                    if plugin == "modelfit_CModel":
                        msf = measRecord.get(msfKey)
                        for comp in msf.getComponents():
                            if verbose:
                                print "PSF: ", comp.getEllipse(), comp.getCoefficients()
                    task.plugins[plugin].measure(measRecord, exposure)
                    if plugin.find("PsfApprox") > 0:
                        timeSPA = time.time() - startTime
                        diff = self.comparePsfResults(exposure.getPsf(), measRecord.get(msfKey))
                        psfDiffs.append(diff[0])
                        psfStds.append(diff[1])
                        if verbose:
                            print "DIFFS: ", diff
                        #iters = measRecord.get("meas_extensions_ngmix_EmPsfApprox_iterations")
                        #fdiff = measRecord.get("meas_extensions_ngmix_EmPsfApprox_fdiff")
                        #print "EMRunner for %d returned: "%nGauss, iters, fdiff
                    if plugin == "modelfit_CModel":
                        timeCM = time.time() - startTime
                ellipse = measRecord.get(quadKey)
                if verbose:
                    print "CMODEL: ", file[file.find("exp_"):], measRecord.getCentroid(), ellipse
                    print "SPA=", timeSPA, "CM=", timeCM
                timesCM.append(timeCM)
                timesSPA.append(timeSPA)
            except Exception as e:
                print plugin, "EXCEPTION", e.message
                continue
        if len(timesSPA) > 2:
            print "SUMMARY using PsfApprox %s with %d exposures"%(model, len(timesSPA))
            avg, stdev = self.getStats(timesSPA)
            print "   PsfApprox runtime: %.4f sec. +- %.6f sec."%(avg, stdev)
            avg, stdev = self.getStats(timesCM)
            print "   CModel runtime: %.4f sec. +- %.6f sec"%(avg, stdev)
            avg, stdev = self.getStats(psfDiffs)
            print "   PSF avg pixel diff: %.8f"%(avg,)
            avg, stdev = self.getStats(psfStds)
            print "   PSF std pixel dff stdev: %.8f"%(avg, )
        return measCat

#   compare the results of two difference runs of CModel by taking the ratios of the
#   moments returned by CModel_dev_ellipse
def compareResults(cat1, cat2):
    val1Count = 0
    val1Sum = 0.
    val1SS = 0;
    val2Count = 0
    val2Sum = 0.
    val2SS = 0;
    key1 = lsst.afw.table.QuadrupoleKey(cat1.getSchema()["modelfit_CModel_dev_ellipse"])
    key2 = lsst.afw.table.QuadrupoleKey(cat2.getSchema()["modelfit_CModel_dev_ellipse"])
    for i in range(len(cat1)):
        xx = cat1[i].get(key1).getIxx()
        yy = cat1[i].get(key1).getIyy()
        xy = cat1[i].get(key1).getIxy()
        e1_1 = (xx - yy) / (xx + yy)
        e2_1 = 2 * xy / (xx + yy)
        xx = cat2[i].get(key2).getIxx()
        yy = cat2[i].get(key2).getIyy()
        xy = cat2[i].get(key2).getIxy()
        e1_2 = (xx - yy) / (xx + yy)
        e2_2 = 2 * xy / (xx + yy)
        val1 = (e1_1 - e1_2)
        val1Sum += val1
        val1SS += val1 * val1
        val1Count += 1
        val2 = (e2_1 - e2_2)
        val2Sum += val2
        val2SS += val2 * val2
        val2Count += 1
    val1Avg = val1Sum/val1Count
    val1Stdev = math.sqrt((val1Count * val1SS - val1Sum * val1Sum)/ (val1Count * (val1Count - 1)))
    print "e1 differences: %.8f +- %.8f for %d sources"%(val1Avg, val1Stdev, val1Count)
    val2Avg = val2Sum/val2Count
    val2Stdev = math.sqrt((val2Count * val2SS - val2Sum * val2Sum)/ (val2Count * (val2Count - 1)))
    print "e2 differences: %.8f +- %.8f for %d sources"%(val2Avg, val2Stdev, val2Count)

if __name__ == "__main__":
    dataDir = os.path.join(os.environ["GALAXY_SHEAR_EXPERIMENTS_DIR"], "examples/data")
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--template", type=str, help="select files", default="exp_")
    parser.add_argument("-n", "--number", type=int, help="maximum number to do", default=10)
    parser.add_argument("-v", "--verbose", action='store_true', help="print verbose statistics")
    parser.add_argument("-m", "--models", type=str, help="SingleGaussian,...,nGauss#", default=None)
    parser.add_argument("-c", "--compare", action='store_true', help="compare CModels")
    parser.add_argument("-p", "--psf", type=str, help="specify a PSF with one or more Gaussians", default=None)
    args = parser.parse_args()
    test = ShapeCompare()
    psf = None
    if not args.psf is None:
        psf = test.makePsf(args.psf.split(","), 67)
    else:
        psf = None

    if not args.models is None:
        models = args.models.split(",")
        results = {}
        for model in models:
           config = test.getConfig(model)
           results[model] = test.runShapeCompare(os.path.join(dataDir, args.template), config, model, count=args.number, psf=psf, verbose=args.verbose)
        if args.compare:
            print "CModel comparison: for %s vs %s:"%(models[0], models[1])
            results2 = compareResults(results[models[0]], results[models[1]])
