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

"""
analyzeShearTest is a portion of the shear test suite which summarizes the information
from all of the GalSim runs and creates a table which is suitable for plotting.

The data is divided into GalSim runs each with a filter[23], seeing[.05 .07 .09] and shear.
All subfields in the run have the same shear, but each has a randomly chosen shear angle.
In this program, we average each subfield to get an estimator of the shear value.  In our
initial tests, there are 32 subfields with a different shear angle.  There are 1024 galaxies
for each subfield, comprised of 512 galaxy pair rotated at 90 degrees to each other.

These 32 samples are averaged to provide an overall estimator of the shear, and the stddev of
these samples provides an error for this estimator.
"""

import argparse
import math
import os.path

import lsst.daf.persistence as dafPersist
import lsst.afw.table as afwTable

from lsst.galaxy_shear.shearConfig import RunShearConfig
def getStandardAngle(angle):
    while angle > math.pi*2.0:
        angle = angle - math.pi*2.0
    while angle < 0:
        angle = angle + math.pi*2.0
    return angle * 180.0/math.pi

def runAnal(baseDir, outFile, config, test=None):

    #  Create a summary table of the src*.fits tables in baseDir
    schema = afwTable.Schema()
    filterKey = schema.addField("filter", type = int, doc = "filter 2 or 3 used in PhoSim psf generator.")
    seeingKey = schema.addField("seeing", type = float, doc = "raw seeing used by PhoSim.")
    shearValueKey = schema.addField("shear_value", type = float, doc = "constant shear level used in the GalSim run.")
    countKey = schema.addField("nsource", type = int, doc = "total number of galaxies in this run")
    gKey = schema.addField("g", type = float, doc = "g value calcuted from GalSim g1 and g2")
    g1Key = schema.addField("g1", type = float, doc = "g value calcuted from GalSim g1 and g2")
    g2Key = schema.addField("g2", type = float, doc = "g value calcuted from GalSim g1 and g2")
    eAvgKey = schema.addField("eAvg", type = float, doc = "average e for the sources in this subfield")
    eStdKey = schema.addField("eStd", type = float, doc = "stddev e for the sources in this subfield")
    e1AvgKey = schema.addField("e1Avg", type = float, doc = "e1 deviation from g1, for nsub subfields.")
    e1StdKey = schema.addField("e1Std", type = float, doc = "stddev e for the sources in this subfield")

    e2AvgKey = schema.addField("e2Avg", type = float, doc = "e2 deviation from g2, for nsub subfields")
    e2StdKey = schema.addField("e2Std", type = float, doc = "stddev e for the sources in this subfield")
    outCat = afwTable.BaseCatalog(schema)

    #   The shape type can be "moments" or "ellipticity"
    #   If we done have an e_sigma, the signal field  tells us to extract the weighting function
    #   from the SNR of some other field, where the sigma will be the error, and the
    #   The sigma used for weighting is shape_field/(sigma_field/sigma_signal_field)
    if "shape_field" in config.keys():
        shape_field = config.shape_field + "_"
        shape_type = config.shape_type
        if "sigma_field" in config.keys():
            sigma_field = config.sigma_field
            if "sigma_signal_field" in config.keys():
                sigma_signal_field = config.sigma_signal_field
            else:
                sigma_signal_field = None
        else:
            sigma_field = None
            sigma_signal_field = None
    else:
        #   If shape_field is not in params, the ellipticity and sigma field are here:
        shape_field = ""
        shape_type = "ellipticity"
        sigma_field = "e_sigma"
        sigma_signal_field = None

    #   Use the butler to get the src.fits file for each subfield/epoch
    #butler = dafPersist.Butler(root=os.path.join(baseDir, config.exp_type + "/ground/constant"))

    count = 0
    nanCount = 0
    allCount = 0
    eAvgSum = 0.0
    e1AvgSum = 0.0
    e2AvgSum = 0.0
    eAvgSumSq = 0.0
    e1AvgSumSq = 0.0
    e2AvgSumSq = 0.0
    weightSum = 0.0
    nAvgs = 0
    flagKeys = []
    flagNames = []
    flagCount = []
    for subfield in range(config.n_subfields):
        for epoch in range(config.n_epochs):
            catCount = 0
            catESum = 0.0
            catE1Sum = 0.0
            catE2Sum = 0.0
            catWeightSum = 0.0
            catESumSq = 0.0
            catE1SumSq = 0.0
            catE2SumSq = 0.0
            if not test is None:
                sourceCat = afwTable.BaseCatalog.readFits(os.path.join(os.path.join(baseDir, test), "src-%03d.fits"%subfield))
                #butler.get("test_src", {'subfield':subfield, 'epoch':epoch, 'test': test})
            else:
                sourceCat = afwTable.BaseCatalog.readFits(os.path.join(baseDir, "src-%03d.fits"%subfield))
                #sourceCat = butler.get("src", {'subfield':subfield, 'epoch':epoch})
            # make a list of all the error flags
            if len(flagCount) == 0:
                for name in sourceCat.getSchema().getNames():
                    if name.find("_flag") > 0:
                        flagKeys.append(sourceCat.getSchema().find(name).getKey())
                        flagCount.append(0)
                        flagNames.append(name)
            if len(sourceCat) > 0:
                g1 = sourceCat[0].get("g1")
                g2 = sourceCat[0].get("g2")
            #  These are constant values per subfield recorded by galsim

            #   Set up the measurement keys, depending on type (moments or ellip)
            if shape_type == "ellipticity":
                e1Key = sourceCat.getSchema().find(shape_field + "e1").getKey()
                e2Key = sourceCat.getSchema().find(shape_field + "e2").getKey()
            else:
                xxKey = sourceCat.getSchema().find(shape_field + "xx").getKey()
                yyKey = sourceCat.getSchema().find(shape_field + "yy").getKey()
                xyKey = sourceCat.getSchema().find(shape_field + "xy").getKey()

            #   If there is a sigma field, get it.  If it also has a signal field
            #   for signal-to-noise determination, get that too
            if sigma_field == None:
                sigmaKey = None
                sigmaSignalKey = None
            else:
                sigmaKey = sourceCat.getSchema().find(sigma_field).getKey()
                if sigma_signal_field is None:
                    sigmaSignalKey = None
                else:
                    sigmaSignalKey = sourceCat.getSchema().find(sigma_signal_field).getKey()

            #   Now loop through the catalog and summarize the ellipticity measurements
            for source in sourceCat:
                for i, key in enumerate(flagKeys):
                    if source.get(key):
                        flagCount[i] = flagCount[i] + 1
                if shape_type == "ellipticity":
                    e1 = source.get(e1Key)
                    e2 = source.get(e2Key)
                else:
                    xx = source.get(xxKey)
                    yy = source.get(yyKey)
                    xy = source.get(xyKey)
                    e1 = (xx - yy) / (xx + yy)
                    e2 = 2 * xy / (xx + yy)
                e = math.sqrt(e1 * e1 + e2 * e2)

                #   Calculate the rotation in the range [-pi, pi], just to check the pairs
                theta = getStandardAngle(math.atan2(e2,e1)/2.0)
                if sigmaKey == None:
                    weight = 1.0
                else:
                     #   Calculate the ellipticity error if fields are defined to do th
                     #   Use it to calculate the variance, assuming a .25 shape noise
                    sigma = source.get(sigmaKey)
                    if not sigmaSignalKey == None:
                        sigmaSignal = source.get(sigmaSignalKey)
                        sigma = sigma * e / sigmaSignal
                    weight = 1.0/(.25*.25 + sigma*sigma)
                #   Discard nans
                if not e >= 0 and not e <= 0:
                    nanCount = nanCount + 1 
                else:
                    catWeightSum = catWeightSum + weight
                    catESum = catESum + (e * weight)
                    catE1Sum = catE1Sum + (e1 * weight)
                    catE2Sum = catE2Sum + (e2 * weight)
                    catESumSq = catESumSq + (e * e * weight)
                    catE1SumSq = catE1SumSq + (e1 * e1 * weight)
                    catE2SumSq = catE2SumSq + (e2 * e2 * weight)
                    catCount = catCount + 1
            #   For each subfield/epoch, calculate averages
            if catWeightSum > 0:
                 e1Avg = catE1Sum/catWeightSum
                 e2Avg = catE2Sum/catWeightSum
                 eAvg = catESum/catWeightSum
                 e1Stddev = math.sqrt(catE1SumSq/catWeightSum - e1Avg*e1Avg)/math.sqrt(catWeightSum)
                 e2Stddev = math.sqrt(catE2SumSq/catWeightSum - e2Avg*e2Avg)/math.sqrt(catWeightSum)
                 eStddev = math.sqrt(catESumSq/catWeightSum - eAvg*eAvg)
                 print "%d galaxies weight %f.2: e1=%.4f +-%.4f, e2=%.4f +-%.4f, g=%.3f, (%.3f,%.3f)"%(catCount,
                        catWeightSum, e1Avg, e1Stddev, e2Avg, e2Stddev,
                        math.sqrt(g1*g1 + g2*g2), g1, g2)

                 #   And append to the output fits file
                 outrec = outCat.addNew()
                 outrec.set(g1Key, float(g1))
                 outrec.set(g2Key, float(g2))
                 outrec.set(filterKey, config.filter)
                 outrec.set(seeingKey, config.seeing)
                 outrec.set(countKey, catCount)
                 outrec.set(gKey, math.sqrt(g1*g1 + g2*g2))
                 outrec.set(eAvgKey, eAvg)
                 outrec.set(eStdKey, eStddev)
                 outrec.set(e1StdKey, e1Stddev)
                 outrec.set(e2StdKey, e2Stddev)
                 outrec.set(e1AvgKey, e1Avg)
                 outrec.set(e2AvgKey, e2Avg)
                 nAvgs = nAvgs + 1

                 #  Add to the global sums and sums of squares
                 eAvgSum = eAvgSum + eAvg
                 eAvgSumSq = eAvgSumSq + eAvg * eAvg
                 e1AvgSum = e1AvgSum + e1Avg
                 e1AvgSumSq = e1AvgSumSq + e1Avg*e1Avg
                 e2AvgSum = e2AvgSum + e2Avg
                 e2AvgSumSq = e2AvgSumSq + e2Avg*e2Avg
                 count = count + catCount
                 allCount = allCount + len(sourceCat)
                 weightSum = weightSum + catWeightSum
    outCat.writeFits(os.path.join(baseDir, outFile))
    print "Total from all subfields: %d measured out of %d, %d had a nan measuremnt"%(count, allCount, nanCount)
    for i in range(len(flagCount)):
        if flagCount[i] > 0:
            print flagNames[i], ": ", flagCount[i]

if __name__ == "__main__":
#   analyzeShearTest main program:
#   base_dir         Name of the subdirectory containing run_params file
#   out_file         Name of the catalog for output. .fits is appended
#   test             Name of the test directory, if this is a named test

    parser = argparse.ArgumentParser()
    parser.add_argument("base_dir", type=str, help="Name of the directory containing run_params file")
    parser.add_argument("-o", "--out_file", type=str, help="Name of the catalog for output", default=None)
    parser.add_argument("-t", "--test", type=str, help="Type of the test (set of source tables)",
                        default=None)
    args = parser.parse_args()
    if args.out_file is None:
        if args.test is None:
            out_file = "subfields.fits"
        else:
            out_file = "subfields_" + args.test + ".fits"
    #   open the config file for this run
    config = RunShearConfig()
    config.load(os.path.join(args.base_dir, "shear.config"))
    runAnal(args.base_dir, out_file, config, test=args.test)
