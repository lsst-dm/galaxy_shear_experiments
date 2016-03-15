#/!/usr/bin/env python
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
import random

import lsst.daf.persistence as dafPersist
import lsst.afw.table as afwTable

from lsst.galaxy_shear.shearConfig import RunShearConfig

# do the analysis on a single source catalog, and return the statistics
# if flagCounts is defined, an accumulator is kept of

def analyzeCat(sourceCat, config, flagCounts=None, flagKeys=None, flagNames=None):
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

    # These are the "catalog" accumulators.  Each catalog represents a single
    # subfield/epoch, which is essentially a single pointing with constant (g1,g2)
    catCount = 0
    catESum = 0.0
    catE1Sum = 0.0
    catE2Sum = 0.0
    catWeightSum = 0.0
    catESumSq = 0.0
    catE1SumSq = 0.0
    catE2SumSq = 0.0
    catNanCount = 0
    #   Set up the measurement keys, depending on type (moments or ellipticity)
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
        if not flagCounts is None:
            for i, key in enumerate(flagKeys):
                if source.get(key):
                    flagCounts[i] = flagCounts[i] + 1
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
            catNanCount = catNanCount + 1
        else:
            catWeightSum = catWeightSum + weight
            catESum = catESum + (e * weight)
            catE1Sum = catE1Sum + (e1 * weight)
            catE2Sum = catE2Sum + (e2 * weight)
            catESumSq = catESumSq + (e * e * weight)
            catE1SumSq = catE1SumSq + (e1 * e1 * weight)
            catE2SumSq = catE2SumSq + (e2 * e2 * weight)
            catCount = catCount + 1

    if catCount == 0:
        raise Exception("catalog found with no good records")
    e1Avg = catE1Sum/catWeightSum
    e2Avg = catE2Sum/catWeightSum
    eAvg = catESum/catWeightSum
    eStddev = math.sqrt(catE1SumSq/catWeightSum - e1Avg*e1Avg)/math.sqrt(catWeightSum)
    e1Stddev = math.sqrt(catE1SumSq/catWeightSum - e1Avg*e1Avg)/math.sqrt(catWeightSum)
    e2Stddev = math.sqrt(catE2SumSq/catWeightSum - e2Avg*e2Avg)/math.sqrt(catWeightSum)
    print "%d good, weight %.2f: e1=%.4f +-%.4f, e2=%.4f +-%.4f"%(catCount,
             catWeightSum, e1Avg, e1Stddev, e2Avg, e2Stddev)
    return e1Avg, e1Stddev, e2Avg, e2Stddev, eAvg, eStddev, catCount, catNanCount

# Resample the source catalog and return as a newCat
# Note that the resampling has to be done in pairs, so we create an index of each
# galaxy and its pair using information from the epoch_catalog created by great3sims
# When we select a galaxy from the sampling file, we must also take its pair.
def resampleCat(sourceCat, epochFile, lookupFile):
    if not os.path.isfile(lookupFile):
        makeLookup(epochFile, lookupFile)
    lookupCat = afwTable.BaseCatalog.readFits(lookupFile)
    newCat = afwTable.SourceCatalog(sourceCat)
    length = len(sourceCat)
    lookupDict = dict()
    galKey = lookupCat.getSchema().find("gal").key
    pairKey = lookupCat.getSchema().find("pair").key
    # create a has from the lookup catalog, saving both directions of the pairing.
    # some catalogs only have one of the pair.
    for lookup in lookupCat:
        lookupDict[lookup.get(galKey)] = lookup.get(pairKey)
        lookupDict[lookup.get(pairKey)] = lookup.get(galKey)
    for i in range(0, length, 2):
        index = random.randint(0, length - 1)
        pair_index = lookupDict[index]
        newCat[i] = sourceCat[index]
        newCat[i+1] = sourceCat[pair_index]
    return newCat

def isNear(value1, value2, diff, tol):
     return abs(abs(value1 - value2) - diff) < tol

# Make a lookup of galaxies and their pairs.  Save to disk as it may be needed later
def makeLookup(epochCatName, lookupCatName):
    epochCat = afwTable.BaseCatalog.readFits(epochCatName)
    schema = epochCat.getSchema()
    lookup = dict()
    keys = []
    for name in ("cosmos_ident", "gal_sn", "bulge_n", "bulge_hlr", "bulge_q", "bulge_flux", "disk_hlr", "disk_q", "disk_flux"):
        keys.append(schema.find(name).key)
    bulgeThetaKey = schema.find("bulge_beta_radians").key
    diskThetaKey = schema.find("disk_beta_radians").key

    for i in range(len(epochCat)):
        if i in lookup.values():
            continue
        source = epochCat[i]
        for j in range(i+1, len(epochCat)):
            fail = False
            for key in keys:
                if not source.get(key) == epochCat[j].get(key):
                    fail = True
                    break
            if fail:
                continue
            thetai = epochCat[i].get(bulgeThetaKey)
            thetaj = epochCat[j].get(bulgeThetaKey)
            if not isNear(thetai, thetaj, math.pi/2.0, .05):
                continue
            thetai = epochCat[i].get(diskThetaKey)
            thetaj = epochCat[j].get(diskThetaKey)
            if not thetai == 0.0 and not isNear(thetai, thetaj, math.pi/2.0, .05):
                continue
            if j in lookup.values():
                continue
            lookup[i] = j
            lookup[j] = i
            break
        if not i in lookup.keys():
            print lookupCatName, "Failed to find match", i
            lookup[i] = -1
    schema = afwTable.Schema()
    galKey = schema.addField("gal", type = int, doc = "index of record in a source catalog")
    pairKey = schema.addField("pair", type = int, doc = "index of record in a source catalog")
    lookupCat = afwTable.BaseCatalog(schema)
    for key in lookup.keys():
        pair = lookup[key]
        rec = lookupCat.addNew()
        rec.set(galKey, key)
        rec.set(pairKey, pair)
    lookupCat.writeFits(lookupCatName)

def runAnal(baseDir, outFile, config, test=None, resample=False):
    #  Create a summary table of the src*.fits tables in baseDir
    schema = afwTable.Schema()
    filterKey = schema.addField("filter", type = int, doc = "filter 2 or 3 used in PhoSim psf generator.")
    seeingKey = schema.addField("seeing", type = float, doc = "raw seeing used by PhoSim.")
    shearValueKey = schema.addField("shear_value", type = float,
                                    doc = "constant shear level used in the GalSim run.")
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

    #  These are the "overall" accumulators (sums and counts for all catalogs)
    count = 0
    nanCount = 0
    allCount = 0
    flagKeys = []
    flagNames = []
    flagCounts = []

    for subfield in range(config.n_subfields):
        for epoch in range(config.n_epochs):
            # Open the src.fits file for the current subfield and epoch.
            if not test is None:
                sourceFile = os.path.join(os.path.join(baseDir, test), "src-%03d.fits"%subfield)
            else:
                sourceFile = os.path.join(baseDir, "src-%03d.fits"%subfield)
            #skip any subfields where the src.fits file does not exist
            if not os.path.exists(sourceFile):
                print "src.fits file does not exist for %03d"%subfield
                continue
            sourceCat = afwTable.SourceCatalog.readFits(sourceFile)
            # save g1,g2 for later printouts
            #  These are constant values per subfield recorded by galsim
            if len(sourceCat) > 0:
                g1 = sourceCat[0].get("g1")
                g2 = sourceCat[0].get("g2")

            # make a list of all the error flags
            if len(flagCounts) == 0:
                for name in sourceCat.getSchema().getNames():
                    if name.find("_flag") > 0:
                        flagKeys.append(sourceCat.getSchema().find(name).getKey())
                        flagCounts.append(0)
                        flagNames.append(name)
            # bootsitrap resample the source catalog using if requested
            if resample:
                epochFile = os.path.join(os.path.join(baseDir), "control", "ground", "constant",
                                         "epoch_catalog-%03d-%01d.fits"%(subfield, epoch))
                lookupFile = os.path.join(os.path.join(baseDir), "control", "ground", "constant",
                                          "lookup-%03d-%01d.fits"%(subfield, epoch))
                sourceCat = resampleCat(sourceCat, epochFile, lookupFile)
            (e1Avg, e1Stddev, e2Avg, e2Stddev, eAvg, eStddev, catCount, catNanCount) = \
                analyzeCat(sourceCat, config, flagCounts, flagKeys, flagNames)

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

            nanCount = nanCount + catNanCount
            count = count + catCount
            allCount = allCount + len(sourceCat)

    outCat.writeFits(os.path.join(baseDir, outFile))
    print "Total from all subfields: %d measured out of %d, %d had a nan measuremnt"%(count,
           allCount, nanCount)
    for i in range(len(flagCounts)):
        if flagCounts[i] > 0:
            print flagNames[i], ": ", flagCounts[i]

if __name__ == "__main__":
#   analyzeShearTest main program:
#   base_dir         Name of the subdirectory containing run_params file
#   out_file         Name of the catalog for output. .fits is appended
#   test             Name of the test directory, if this is a named test

    parser = argparse.ArgumentParser()
    parser.add_argument("base_dir", type=str, help="Name of the directory containing run_params file")
    parser.add_argument("-o", "--out_file", type=str, help="Name of the catalog for output", default=None)
    parser.add_argument("-r", "--resample", type=int,
                         help="Resample source catalog, append this number to anal filename",
                         default=None)
    parser.add_argument("-t", "--test", type=str, help="Type of the test (set of source tables)",
                        default=None)
    args = parser.parse_args()

    # if the output file  name is given, use it.  Otherwise, name it anal.fits or anal_test.fits
    if args.out_file is None:
        if args.test is None:
            out_file = "anal.fits"
        else:
            out_file = "anal_" + args.test + ".fits"
    #   open the config file for this run
    config = RunShearConfig()
    config.load(os.path.join(args.base_dir, "shear.config"))
    runAnal(args.base_dir, args.out_file, config, test=args.test, resample=args.resample)
