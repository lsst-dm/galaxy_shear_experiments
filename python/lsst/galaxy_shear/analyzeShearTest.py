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



def runAnal(base_dir, out_fits, out_sum_fits, params, test=None):

    schema = afwTable.Schema()
    filterKey = schema.addField("filter", type = int, doc = "filter 2 or 3 used in PhoSim psf generator.")
    seeingKey = schema.addField("seeing", type = float, doc = "raw seeing used by PhoSim.")
    shearKey = schema.addField("shear", type = float, doc = "constant shear level used in the GalSim run.")
    countKey = schema.addField("nsource", type = int, doc = "total number of galaxies in this run")
    gKey = schema.addField("g", type = float, doc = "g value calcuted from GalSim g1 and g2")
    g1Key = schema.addField("g1", type = float, doc = "g value calcuted from GalSim g1 and g2")
    g2Key = schema.addField("g2", type = float, doc = "g value calcuted from GalSim g1 and g2")
    eAvgKey = schema.addField("eAvg", type = float, doc = "average e for the sources in this subfield")
    eStdKey = schema.addField("eStd", type = float, doc = "stddev e for the sources in this subfield")
    e1DevKey = schema.addField("e1Dev", type = float, doc = "e1 deviation from g1, for nsub subfields.")
    e1StdKey = schema.addField("e1Std", type = float, doc = "stddev e for the sources in this subfield")

    e2DevKey = schema.addField("e2Dev", type = float, doc = "e2 deviation from g2, for nsub subfields")
    e2StdKey = schema.addField("e2Std", type = float, doc = "stddev e for the sources in this subfield")
    eSumKey = schema.addField("eSum", type = float, doc = "weighted sum of e for all sources")
    e1SumKey = schema.addField("e1Sum", type = float, doc = "weighted sum of e1 for all sources")
    e2SumKey = schema.addField("e2Sum", type = float, doc = "weighted sum of e2 for all sources")
    weightSumKey = schema.addField("weightSum", type = float, doc = "sum of weights for all sources")
    outSourceCat = afwTable.BaseCatalog(schema)
    shear_value = float(params["shear_value"])
    filter =  int(params["filter"])
    seeing  = float(params["seeing"])
    n_subfields = int(params["n_subfields"])
    n_epochs = int(params["n_epochs"])
    galaxy_stamp_size = int(params["galaxy_stamp_size"])
    output_dir = params["output_dir"]

    #   The shape type can be "moments" or "ellipticity"
    #   If we done have an e_sigma, the signal field  tells us to extract the weighting function
    #   from the SNR of some other field, where the sigma will be the error, and the
    #   The sigma used for weighting is shape_field/(sigma_field/sigma_signal_field)
    if "shape_field" in params.keys():
        shape_field = params["shape_field"] + "_"
        shape_type = params["shape_type"]
        if "sigma_field" in params.keys():
            sigma_field = params["sigma_field"]
            if "sigma_signal_field" in params.keys():
                sigma_signal_field = params["sigma_signal_field"]
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
    butler = dafPersist.Butler(root=os.path.join(base_dir, output_dir))

    count = 0
    eSum = 0.0
    e1DevSum = 0.0
    e2DevSum = 0.0
    eSumSq = 0.0
    e1DevSumSq = 0.0
    e2DevSumSq = 0.0
    weightSum = 0.0
    nAvgs = 0
    for subfield in range(n_subfields):
        for epoch in range(n_epochs):
            catCount = 0
            catESum = 0.0
            catE1Sum = 0.0
            catE2Sum = 0.0
            catWeightSum = 0.0
            if not test is None:
                sourceCat = butler.get("test_src", {'subfield':subfield, 'epoch':epoch, 'test': test})
            else:
                sourceCat = butler.get("src", {'subfield':subfield, 'epoch':epoch})
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
                if shape_type == "ellipticity":
                    e1 = source.get(e1Key)
                    e2 = source.get(e2Key)
                else:
                    xx = source.get(xxKey)
                    yy = source.get(yyKey)
                    xy = source.get(xyKey)
                    e1 = .5 * (xx - yy) / (xx + yy)
                    e2 = .5 * 2 * xy / (xx + yy)
                e = math.sqrt(e1 * e1 + e2 * e2)

                #   Calculate the rotation in the range [-pi, pi], just to check the pairs
                theta = math.atan2(e2,e1)/2.0
                thetaApplied = math.atan2(g2,g1)/2.0
                if theta >= math.pi:
                    theta = theta - math.pi
                if theta <= -math.pi:
                    theta = theta + math.pi
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
                    continue
                catWeightSum = catWeightSum + weight
                catESum = catESum + (e * weight)
                catE1Sum = catE1Sum + (e1 * weight)
                catE2Sum = catE2Sum + (e2 * weight)
                catCount = catCount + 1

            #   For each subfield/epoch, calculate averages
            if catWeightSum > 0:
                 e1Avg = catE1Sum/catWeightSum
                 e2Avg = catE2Sum/catWeightSum
                 eAvg = math.sqrt(e1Avg*e1Avg + e2Avg*e2Avg)
                 e1Dev = e1Avg
                 e2Dev = e2Avg
                 #  Add to the global sums and sums of squares
                 eSum = eSum + eAvg
                 eSumSq = eSumSq + eAvg * eAvg
                 e1DevSum = e1DevSum + e1Dev
                 e1DevSumSq = e1DevSumSq + e1Dev*e1Dev
                 e2DevSum = e2DevSum + e2Dev
                 e2DevSumSq = e2DevSumSq + e2Dev*e2Dev
                 count = count + catCount
                 weightSum = weightSum + catWeightSum
                 #   And append to the output fits file
                 outrec = outSourceCat.addNew()
                 outrec.set(g1Key, float(g1))
                 outrec.set(g2Key, float(g2))
                 outrec.set(filterKey, filter)
                 outrec.set(seeingKey, seeing)
                 outrec.set(shearKey, shear_value)
                 outrec.set(countKey, catCount)
                 outrec.set(gKey, math.sqrt(g1*g1 + g2*g2))
                 outrec.set(eAvgKey, eAvg)
                 outrec.set(e1DevKey, e1Avg-g1)
                 outrec.set(e2DevKey, e2Avg-g2)
                 outrec.set(eSumKey, catESum)
                 outrec.set(e1SumKey, catE1Sum)
                 outrec.set(e2SumKey, catE2Sum)
                 outrec.set(weightSumKey, catWeightSum)
                 nAvgs = nAvgs + 1
    outSourceCat.writeFits(out_fits)
    #   Now summarize the averages and stddevs for all the subfield/epochs
    if nAvgs > 2:
        eStddev = math.sqrt( (nAvgs * eSumSq - eSum * eSum)/(nAvgs * (nAvgs - 1)))
        e1Stddev = math.sqrt((nAvgs * e1DevSumSq - e1DevSum * e1DevSum)/(nAvgs * (nAvgs - 1)))
        e2Stddev = math.sqrt((nAvgs * e2DevSumSq - e2DevSum * e2DevSum)/(nAvgs * (nAvgs - 1)))
        print "%d subfields: e=%.4f +-%.4f, e1=%.4f +-%.4f, e2=%.4f +-%.4f, g=%.3f, (%.3f,%.3f)"%(nAvgs,
               eSum/nAvgs, eStddev, e1DevSum/nAvgs, e1Stddev, e2DevSum/nAvgs, e2Stddev,
               math.sqrt(g1*g1 + g2*g2), g1, g2)
    else:
        print "Not enough subfields to calculate error.  Subfield not added to catalog"
        print "%d subfields: e = %.4f, e1 = %.4f, e2 = %.4f, g = %.3f"%(nAvgs,
               eSum/nAvgs, e1DevSum/nAvgs, e2DevSum/nAvgs,
               math.sqrt(g1*g1 + g2*g2))
    
    sumSourceCat = afwTable.BaseCatalog(outSourceCat.getSchema())
    sumrec = sumSourceCat.addNew()
    sumrec.set(g1Key, float(g1))
    sumrec.set(g2Key, float(g2))
    sumrec.set(filterKey, filter)
    sumrec.set(seeingKey, seeing)
    sumrec.set(shearKey, shear_value)
    sumrec.set(countKey, count)
    sumrec.set(gKey, math.sqrt(g1*g1 + g2*g2))
    sumrec.set(eAvgKey, eSum/nAvgs)
    sumrec.set(eStdKey, eStddev)
    sumrec.set(e1DevKey, e1DevSum/nAvgs)
    sumrec.set(e1StdKey, e1Stddev)
    sumrec.set(e2DevKey, e2DevSum/nAvgs)
    sumrec.set(e2StdKey, e2Stddev)
    sumrec.set(eSumKey, eSum)
    sumrec.set(e1SumKey, e1DevSum)
    sumrec.set(e2SumKey, e2DevSum)
    sumrec.set(e2SumKey, catE2Sum)
    sumrec.set(weightSumKey, weightSum)
    sumSourceCat.writeFits(out_sum_fits)
    count = 0
    eSum = 0.0
    e1DevSum = 0.0
    e2DevSum = 0.0
    eSumSq = 0.0
    e1DevSumSq = 0.0
    e2DevSumSq = 0.0

if __name__ == "__main__":
    """

    base_dir         Name of the subdirectory containing run_params file
    out_cat          Name of the catalog for output. .fits is appended
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("basedir", type=str, help="Name of the directory containing run_params file")
    parser.add_argument("-o", "--outfile", type=str, help="Name of the catalog for output", default=None)
    parser.add_argument("-t", "--test", type=str, help="Type of the test (set of source tables)", default=None)
    args = parser.parse_args()
    if args.outfile is None:
        if args.test is None:
            outfile = os.path.join(args.basedir, "subfields_" + args.test + ".fits")
        else:
            outfile = os.path.join(args.basedir, "subfields.fits")

    #   open the run_params file for this run and extract what we need.
    fin = open(os.path.join(args.base_dir, "run_params"), 'r')
    params = dict()
    for line in fin:
        if line.isspace():
            continue
        inp = line.split()
        if inp[0][0] == '#':
            continue
        params[inp[0]] = inp[1]
    runAnal(args.base_dir, args.out_cat, params, test=args.test)
