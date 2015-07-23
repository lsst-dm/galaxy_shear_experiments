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
import os.path

import lsst.daf.persistence as dafPersist
import lsst.afw.table as afwTable
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
def getCatalog(catName):
    #  This routine processes one run at a time, and appends to an existing catalog.
    if os.path.exists(catName):
        outCat = afwTable.BaseCatalog.readFits(catName)
    else:
    #   if the catalog has not been created yet, create an empy one with the right schema.
        schema = afwTable.Schema()
        schema.addField("filter", type = int, doc = "")
        schema.addField("seeing", type = float, doc = "")
        schema.addField("shear", type = float, doc = "")
        schema.addField("nsub", type = int, doc = "")
        schema.addField("nsource", type = int, doc = "")
        schema.addField("g", type = float, doc = "")
        schema.addField("eAvg", type = float, doc = "")
        schema.addField("eStd", type = float, doc = "")
        schema.addField("e1Dev", type = float, doc = "")
        schema.addField("e1Std", type = float, doc = "")
        schema.addField("e2Dev", type = float, doc = "")
        schema.addField("e2Std", type = float, doc = "")
        outCat = afwTable.BaseCatalog(schema)
    return outCat

if __name__ == "__main__":
    """

    base_dir         Name of the subdirectory containing run_params file
    out_cat          Name of the catalog for output. .fits is appended
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("base_dir", type=str, help="Name of the directory containing run_params file")
    parser.add_argument("out_cat", type=str, help="Name of the catalog for output")
    args = parser.parse_args()

    outSourceCat = getCatalog(args.out_cat + ".fits")
    schema = outSourceCat.getSchema()

    #   find all the keys we need to output the summary of each subfield
    filterKey = schema.find("filter").getKey()
    seeingKey = schema.find("seeing").getKey()
    shearKey = schema.find("shear").getKey()
    nsubKey = schema.find("nsub").getKey()
    nsourceKey = schema.find("nsource").getKey()
    gKey = schema.find("g").getKey()
    eAvgKey = schema.find("eAvg").getKey()
    eStdKey = schema.find("eStd").getKey()
    e1DevKey = schema.find("e1Dev").getKey()
    e1StdKey = schema.find("e1Std").getKey()
    e2DevKey = schema.find("e2Dev").getKey()
    e2StdKey = schema.find("e2Std").getKey()

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

    shear_value = float(params["shear_value"])
    psf_file = os.path.basename(params["psf_dir"])
    filter =  int(params["filter"])
    seeing  = float(params["seeing"])
    n_subfields = int(params["n_subfields"])
    n_epochs = int(params["n_epochs"])
    galaxy_stamp_size = int(params["galaxy_stamp_size"])
    output_dir = params["output_dir"]

    #   This tells us the algorithm fields we need to extract the ellipticity and error
    #   The shape type can be "moments" or "ellipticity"
    #   The error can be extracted from an ellipticity error field or based on the
    #   signal-to-noise of some other field such as flux and flux error.  In that
    #   case, the weighting value used is shape_field/(sigma_field/sigma_signal_field)
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
    butler = dafPersist.Butler(root=os.path.join(args.base_dir, output_dir))

    count = 0
    eSum = 0.0
    e1DevSum = 0.0
    e2DevSum = 0.0
    eSumSq = 0.0
    e1DevSumSq = 0.0
    e2DevSumSq = 0.0
    thetaSum = 0.0
    for subfield in range(n_subfields):
        for epoch in range(n_epochs):
            catCount = 0
            catE1DevSum = 0.0
            catE2DevSum = 0.0
            catWeightSum = 0.0
            sourceCat = butler.get("src", {'subfield':subfield, 'epoch':epoch})

            #  These are constant values per subfield recorded by galsim
            g1Key = sourceCat.getSchema().find("g1").getKey()
            g2Key = sourceCat.getSchema().find("g2").getKey()

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
                if sigma_signal_field == None:
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
                    e1 = (xx - yy) / (xx + yy)
                    e2 = 2 * xy / (xx + yy)
                e = math.sqrt(e1 * e1 + e2 * e2)
                g1 = source.get(g1Key)
                g2 = source.get(g2Key)

                #   Calculate the rotation in the range [-pi, pi], just to check the pairs
                theta = math.atan2(e2,e1)/2.0
                if theta >= math.pi:
                    theta = theta - math.pi
                if theta <= -math.pi:
                    theta = theta + math.pi
                if sigmaKey == None:
                    weight = 1.0
                else:
                    #   Calculate the ellipticity error if fields are defined to do this
                    #   Use it to calculate the variance, assuming a .25 shape noise
                    sigma = source.get(sigmaKey)
                    if not sigmaSignalKey == None:
                        sigmaSignal = source.get(sigmaSignalKey)
                        sigma = sigma * e / sigmaSignal
                    weight = 1.0/(.25*.25 + sigma*sigma)

                #   Discard nans
                if not e >= 0 and not e <= 0:
                    continue
                thetaSum = thetaSum + theta
                catWeightSum = catWeightSum + weight
                catE1DevSum = catE1DevSum + (e1 * weight)
                catE2DevSum = catE2DevSum + (e2 * weight)
                catCount = catCount + 1

            #   For each subfield/epoch, calculate averages
            e1Dev = (catE1DevSum - g1)/catWeightSum
            e2Dev = (catE2DevSum - g2)/catWeightSum
            eavg = math.sqrt((catE1DevSum*catE1DevSum + catE2DevSum*catE2DevSum))/catWeightSum
            print "    cat average of %d: e = %.4f, de1 = %.4f, de2 = %.4f, with g1 = %.4f, g2 = %.4f"%(catCount, eavg, e1Dev, e2Dev, g1, g2)
            #  Add to the global sums and sums of squares
            eSum = eSum + eavg
            e1DevSum = e1DevSum + e1Dev
            e2DevSum = e2DevSum + e2Dev
            eSumSq = eSumSq + eavg * eavg
            e1DevSumSq = e1DevSumSq + e1Dev*e1Dev
            e2DevSumSq = e2DevSumSq + e2Dev*e2Dev
            count = count + catCount

    #   Now summarize the averages and stddevs for all the subfield/epochs
    nAvgs = n_subfields * n_epochs
    if nAvgs > 2:
        eStddev = math.sqrt( (nAvgs * eSumSq - eSum * eSum)/(nAvgs * (nAvgs - 1)))
        e1Stddev = math.sqrt((nAvgs * e1DevSumSq - e1DevSum * e1DevSum)/(nAvgs * (nAvgs - 1)))
        e2Stddev = math.sqrt((nAvgs * e2DevSumSq - e2DevSum * e2DevSum)/(nAvgs * (nAvgs - 1)))
        print "%d subfields: e = %.4f +- %.4f, e1 = %.4f +- %.4f, e2 = %.4f +- %.4f, g = %.3f %.3f"%(nAvgs,
               eSum/nAvgs, eStddev, e1DevSum/nAvgs, e1Stddev, e2DevSum/nAvgs, e2Stddev,
               math.sqrt(g1*g1 + g2*g2), thetaSum)
        #   And append to the output fits file
        outrec = outSourceCat.addNew()
        outrec.set(filterKey, filter)
        outrec.set(seeingKey, seeing)
        outrec.set(shearKey, shear_value)
        outrec.set(nsubKey, nAvgs)
        outrec.set(nsourceKey, count)
        outrec.set(gKey, math.sqrt(g1*g1 + g2*g2))
        outrec.set(eAvgKey, eSum/nAvgs)
        outrec.set(eStdKey, eStddev)
        outrec.set(e1DevKey, e1DevSum/nAvgs)
        outrec.set(e1StdKey, e1Stddev)
        outrec.set(e2DevKey, e2DevSum/nAvgs)
        outrec.set(e2StdKey, e2Stddev)
    else:
        print "Not enough subfields to calculate and error.  Not added to catalog"
        print "%d subfields: e = %.4f, e1 = %.4f, e2 = %.4f, g = %.3f"%(nAvgs,
               eSum/nAvgs, e1DevSum/nAvgs, e2DevSum/nAvgs,
               math.sqrt(g1*g1 + g2*g2))
    outSourceCat.writeFits(args.out_cat + ".fits")

