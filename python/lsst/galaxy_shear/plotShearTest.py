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
import matplotlib.pyplot as plot
import math
import sys
import os.path
import numpy
import pyfits

import lsst.afw.math as afwMath
import lsst.afw.table as afwTable

"""
plotShearTest is a python program which is used make a applied vs. measured shear plot for
the shear_test_experiments simulations, and to also display error bars for the stddev's
calculated doing mulitple subfields with the same applied shear.

It is assumed that previous scripts have created an analysis_file which has statistical
summaries of the shear and shear error as a function of seeing, filter, and applied shear.

This program simply takes the values for a given filter and seeing, fits the points
to a line, and graphs the measured values and errors and the regression line.
"""

def runPlot(useFilter, useSeeing, useIndex=None):

    outSourceCat = afwTable.BaseCatalog.readFits(args.analysis_file)

    # Keys for fetching data from the output table
    schema = outSourceCat.getSchema()
    filterKey = schema.find("filter").getKey()
    seeingKey = schema.find("seeing").getKey()
    shearKey = schema.find("shear_value").getKey()
    nsourceKey = schema.find("nsource").getKey()
    eAvgKey = schema.find("eAvg").getKey()
    eStdKey = schema.find("eStd").getKey()
    e1AvgKey = schema.find("e1Avg").getKey()
    e1StdKey = schema.find("e1Std").getKey()
    e2AvgKey = schema.find("e2Avg").getKey()
    e2StdKey = schema.find("e2Std").getKey()
    countKey = schema.find("nsource").getKey()
    g1Key = schema.find("g1").getKey()
    g2Key = schema.find("g2").getKey()
    g1 = []
    g2 = []
    e = []
    eStd = []
    e1Avg = []
    e1Std = []
    e2Avg = []
    e2Std = []
    nSources = 0
    # Go through the catalog and get the values requested in the plot
    for rec in outSourceCat:
        filter = rec.get(filterKey)
        seeing = rec.get(seeingKey)
        if filter == useFilter and seeing == useSeeing:
            g1.append(rec.get(g1Key))
            g2.append(rec.get(g2Key))
            nSources = nSources + rec.get(countKey)
            e1Avg.append(rec.get(e1AvgKey))
            e1Std.append(rec.get(e1StdKey))
            e2Avg.append(rec.get(e2AvgKey))
            e2Std.append(rec.get(e2StdKey))
            print "Stddev e1, e2: ", rec.get(e1StdKey), rec.get(e2StdKey)
            #print rec.get(eAvgKey), rec.get(e1AvgKey), rec.get(e2AvgKey)
            #print rec.get(eStdKey), rec.get(e1StdKey), rec.get(e2StdKey)
    g1 = numpy.array(g1, dtype=float)
    g2 = numpy.array(g2, dtype=float)
    eStd = numpy.array(eStd, dtype=float)
    e1Avg = numpy.array(e1Avg, dtype=float)
    e1Std = numpy.array(e1Std, dtype=float)
    e2Avg = numpy.array(e2Avg, dtype=float)
    e2Std = numpy.array(e2Std, dtype=float)
    
    m1, b1 = numpy.polyfit(g1, e1Avg, 1)
    m2, b2 = numpy.polyfit(g2, e2Avg, 1)
    label = "%d gals, f%d_%.1f, e1m: %.4f b: %.4f) e2m: %.4f b:%.4f)"%(nSources, filter, seeing, m1, b1, m2, b2)
    figure = plot.figure()
    figure.suptitle(label)
    ax = figure.add_subplot(111)
    ax.set_xlabel("applied shear")
    ax.set_ylabel("measured shear")
    # plot the data with error bars
    ax.errorbar(g1, e1Avg, yerr=e1Std, marker = "o", linestyle='None', label =label, color='red')
    ax.errorbar(g2, e2Avg, yerr=e2Std, marker = "o", linestyle='None', label =label)
    # plot a first order fit
    ax.plot(g1, m1*g1 + b1)
    ax.plot(g2, m2*g2 + b2)
    plot.show()

if __name__ == "__main__":

#   This main program reads data for multiple runs in filter, seeing, and shear value, and
#
#   analysis_file         Name of the fits file containing the analysis data
#   seeing                plot only runs with the specified seeing
#   filter                plot only runs with the specified filter

    parser = argparse.ArgumentParser()
    parser.add_argument("analysis_file", type=str, help="Analysis fits file")
    parser.add_argument("-f", "--filter", type=int, help="filter to be used for plotting", default=3)
    parser.add_argument("-b", "--index", type=int, help="indexs of filter to be used for plotting", default=None)
    parser.add_argument("-s", "--seeing", type=float, help="seeing to be used for plotting", default=0.5)
    args = parser.parse_args()
    runPlot(args.filter, args.seeing, useIndex=None)
