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
if __name__ == "__main__":
    """
    This main program reads data for multiple runs in filter, seeing, and shear value, and

    analysis_file         Name of the fits file containing the analysis data
    seeing                plot only runs with the specified seeing
    filter                plot only runs with the specified filter
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("analysis_file", type=str, help="Analysis fits file")
    parser.add_argument("-f", "--filter", type=int, help="filter to be used for plotting", default=2)
    parser.add_argument("-s", "--seeing", type=float, help="seeing to be used for plotting", default=0.7)
    args = parser.parse_args()
    outSourceCat = afwTable.BaseCatalog.readFits(args.analysis_file)

    # Keys for fetching data from the output table
    schema = outSourceCat.getSchema()
    filterKey = schema.find("filter").getKey()
    seeingKey = schema.find("seeing").getKey()
    shearKey = schema.find("shear").getKey()
    nsourceKey = schema.find("nsource").getKey()
    eAvgKey = schema.find("eAvg").getKey()
    eStdKey = schema.find("eStd").getKey()
    shear = []
    e = []
    eStd = []

    # Go through the catalog and get the values requested in the plot
    for rec in outSourceCat:
        filter = rec.get(filterKey)
        seeing = rec.get(seeingKey)
        if filter == args.filter and seeing == args.seeing:
            shear.append(rec.get(shearKey))
            e.append(rec.get(eAvgKey))
            eStd.append(rec.get(eStdKey))
    shear = numpy.array(shear, dtype=float)
    e = numpy.array(e, dtype=float)
    eStd = numpy.array(eStd, dtype=float)

    m, b = numpy.polyfit(shear, e, 1)
    label = "seeing = %.4f, m = %.3f, b = %.3f"%(seeing, m, b)
    figure = plot.figure()
    figure.suptitle(label)
    ax = figure.add_subplot(111)
    ax.set_xlabel("applied shear")
    ax.set_ylabel("measured shear")
    # plot the data with error bars
    ax.errorbar(shear, e, yerr=eStd, marker = "o", linestyle='None', label =label)
    # plot a first order fit
    ax.plot(shear, m*shear + b)
    plot.show()