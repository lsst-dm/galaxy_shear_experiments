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
import glob
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

def runPlot(analysis_file, useFilter, useSeeing, useIndex=None, display=1):

    # input file is the analyis of a set of subfields, all from the same measurement run
    subfieldCat = afwTable.BaseCatalog.readFits(analysis_file)
    # Keys for fetching data from the output table
    schema = subfieldCat.getSchema()
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

    # accumulate vectors of values for all the subfields in this run
    g1 = []
    g2 = []
    e = []
    eStd = []
    e1Avg = []
    e1Std = []
    e2Avg = []
    e2Std = []
    nSources = 0
    # Accumulate sums for likelikelihood calculation, which can then be used to get error
    # for the slope and intercept of the 1D fit.  This is using the formulas in Numerical
    # Recipes Chapter 15.2
    s1 = 0
    s1x = 0
    s1y = 0
    s1xx = 0
    s1xy = 0
    s2 = 0
    s2x = 0
    s2y = 0
    s2xx = 0
    s2xy = 0
 
    # Go through the catalog and get the values requested in the plot
    for rec in subfieldCat:
        filter = rec.get(filterKey)
        seeing = rec.get(seeingKey)
        if (useFilter is None or filter == useFilter) and (useSeeing is None or seeing == useSeeing):
            x1 = rec.get(g1Key)
            g1.append(x1)
            x2 = rec.get(g2Key)
            g2.append(x2)
            nSources = nSources + rec.get(countKey)
            # each subfield catalog contains and estimated mean and std error of the mean for e1 and e2
            y1 = rec.get(e1AvgKey)
            sigma1 = rec.get(e1StdKey)
            e1Avg.append(y1)
            e1Std.append(sigma1)
            w1 = 1.0/(sigma1 * sigma1)
            s1 = s1 + w1
            s1x = s1x + x1 * w1
            s1y = s1y + y1 * w1
            s1xx = s1xx + x1 * x1 * w1
            s1xy = s1xy + x1 * y1 * w1

            y2 = rec.get(e2AvgKey)
            sigma2 = rec.get(e2StdKey)
            e2Avg.append(y2)
            e2Std.append(sigma2)
            w2 = 1.0/(sigma2 * sigma2)
            s2 = s2 + w2
            s2x = s2x + x2 * w2
            s2y = s2y + y2 * w2
            s2xx = s2xx + x2 * x2 * w2
            s2xy = s2xy + x2 * y2 * w2

    #  Now calculate the slope and intercept, as well as their errors
    delta1 = (s1 * s1xx) - (s1x * s1x) 
    b1 = (s1xx * s1y - s1x * s1xy)/delta1
    m1 = (s1 * s1xy - s1x * s1y)/delta1
    m1err = math.sqrt(s1/delta1)
    b1err = math.sqrt(s1xx/delta1)

    delta2 = (s2 * s2xx) - (s2x * s2x)
    b2 = (s2xx * s2y - s2x * s2xy)/delta2
    m2 = (s2 * s2xy - s2x * s2y)/delta2
    m2err = math.sqrt(s2/delta2)
    b2err = math.sqrt(s2xx/delta2)

    nPoints = len(g1)
    
    # print the results to the console
    print "Fits for %s %d/%d sources, m1:%.6f+-%.6f, b1: %.6f+-%.6f m2:%.6f+-%.6f b2:%.6f+-%.6f"%(analysis_file,
        nSources, nPoints, m1, m1err, b1, b1err, m2, m2err, b2, b2err)
    # plot a first order fit
    if display:
        import matplotlib.pyplot as plot
        label = "Avg %.1f gals/subfield, m1=%.4f b1=%.4f, m2=%.4f b2=%.4f"%(float(nSources)/nPoints, m1, b1, m2, b2)
        g1 = numpy.array(g1, dtype=float)
        g2 = numpy.array(g2, dtype=float)
        eStd = numpy.array(eStd, dtype=float)
        e1Avg = numpy.array(e1Avg, dtype=float)
        e1Std = numpy.array(e1Std, dtype=float)
        e2Avg = numpy.array(e2Avg, dtype=float)
        e2Std = numpy.array(e2Std, dtype=float)
        zeros = numpy.zeros(nPoints, dtype=float)
        figure = plot.figure()
        figure.suptitle(label)
        ax = figure.add_subplot(111)
        ax.set_xlabel("applied shear")
        ax.set_ylabel("measured shear")
        # plot the data with error bars
        ax.errorbar(g1, e1Avg, yerr=e1Std, marker = ".", linestyle='None', label =label, color='red')
        ax.errorbar(g2, e2Avg, yerr=e2Std, marker = ".", linestyle='None', label =label, color="blue")
        ax.plot(g1, m1*g1 + b1, color="black")
        ax.plot(g2, m2*g2 + b2, color="cyan")
        plot.show()

    return nPoints, m1, m1err, b1, b1err, m2, m2err, b2, b2err

if __name__ == "__main__":

#   This main program reads data for multiple runs in filter, seeing, and shear value, and
#
#   analysis_file         Name of the fits file containing the analysis data
#   seeing                plot only runs with the specified seeing
#   filter                plot only runs with the specified filter

    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--analysis_file", type=str, help="Analysis fits file", default=None) 
    parser.add_argument("-t", "----test", type=str, help="type of test: s, n, r", default=None)
    parser.add_argument("-f", "--filter", type=int, help="filter to be used for plotting", default=None)
    parser.add_argument("-b", "--index", type=int, help="index of filter to be used for plotting", default=None)
    parser.add_argument("-s", "--seeing", type=float, help="seeing to be used for plotting", default=None)
    parser.add_argument("-d", "--display", type=int, help="display the resulting plot", default=1)
    args = parser.parse_args()

    #  If this is must a single analysis file, just display the fit
    if args.test is None:
        nPoints, m1, m1err, b1, b1err, m2, m2err, b2, b2err = runPlot(args.analysis_file, args.filter, args.seeing, useIndex=None, display=args.display)
    #  But if it is a test involving multiple runs, call runPlot multiple times, then display a graph of all the fits.
    else:
        axis_type = "number"
        if args.test == "s":
            test_type = "stampsize" 
        if args.test == "n":
            test_type = "nGrowFootprint" 
        if args.test == "r":
            test_type = "nInitialRadii"
        if args.test == "a":
            axis_type = "string"
            test_type = "PsfShapeletApprox"
        files = glob.glob("%s*/subfields.fits"%(args.test))
                 
        m1s = []
        m2s = []
        b1s = []
        b2s = []
        m1errs = []
        m2errs = []
        b1errs = []
        b2errs = []
        testnumbers = []
        testlabels = []
        for i, file in enumerate(files):
            if args.display == 'all':
                display = True
            else:
                display = False 
            nPoints, m1, m1err, b1, b1err, m2, m2err, b2, b2err = runPlot(file, args.filter, args.seeing, useIndex=None, display=display)
            m1s.append(m1)
            m2s.append(m2)
            b1s.append(b1)
            b2s.append(b2)
            m1errs.append(m1err)
            m2errs.append(m2err)
            b1errs.append(b1err)
            b2errs.append(b2err)
            start = 0
            end = file.find("/")
            name = file[start:end] 
            
            if axis_type == "number":
                testnumbers.append(float(name))
                testlabels.append(float(name))
            else:
                testnumbers.append(float(i))
                testlabels.append("")
                testlabels.append(name)

        if args.display:
            import matplotlib.pyplot as plot
            figure, (ax1, ax2, ax3, ax4) = plot.subplots(nrows=4, ncols=1)

            figure.suptitle("Shear estimation bias vs. %s"%test_type)

            ax1.set_title("m1", size="small") 
            ax2.set_title("c1", size="small") 
            ax3.set_title("m2", size="small") 
            ax4.set_title("c2", size="small") 
            testnumbers = numpy.array(testnumbers, dtype=float)
            margin = (testnumbers.max() - testnumbers.min())/20.0

            m1s = numpy.array(m1s, dtype=float)
            m2s = numpy.array(m2s, dtype=float)
            m1errs = numpy.array(m1errs, dtype=float)
            m2errs = numpy.array(m2errs, dtype=float)
            # put both graphs on the same units scale
            yrange = 2.0 * max((m1s.max() - m1s.min() + m1errs.max(), m2s.max() - m2s.min() + m2errs.max())) 
            # plot the data with error bars
            ax1.set_xlim(testnumbers.min() - margin, testnumbers.max() + margin)
            ax1.set_ylim((m1s.min() + m1s.max() - yrange)/2.0, (m1s.min() + m1s.max() + yrange)/2.0)
            ax1.set_xticks([])
            ax1.set_xticklabels(testlabels)
            ax1.errorbar(testnumbers, m1s, yerr=m1errs, marker = ".", markersize=10, linestyle='None', color='red')

            ax3.set_xlim(testnumbers.min() - margin, testnumbers.max() + margin)
            ax3.set_ylim((m2s.min() + m2s.max() - yrange)/2.0, (m2s.min() + m2s.max() + yrange)/2.0)
            ax3.set_xticks([])
            ax3.set_xticklabels(testlabels)
            ax3.errorbar(testnumbers, m2s, yerr=m2errs, marker = ".", markersize=10, linestyle='None')

            # plot the data with error bars
            b1s = numpy.array(b1s, dtype=float)
            b2s = numpy.array(b2s, dtype=float)
            b1errs = numpy.array(b1errs, dtype=float)
            b2errs = numpy.array(b2errs, dtype=float)
            # put both graphs on the same units scale
            yrange = 2.0 * max((b1s.max() - b1s.min() + b1errs.max(), b2s.max() - b2s.min() + b2errs.max())) 
            ax2.set_xlim(testnumbers.min() - margin, testnumbers.max() + margin)
            ax2.set_ylim((b1s.min() + b1s.max() - yrange)/2.0, (b1s.min() + b1s.max() + yrange)/2.0)
            ax2.set_xticks([])
            ax2.set_xticklabels(testlabels)
            ax2.errorbar(testnumbers, b1s, yerr=b1errs, marker = ".", markersize=10, linestyle='None', color='red')

            ax4.set_xlim(testnumbers.min() - margin, testnumbers.max() + margin)
            ax4.set_ylim((b2s.min() + b2s.max() - yrange)/2.0, (b2s.min() + b2s.max() + yrange)/2.0)
            ax4.set_xticklabels(testlabels)
            ax4.errorbar(testnumbers, b2s, yerr=b2errs, marker = ".", markersize=10, linestyle='None')
            plot.show()
