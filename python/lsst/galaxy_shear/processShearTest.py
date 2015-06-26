#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 20082014 LSST Corporation.
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

import numpy
import lsst.daf.base
import lsst.pex.config
import lsst.pipe.base
import lsst.afw.table
import lsst.afw.image
import lsst.afw.math
import lsst.meas.algorithms

from lsst.obs.great3.processBase import *

class ProcessShearTestConfig(ProcessBaseConfig):
    pass
"""
   This subclass obs.great3.ProcessBaseTask is used to call a measurement algorithm,
   assumed to be the only algorithm which is not a Centroid algorithm, using as its
   sources a GalSim catalog provided by the "epoch_catalog" for the dataId
"""
class ProcessShearTestTask(ProcessBaseTask):

    ConfigClass = ProcessShearTestConfig

    _DefaultName = "processShearTest"

    def __init__(self, **kwds):
        ProcessBaseTask.__init__(self, **kwds)
        self.schema.addField("psf_index", type = int, doc = "index of psf within library")
        self.schema.addField("psf_library", type = int, doc = "number of psf_libary_nn.fits")
        self.centroidPlugin = None
        self.measPlugin = None
        for plugin in self.measurement.plugins.keys():
            if plugin.find("Centroid") > 0:
                self.centroidPlugin = self.measurement.plugins[plugin]
            else:
                self.measPlugin = self.measurement.plugins[plugin]

    #   Transfer information already know about the sources from the GalSim epoch catalog
    #   to an afwTable catalog.  The Psf library and index in particular are needed.
    def buildSourceCatalog(self, imageBBox, dataRef):
        """Build an empty source catalog, using the provided sim catalog's position to generate
        square Footprints and its ID to set that of the source catalog.
        """
        sourceCat = lsst.afw.table.SourceCatalog(self.schema)
        sourceCat.getTable().setMetadata(self.algMetadata)
        simCat = dataRef.get("epoch_catalog", immediate=True)
        xKey = simCat.schema.find('x').key
        yKey = simCat.schema.find('y').key
        idKey = simCat.schema.find('ID').key
        nGals = imageBBox.getWidth() / self.config.galaxyStampSize
        assert nGals * self.config.galaxyStampSize == imageBBox.getWidth()
        n = imageBBox.getWidth() / nGals
        self.dims = lsst.afw.geom.Extent2I(n, n)
        offset = lsst.afw.geom.Extent2I(simCat[0][xKey], simCat[0][yKey])
        for simRecord in simCat:
            sourceRecord = sourceCat.addNew()
            sourceRecord.setId(simRecord.get(idKey))
            sourceRecord.set("psf_index", simRecord.get("psf_index"))
            sourceRecord.set("psf_library", simRecord.get("psf_library"))
            position = lsst.afw.geom.Point2I(simRecord.get(xKey), simRecord.get(yKey))
            bbox = lsst.afw.geom.Box2I(position - offset, self.dims)
            footprint = lsst.afw.detection.Footprint(bbox, imageBBox)
            peakRecord = footprint.getPeaks().addNew()
            peakRecord.setFx(position.getX())
            peakRecord.setFy(position.getY())
            sourceRecord.setFootprint(footprint)
        return sourceCat

    #   Run a measurement algorithm on the set of sources provided by
    #   a GalSim catalog.  Also run a centroid algorithm if one is provided
    def run(self, dataRef):
        exposure = self.buildExposure(dataRef)
        sourceCat = self.buildSourceCatalog(exposure.getBBox(lsst.afw.image.PARENT), dataRef)
        images_file = dataRef.get("image", immediate=True)
        #   the GalSim catalog will provide the number of the psf_library
        #   and the index of the hdu with the psf for each source.
        for source in sourceCat:
            #   Locate the psb_library and get the indexed psf
            dataRef.dataId["psf_index"] = source.get("psf_index")
            dataRef.dataId["psf_library"] = source.get("psf_library")
            psf_file = dataRef.get("psf_library", immediate=True)
            data = psf_file.getArray().astype(numpy.float64)
            kernel = lsst.afw.math.FixedKernel(lsst.afw.image.ImageD(data))
            psf = lsst.meas.algorithms.KernelPsf(kernel)
            # Create a bounding box around the galaxy image of self.dims
            x = int(source.getFootprint().getCentroid().getX()+1) - (self.dims.getX()/2)
            y = int(source.getFootprint().getCentroid().getY()+1) - (self.dims.getY()/2)
            bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(x,y), self.dims)
            exp = lsst.afw.image.ExposureF(exposure, bbox)
            exp.setPsf(psf)
            #  Now do the measurements, calling the measure algorithms to increase speed
            if self.centroidPlugin:
                try:
                    self.centroidPlugin.measure(source, exp)
                except:
                    self.centroidPlugin.fail(source)
                    continue
            if self.measPlugin:
                try:
                    self.measPlugin.measure(source, exp)
                except:
                    self.measPlugin.fail(source)
        dataRef.put(sourceCat, "src")

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument(name="--id", datasetType="image", level="image",
                               help="data ID, e.g. --id subfield=0")
        return parser

    def writeConfig(self, butler, clobber=False):
        pass

    def writeSchemas(self, butler, clobber=False):
        pass

    def writeMetadata(self, dataRef):
        pass

    def writeEupsVersions(self, butler, clobber=False):
        pass
