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

class ProcessShearTestTask(ProcessBaseTask):

    ConfigClass = ProcessShearTestConfig

    _DefaultName = "processShearTest"

    def __init__(self, **kwds):
        ProcessBaseTask.__init__(self, **kwds)
        self.schema.addField("psf_index", type = long, doc = "used for psf indexing")
        self.schema.addField("psf_file_number", type = int, doc = "used for psf indexing")

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
        psfKey = simCat.schema.find('psf_index').key
        psfFileKey = simCat.schema.find('psf_file_number').key
        nGals = imageBBox.getWidth() / self.config.galaxyStampSize
        assert nGals * self.config.galaxyStampSize == imageBBox.getWidth()
        n = imageBBox.getWidth() / nGals
        dims = lsst.afw.geom.Extent2I(n, n)
        offset = lsst.afw.geom.Extent2I(simCat[0][xKey], simCat[0][yKey])
        for simRecord in simCat:
            sourceRecord = sourceCat.addNew()
            sourceRecord.setId(simRecord.get(idKey))
            sourceRecord.set("psf_index", simRecord.get(psfKey))
            sourceRecord.set("psf_file_number", simRecord.get(psfFileKey))
            position = lsst.afw.geom.Point2I(simRecord.get(xKey), simRecord.get(yKey))
            bbox = lsst.afw.geom.Box2I(position - offset, dims)
            footprint = lsst.afw.detection.Footprint(bbox, imageBBox)
            peakRecord = footprint.getPeaks().addNew()
            peakRecord.setFx(position.getX())
            peakRecord.setFy(position.getY())
            sourceRecord.setFootprint(footprint)
        return sourceCat

    def run(self, dataRef):

        exposure = self.buildExposure(dataRef)
        sourceCat = self.buildSourceCatalog(exposure.getBBox(lsst.afw.image.PARENT), dataRef)
        images_file = dataRef.get("image", immediate=True)

        runCat = lsst.afw.table.SourceCatalog(self.schema)
        runCat.addNew()
        for source in sourceCat:
            runCat[0] = source
            dataRef.dataId["psf_index"] = source.get("psf_index")
            dataRef.dataId["psf_file_number"] = source.get("psf_file_number")
            psf_file = dataRef.get("psf_files", immediate=True)
            data = psf_file.getArray().astype(numpy.float64) 
            kernel = lsst.afw.math.FixedKernel(lsst.afw.image.ImageD(data))
            psf = lsst.meas.algorithms.KernelPsf(kernel)
            exposure.setPsf(psf)
            self.measurement.run(exposure, runCat)
            print source.get("base_PsfFlux_flux"), source.get("base_PsfFlux_fluxSigma")
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
