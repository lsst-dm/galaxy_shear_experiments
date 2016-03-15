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

import math
import numpy
import time
import lsst.daf.base
import lsst.pex.config
import lsst.pipe.base
import lsst.afw.table
import lsst.afw.image
import lsst.afw.math
import lsst.meas.algorithms
import lsst.meas.base
from lsst.meas.base.baseLib import MeasurementError
from lsst.meas.base.baseMeasurement import FATAL_EXCEPTIONS
from lsst.obs.great3.processBase import *

class ProcessShearTestConfig(ProcessBaseConfig):
    test = lsst.pex.config.Field(
        dtype=str,
        default=None,
        optional=True,
        doc="Name of the measurement test"
    )
    measPlugin = lsst.pex.config.Field(
        dtype=str,
        default=None,
        optional=True,
        doc="registered name of measurement plugin"
    )
    footprintSize = lsst.pex.config.Field(
        dtype=int,
        default=None,
        optional=True,
        doc="size of footprint, if a square is used"
    )
    stampSize = lsst.pex.config.Field(
        dtype=int,
        default=None,
        optional=True,
        doc="size of cutout, if reduction is desired"
    )
    psfSize = lsst.pex.config.Field(
        dtype=int,
        default=None,
        optional=True,
        doc="size of cutout, if reduction is desired"
    )
    noClobber = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="Delete existing output src table?"
    )

class ProcessShearTestTask(ProcessBaseTask):
    """This subclass obs.great3.ProcessBaseTask is used to call a measurement algorithm,
    assumed to be the only algorithm which is not a Centroid algorithm, using as its
    sources a GalSim catalog provided by the "epoch_catalog" for the dataId
    """

    ConfigClass = ProcessShearTestConfig

    _DefaultName = "processShearTest"

    def __init__(self, **kwds):
        ProcessBaseTask.__init__(self, **kwds)
        self.xKey = self.schema.addField("x", type = float,
            doc = "x position from epoch catalog")
        self.yKey = self.schema.addField("y", type = float,
            doc = "y position from epoch catalog")
        self.psfLibraryKey = self.schema.addField("psf_library",
            type = int, doc = "number of psf library")
        self.psfIndexKey = self.schema.addField("psf_index",
            type = int, doc = "index of psf within library")
        self.psfNumberKey = self.schema.addField("psf_number",
            type = int, doc = "number of psf_nnnn.fits")
        self.bulgeThetaKey = self.schema.addField("bulge_theta", type = float,
            doc = "bulge rotation angle")
        self.diskThetaKey = self.schema.addField("disk_theta", type = float,
            doc = "disk rotation angle")
        self.g1Key = self.schema.addField("g1", type = float,
            doc = "GalSim constant g1 value")
        self.g2Key = self.schema.addField("g2", type = float,
            doc = "GalSim constant g2 value")
        self.footprintCountKey = self.schema.addField("footprintCount", type = int,
            doc = "Number of footprint detected at 5 sigma")

    # Compare the identifiers in the catalog to see if two records have the same Psf
    def comparePsfIds(self, record1, record2):
        if not record1 is None and not record2 is None:
            if record1.get(self.psfLibraryKey) == record2.get(self.psfLibraryKey):
                if record1.get(self.psfIndexKey) == record2.get(self.psfIndexKey):
                    if record1.get(self.psfNumberKey) == record2.get(self.psfNumberKey):
                        return True
        return False

    def buildSourceCatalog(self, imageBBox, dataRef):
        """Build an empty source catalog, using the provided sim catalog's position to generate
        square Footprints and its ID to set that of the source catalog. Also Transfer information
        already known about the sources from the GalSim epoch catalog to an afwTable catalog.
        The Psf library and index in particular are needed.
        """
        sourceCat = lsst.afw.table.SourceCatalog(self.schema)
        measPlugin = self.config.measPlugin
        fieldNames = self.schema.extract(measPlugin + "*").keys()
        self.e1Key = None
        self.e2Key = None
        self.sigmaKey = None
        self.xxKey = None
        self.yyKey = None
        self.xyKey = None
        for name in fieldNames:
            if name.endswith("_e1"):
                self.e1Key = self.schema.find(name).getKey()
            if name.endswith("_e2"):
                self.e2Key = self.schema.find(name).getKey()
            if name.endswith("_sigma"):
                self.sigmaKey = self.schema.find(name).getKey()
            if name.endswith("_xx"):
                self.xxKey = self.schema.find(name).getKey()
            if name.endswith("_yy"):
                self.yyKey = self.schema.find(name).getKey()
            if name.endswith("_xy"):
                self.xyKey = self.schema.find(name).getKey()
        sourceCat.getTable().setMetadata(self.algMetadata)
        simCat = dataRef.get("epoch_catalog", immediate=True)
        xKey = simCat.schema.find('x').key
        yKey = simCat.schema.find('y').key
        psfIndexKey = simCat.schema.find('psf_index').key
        psfLibraryKey = simCat.schema.find('psf_library').key
        psfNumberKey = simCat.schema.find('psf_number').key
        bulgeBetaRadiansKey = simCat.schema.find('bulge_beta_radians').key
        diskBetaRadiansKey = simCat.schema.find('disk_beta_radians').key
        idKey = simCat.schema.find('ID').key
        g1Key = simCat.schema.find('g1').key
        g2Key = simCat.schema.find('g2').key
        nGals = imageBBox.getWidth() // self.config.galaxyStampSize
        assert nGals * self.config.galaxyStampSize == imageBBox.getWidth()
        n = imageBBox.getWidth() / nGals
        self.dims = lsst.afw.geom.Extent2I(n, n)
        offset = lsst.afw.geom.Extent2I(simCat[0][xKey], simCat[0][yKey])
        for simRecord in simCat:
            sourceRecord = sourceCat.addNew()
            sourceRecord.set(self.xKey, simRecord.get(xKey))
            sourceRecord.set(self.yKey, simRecord.get(yKey))
            sourceRecord.setId(simRecord.get(idKey))
            sourceRecord.set(self.psfIndexKey, simRecord.get(psfIndexKey))
            sourceRecord.set(self.psfLibraryKey, simRecord.get(psfLibraryKey))
            sourceRecord.set(self.psfNumberKey, simRecord.get(psfNumberKey))
            sourceRecord.set(self.bulgeThetaKey, simRecord.get(bulgeBetaRadiansKey))
            sourceRecord.set(self.diskThetaKey, simRecord.get(diskBetaRadiansKey))
            sourceRecord.set(self.g1Key, simRecord.get(g1Key))
            sourceRecord.set(self.g2Key, simRecord.get(g2Key))
            position = lsst.afw.geom.Point2I(simRecord.get(xKey), simRecord.get(yKey))
            # the default footprint is the entire rectangle produced by GalSim
            bbox = lsst.afw.geom.Box2I(position - offset, self.dims)
            footprint = lsst.afw.detection.Footprint(bbox, imageBBox)
            peakRecord = footprint.getPeaks().addNew()
            peakRecord.setFx(position.getX())
            peakRecord.setFy(position.getY())
            sourceRecord.setFootprint(footprint)
        return sourceCat

    def run(self, dataRef):
        """Run a measurement algorithm on the set of sources provided by
        a GalSim catalog.  Also run a centroid algorithm if one is provided
        This run method has a measurement task (which makes it easier to
        setup the plugins and schema, but it does not actually use the measurement
        task, instead calling the plugins directly to improve speed and remove
        noise replacement.
        """
        startTime = time.time()
        # The noClobber flag protects data which was already run from deletion
        if self.config.noClobber:
            if dataRef.datasetExists("src"):
                return
        exposure = self.buildExposure(dataRef)
        sourceCat = self.buildSourceCatalog(exposure.getBBox(lsst.afw.image.PARENT), dataRef)
        images_file = dataRef.get("image", immediate=True)
        #   the GalSim catalog will provide the number of the psf_library
        #   and the index of the hdu with the psf for each source, and the
        #   internal number (identifier) of the psf.
        count = 0
        e1_sum = 0.0
        e2_sum = 0.0
        weight_sum = 0.0
        if len(sourceCat) > 0:
            g1 = sourceCat[0].get(self.g1Key)
            g2 = sourceCat[0].get(self.g2Key)
            theta0 = 180.0*math.atan2(g2,g1)/2.0*math.pi
        lastRecord = None
        lastSPARecord = None
        measErrors = 0
        # Dictionary of keys for Psf Modeling Fields
        psfApproxKeys = sourceCat.getSchema().extract("*modelfit_ShapeletPsfApprox*").keys()
        for source in sourceCat:
            #   Locate the psf_library and get the indexed psf but only if the Psf has changed
            #   This will make psf loading more efficient
            if not self.comparePsfIds(lastRecord, source):
                dataRef.dataId["psf_index"] = source.get(self.psfIndexKey)
                dataRef.dataId["psf_library"] = source.get(self.psfLibraryKey)
                dataRef.dataId["psf_number"] = source.get(self.psfNumberKey)
                if dataRef.datasetExists("psf_file"):
                    psf_file = dataRef.get("psf_file", immediate=True)
                else:
                    print "loading Library %d, for psfnumber %d"%(source.get("psf_library"),
                        source.get("psf_number"))
                    psf_file = dataRef.get("psf_library", immediate=True)
                data = psf_file.getArray().astype(numpy.float64)
                # Resize the Psf image if requested
                if not self.config.psfSize is None:
                    if data.shape[0] < self.config.psfSize:
                        raise Exception("psfSize %d requested is larger than actual psf"%self.config.psfSize)
                    clip = (data.shape[0] - self.config.psfSize)/2
                    data = data[clip:data.shape[1]-clip, clip:data.shape[0]-clip]
                    if lastRecord is None:
                        self.log.warn("Psf clipped to %d x %d"%data.shape)
                kernel = lsst.afw.math.FixedKernel(lsst.afw.image.ImageD(data))
                psf = lsst.meas.algorithms.KernelPsf(kernel)
            lastRecord = source

            # Create a bounding box around the galaxy image of self.dims
            x = int(source.getFootprint().getCentroid().getX()+1) - (self.dims.getX()/2)
            y = int(source.getFootprint().getCentroid().getY()+1) - (self.dims.getY()/2)
            bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(x,y), self.dims)
            exp = lsst.afw.image.ExposureF(exposure, bbox)
            exp.setPsf(psf)
            # Add a fake WCS and calib, just to keep all the meas_base algorithms happy
            CD = numpy.array([[5.55E-5, 0.0], [0.0, 5.55E-5]])
            crpix = lsst.afw.geom.Point2D(0.0,0.0)
            crval = lsst.afw.geom.Point2D(0.0,0.0)
            exp.setWcs(lsst.afw.image.Wcs(crval, crpix, CD))
            calib = lsst.afw.image.Calib()
            calib.setFluxMag0((3531360874589.57, 21671681149.139))
            exp.setCalib(calib)

            #  Add a real footprint unless a rectangular footprint has been requested
            if self.config.footprintSize is None:
                try:
                    task = lsst.meas.algorithms.SourceDetectionTask()
                    task.log.setThreshold(task.log.WARN)
                    footprints = task.detectFootprints(exp,
                             sigma=4.0).positive.getFootprints()
                    if len(footprints) == 1:
                        source.setFootprint(footprints[0])
                    # find the one closest to the original centroid
                    else:
                        continue
                        ds = 1000000  # absurdly large distance
                        for i in range(len(footprints)):
                            fpi = footprints[i]
                            dsi = source.getFootprint().getCentroid().distanceSquared(fp.getCentroid())
                            if dsi < ds:
                                fp = fpi
                                ds = dsi
                        if math.sqrt(ds) > 2:      #multiple detection too far away
                            continue
                        source.setFootprint(fp)
                    source.set(self.footprintCountKey, len(footprints))
                except:
                    source.set(self.footprintCountKey, -1)
            #  else make the footprint a square of the requested size
            else:
                bbox = source.getFootprint().getBBox()
                shrinkAmount = (bbox.getDimensions().getX() - self.config.footprintSize)/2
                bbox.grow(-shrinkAmount)
                source.getFootprint().clipTo(bbox)

            #  shrink the bounds of the exposure if requested:
            if not self.config.stampSize is None:
                x = int(source.getFootprint().getCentroid().getX()+1) - (self.config.stampSize/2)
                y = int(source.getFootprint().getCentroid().getY()+1) - (self.config.stampSize/2)
                bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(x,y),
                           lsst.afw.geom.Extent2I(self.config.stampSize, self.config.stampSize))
                exp = lsst.afw.image.ExposureF(exp, bbox)

            #  Now do the measurements, calling the measure algorithms individually to increase speed
            sigma = None
            if count % 1000 == 0:
                self.log.warn("Measuring subfield %d, count = %d"%(dataRef.dataId["subfield"], count))
            count = count + 1
            try:
                for plugin in self.measurement.plugins.keys():
                    # optmize by not calling SPA if Psf same as lastSPARecord
                    if plugin.find("modelfit_ShapeletPsfApprox") >= 0:
                        if self.comparePsfIds(source, lastSPARecord):
                            self.log.warn("Re-using the last SPA value for %d"%count)
                            for key in psfApproxKeys:
                                source.set(key, lastSPARecord.get(key))
                            continue
                    self.measurement.plugins[plugin].measure(source, exp)

                    # Save the last successfule SPA record
                    if plugin.find("modelfit_ShapeletPsfApprox") >= 0:
                        lastSPARecord = source

            except FATAL_EXCEPTIONS:
                raise
            except MeasurementError as error:
                measErrors = measErrors + 1
                self.measurement.plugins[plugin].fail(source, error)
            except Exception as error:
                self.log.warn("Error in %s.measure on record %s: %s"
                          % (self.measurement.plugins[plugin].name, source.getId(), error))
            if not self.e1Key is None:
                e1 = source.get(self.e1Key)
                e2 = source.get(self.e2Key)
            else:
                xx = source.get(self.xxKey)
                yy = source.get(self.yyKey)
                xy = source.get(self.xyKey)
                e1 = (xx - yy) / (xx + yy)
                e2 = 2 * xy / (xx + yy)
            e = math.sqrt(e1 * e1 + e2 * e2)

            #   Calculate the rotation in the range [-pi, pi], just to check the pairs
            theta = 180.0*math.atan2(e2,e1)/(2.0*math.pi)
            while theta < 0:
                theta = 360 + theta
            while theta > 90.0:
                theta = theta-90.0

        dataRef.put(sourceCat, "src")
        self.log.warn("Elapsed time = %f secs."%(time.time()-startTime))
        self.log.warn("measError count = %d"%measErrors)

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument(name="--id", datasetType="image", level="image",
                               help="data ID, e.g. --id subfield=0")
        return parser

    def writeConfig(self, butler, clobber=False, doBackup=False):
        pass

    def writeSchemas(self, butler, clobber=False, doBackup=False):
        pass

    def writeMetadata(self, dataRef, doBackup=False):
        pass

    def writeEupsVersions(self, butler, clobber=False, doBackup=False):
        pass
