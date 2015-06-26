from lsst.galaxy_shear.processShearTest import ProcessShearTestTask
import lsst.meas.extensions.shapeHSM
root.measurement.algorithms.names = ["base_GaussianCentroid", "ext_shapeHSM_HsmShapeRegauss"]
root.measurement.slots.centroid = 'base_GaussianCentroid'
root.measurement.slots.shape = None
root.measurement.slots.apFlux = None
root.measurement.slots.instFlux = None
root.measurement.slots.psfFlux = None

