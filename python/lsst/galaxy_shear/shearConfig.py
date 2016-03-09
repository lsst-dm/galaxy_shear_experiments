import lsst.pex.config as pexConfig

class RunShearConfig(pexConfig.Config):

    # values used for the great3/galsim runs
    exp_type = pexConfig.Field(dtype=str, optional=True,
        doc="experiment type: control, real_galaxy, etc.", default="control")
    psf_size = pexConfig.Field(dtype=int, optional=True,
        doc="dimension of psf images", default=33)
    n_subfields = pexConfig.Field(dtype=int, optional=True,
        doc="number of galsim subfields", default=512)
    n_epochs = pexConfig.Field(dtype=int, optional=True,
        doc="number of galsim epochs", default=1)
    ndims = pexConfig.Field(dtype=int, optional=True,
        doc="size of grid for galsim image files", default=32)
    galaxy_stamp_size = pexConfig.Field(dtype=int, optional=True,
        doc="size of each galaxy image", default=96)
    fov = pexConfig.Field(dtype=float, optional=True,
        doc="galsim field of view", default=1.500000)
    gal_dir = pexConfig.Field(dtype=str, optional=True,
        doc="directory containing galsim real galaxy sample",
        default="/sandbox/lsstshared/pgee/mylsst9/COSMOS_23.5_training_sample")
    psf_lib_dir = pexConfig.Field(dtype=str, optional=True,
        doc="directory containing psf libraries",
        default="/sandbox/lsstshared/pgee/mylsst9/psflibs")
    psf_dir = pexConfig.Field(dtype=str, optional=True,
        doc="directory", default=None)

    # fields used in measurement algorithm
    shape_field = pexConfig.Field(dtype=str, optional=True,
        doc="base name of fields used for ellipticity measurements",
        default="ext_shapeHSM_HsmShapeRegauss")
    shape_type = pexConfig.Field(dtype=str, optional=True,
        doc="type of shape measurement: moments or ellipticity", default="ellipticity")
    sigma_field = pexConfig.Field(dtype=str, optional=True,
        doc="field used for measurement weighting",
        default="ext_shapeHSM_HsmShapeRegauss_sigma")
    sigma_signal_field = pexConfig.Field(dtype=str, optional=True,
       doc="if defined, sigma_field is for this field, not the shape_field", default=None)

    # values for setting of filter, seeing, and shear.  Also used to name run subdirectory
    filter = pexConfig.Field(dtype=int, optional=True,
        doc="", default=None)
    seeing = pexConfig.Field(dtype=float, optional=True,
        doc="", default=None)
    shear_value = pexConfig.Field(dtype=float, optional=True,
        doc="", default=None)
    shear_angle = pexConfig.Field(dtype=float, optional=True,
        doc="", default=None)


