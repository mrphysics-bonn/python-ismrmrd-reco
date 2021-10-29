""" Helper functions for writing protocol data via Ismrmrd or (deprecated) h5py.
"""

def create_hdr(hdr, params):
    import ismrmrd

    # initialize some optional parameters
    if "dwelltime" not in params:
        params["dwelltime"] = 0
    if "slices" not in params:
        params["slices"] = 1
    if "nsegments" not in params:
        params["nsegments"] = 1
    if "npartitions" not in params:
        params["npartitions"] = 1
    if "nintl" not in params:
        params["nintl"] = 1
    if "ncontrast" not in params:
        params["ncontrast"] = 1
    if "avg" not in params:
        params["avg"] = 1
    if "traj_delay" not in params:
        params["traj_delay"] = 0
    if "t_min" not in params:
        params["t_min"] = 0
    if "os_region" not in params:
        params["os_region"] = 1
    if "redfac" not in params:
        params['redfac'] = 1
    if "sms_factor" not in params:
        params['sms_factor'] = 1
    if "trajtype" not in params:
        raise ValueError("Specify trajectory type: spiral or cartesian")
        
    # experimental conditions
    exp = ismrmrd.xsd.experimentalConditionsType()
    exp.H1resonanceFrequency_Hz = 298060000 # 7T
    hdr.experimentalConditions = exp

    # user parameters
    dtime = ismrmrd.xsd.userParameterDoubleType()
    dtime.name = 'dwellTime_us'
    dtime.value = params["dwelltime"] * 1e6 
    tdel = ismrmrd.xsd.userParameterDoubleType()
    tdel.name = 'traj_delay'
    tdel.value = params["traj_delay"] # this field can be used for any additional trajectory delay [s]
    nseg = ismrmrd.xsd.userParameterDoubleType()
    nseg.name = 'nsegments'
    nseg.value = params["nsegments"] # deprecated, but kept for compatibility
    tmin = ismrmrd.xsd.userParameterDoubleType()
    tmin.name = 't_min'
    tmin.value = params["t_min"] # starting time of t-vector for B0-correction
    os_region = ismrmrd.xsd.userParameterDoubleType()
    os_region.name = 'os_region'
    os_region.value = params["os_region"] # factor for the spiral oversampling region - phase maps will be reconstructed up to this factor of kspace

    hdr.userParameters = ismrmrd.xsd.userParametersType()
    hdr.userParameters.userParameterDouble.append(dtime)
    hdr.userParameters.userParameterDouble.append(tdel)
    hdr.userParameters.userParameterDouble.append(nseg)
    hdr.userParameters.userParameterDouble.append(tmin)
    hdr.userParameters.userParameterDouble.append(os_region)

    # encoding
    encoding = ismrmrd.xsd.encodingType()
    if params["trajtype"] == "spiral":
        encoding.trajectory = ismrmrd.xsd.trajectoryType('spiral')
    if params["trajtype"] == "cartesian":
        encoding.trajectory = ismrmrd.xsd.trajectoryType('cartesian')

    # encoded and recon spaces
    efov = ismrmrd.xsd.fieldOfViewMm()
    efov.x = params["fov"]
    efov.y = params["fov"]
    efov.z = 1e3 * params["slice_res"] * params["slices"] * params["npartitions"]
    rfov = ismrmrd.xsd.fieldOfViewMm()
    rfov.x = params["fov"]
    rfov.y = params["fov"]
    rfov.z = 1e3 * params["slice_res"] * params["slices"] * params["npartitions"]
    
    ematrix = ismrmrd.xsd.matrixSizeType()
    ematrix.x = int(params["fov"]/params["res"]+0.5)
    ematrix.y = int(params["fov"]/params["res"]+0.5)
    ematrix.z = params["npartitions"]
    rmatrix = ismrmrd.xsd.matrixSizeType()
    rmatrix.x = int(params["fov"]/params["res"]+0.5)
    rmatrix.y = int(params["fov"]/params["res"]+0.5)
    rmatrix.z = params["npartitions"]
    
    espace = ismrmrd.xsd.encodingSpaceType() # encoded space incl oversampling
    espace.matrixSize = ematrix
    espace.fieldOfView_mm = efov
    rspace = ismrmrd.xsd.encodingSpaceType()
    rspace.matrixSize = rmatrix
    rspace.fieldOfView_mm = rfov
    
    # set encoded and recon spaces
    encoding.encodedSpace = espace
    encoding.reconSpace = rspace
    
    # encoding limits
    limits = ismrmrd.xsd.encodingLimitsType()
    limits.slice = ismrmrd.xsd.limitType()
    limits.slice.minimum = 0
    limits.slice.maximum = params["slices"] - 1
    limits.slice.center = int(params["slices"]/2)

    limits.kspace_encoding_step_1 = ismrmrd.xsd.limitType()
    limits.kspace_encoding_step_1.minimum = 0
    limits.kspace_encoding_step_1.maximum = params["nintl"] - 1
    limits.kspace_encoding_step_1.center = 0

    limits.average = ismrmrd.xsd.limitType()
    limits.average.minimum = 0
    limits.average.maximum = params["avg"] - 1
    limits.average.center = 0

    limits.phase = ismrmrd.xsd.limitType()
    limits.phase.minimum = 0
    limits.phase.maximum =  0
    limits.phase.center = 0

    limits.contrast = ismrmrd.xsd.limitType()
    limits.contrast.minimum = 0
    limits.contrast.maximum = params["ncontrast"] - 1
    limits.contrast.center = int(params["ncontrast"]/2)

    limits.segment = ismrmrd.xsd.limitType()
    limits.segment.minimum = 0
    limits.segment.maximum = params["nsegments"] - 1 # number of ADC segments
    limits.segment.center = 0

    pi = ismrmrd.xsd.parallelImagingType()
    acc_fac = ismrmrd.xsd.accelerationFactorType()
    acc_fac.kspace_encoding_step_1 = params['redfac']
    acc_fac.kspace_encoding_step_2 = params['sms_factor'] # use this field for multiband factor
    pi.accelerationFactor = acc_fac
    encoding.parallelImaging = pi

    encoding.encodingLimits = limits
    
    hdr.encoding.append(encoding)

################### DEPRECATED ###############################

def add_empty_acq(acqs, acq_ctr, rotmat = None):
    """ Creates an empty acquisition in the acquisitions group of an hdf5 protocol

    acqs: acquisition group of the hdf5 protocol
    acq_ctr: number of the acquisition to be inserted
    rotmat: optionally add a rotation matrix
    """
    import h5py
    import numpy as np

    acq = acqs.create_group('%d' %acq_ctr)
    acq.attrs['trajectory_dimensions'] = 2 # set 2D acq as default

    # rotation matrix
    if rotmat is None:
        acq.create_dataset('phase_dir', data=np.array([1,0,0]), dtype = np.float32)
        acq.create_dataset('read_dir', data=np.array([0,1,0]), dtype = np.float32)
        acq.create_dataset('slice_dir', data=np.array([0,0,1]), dtype = np.float32)
    else:
        acq.create_dataset('phase_dir', data=rotmat[:,0], dtype = rotmat.dtype)
        acq.create_dataset('read_dir', data=rotmat[:,1], dtype = rotmat.dtype)
        acq.create_dataset('slice_dir', data=rotmat[:,2], dtype = rotmat.dtype)

    # flags - WIP: this is not the complete list of flags
    flags = acq.create_group('flags')
    flags.attrs['ACQ_IS_NOISE_MEASUREMENT'] = False
    flags.attrs['ACQ_IS_PHASECORR_DATA'] = False
    flags.attrs['ACQ_IS_DUMMYSCAN_DATA'] = False # sync scans here
    flags.attrs['ACQ_IS_PARALLEL_CALIBRATION'] = False
    flags.attrs['ACQ_LAST_IN_SLICE'] = False
    flags.attrs['ACQ_LAST_IN_REPETITION'] = False

    # encoding counters
    idx = acq.create_group('idx')
    idx.attrs['kspace_encode_step_1'] = 0 # phase enc (intl)
    idx.attrs['kspace_encode_step_2'] = 0 # partition
    idx.attrs['slice'] = 0
    idx.attrs['contrast'] = 0
    idx.attrs['phase'] = 0
    idx.attrs['average'] = 0
    idx.attrs['repetition'] = 0
    idx.attrs['set'] = 0
    idx.attrs['segment'] = 0 # segment number for segmented acquisition


    #### Old h5py header
    # ## set up the protocol file
    # prot = h5py.File("../Protocols/"+date+"/{}.h5".format(filename), 'w')
    # hdr = prot.create_group('hdr')
    # acqs = prot.create_group('acquisitions')

    # ## write header
    # # user parameters
    # user_prm = hdr.create_group('userParameters')
    # user_prm_dbl = user_prm.create_group('userParameterDouble')
    # user_prm_dbl.attrs['dwellTime_us'] = dwelltime * 1e6
    # user_prm_dbl.attrs['traj_delay'] = spiral_delay # this field can be used for any additional trajectory delay [s]

    # # encodings
    # enc = hdr.create_group('encoding')
    # enc1 = enc.create_group('0')
    # enc1.attrs['trajectory'] = 'spiral'

    # enc1_encsp = enc1.create_group('encodedSpace') # encoded space incl oversampling
    # enc1_encsp_matr = enc1_encsp.create_group('matrixSize')
    # enc1_encsp_matr.attrs['x'] = int(fov/res+0.5)
    # enc1_encsp_matr.attrs['y'] = int(fov/res+0.5)
    # enc1_encsp_matr.attrs['z'] = 1 # only 2D atm

    # enc1_encsp_fov = enc1_encsp.create_group('fieldOfView_mm')
    # enc1_encsp_fov.attrs['x'] = fov
    # enc1_encsp_fov.attrs['y'] = fov
    # enc1_encsp_fov.attrs['z'] = 1e3 * slice_res * slices

    # enc1_reconsp = enc1.create_group('reconSpace')
    # enc1_reconsp_matr = enc1_reconsp.create_group('matrixSize')
    # enc1_reconsp_matr.attrs['x'] = int(fov/res+0.5)
    # enc1_reconsp_matr.attrs['y'] = int(fov/res+0.5)
    # enc1_reconsp_matr.attrs['z'] = 1

    # enc1_reconsp_fov = enc1_reconsp.create_group('fieldOfView_mm')
    # enc1_reconsp_fov.attrs['x'] = fov
    # enc1_reconsp_fov.attrs['y'] = fov
    # enc1_reconsp_fov.attrs['z'] = 1e3 * slice_res * slices
    