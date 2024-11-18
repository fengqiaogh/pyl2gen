import numpy as np


class L1_input_t:
    def __init__(self):

        # file reading control params
        self.calfile = 0
        self.xcal_file = 0
        self.btfile = 0
        self.cld_msk_file = ""
        self.viirscalparfile = ""
        self.pversion = "Unspecified"
        self.input_parms = 0
        self.input_files = 0

        self.rad_opt = 1  # radcor switch for MERIS smile correction
        self.geom_per_band = 0  # 0 - use nominal geometry sen/sol_a/a
        # 1 - use band-specific values for instruments
        # that have band specific geometry
        self.xcal_nwave = 0  # number of wavelengths to which xcal applied
        self.xcal_opt = None  # xcal option per band
        self.xcal_wave = None  # sensor wavelengths to which xcal applied
        self.resolution = -1  # process at this nadir pixel res (meters)
        # 250, 500, 1000, -1=native (modis only)
        self.newavhrrcal = 0  # new avhrr calibration equation
        self.sl_pixl = -1  # seawifs straylight pixel limit
        self.sl_frac = 0.25  # seawifs straylight Ltyp fraction
        self.ch22detcor = np.ones(10)  # channel 22 detector corrections
        self.ch23detcor = np.ones(10)  # channel 23 detector corrections
        self.cirrus_thresh = np.ones(2)  # cirrus reflectance thresholds
        self.albedo = -1.0  # cloud reflectance threshold
        self.cloud_wave = 865.0  # cloud test wavelength
        self.cloud_eps = -1.0  # cloud reflectance ratio
        self.glint = 0.005  # glint threshold
        self.sunzen = 75.0  # solar zenith angle threshold
        self.satzen = 60.0  # sensor zenith angle threshold
        self.hipol = 0.50  # high polarization threshold
        self.gain = None  # Vicarious calibration gain
        self.offset = None  # Vicarious calibration offset

        self.spixl = 1  # starting pixel no. of the input (1-rel)
        self.epixl = -1  # ending pixel no. of the input (1-rel)
        self.dpixl = 1  # pixel subsampling increment
        self.sline = 1  # starting line no. of the input (1-rel)
        self.eline = -1  # ending line no. of the input (1-rel)
        self.dline = 1  # line subsampling increment

        self.outband_opt = 99  # 1=apply seawifs out-of-band correction
        self.evalmask = 0
        self.landmask = 1  # 0=off, 1=on
        self.bathmask = 0  # 0=off, 1=on
        self.cloudmask = 0  # 0=off, 1=on
        self.glintmask = 0  # 0=off, 1=on
        self.sunzenmask = 0  # 0=off, 1=on
        self.satzenmask = 0  # 0=off, 1=on
        self.hiltmask = 0  # 0=off, 1=on
        self.stlightmask = 0  # 0=off, 1=on
