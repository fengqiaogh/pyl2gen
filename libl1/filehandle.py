import numpy as np

FORWARD = 0
INVERSE_ZERO = 1
INVERSE_NLW = 2
INVERSE_LW = 3

READ = 0
WRITE = 1

L1_PRODSTRLEN = 2048
L1_MAXPROD = 1000
L1_NFLAGS = 32

NBANDSIR = 8


class Filehandle:
    def __init__(self):
        self.name = ""
        self.format = -1
        self.sensorID = -1
        self.subsensorID = -1
        self.spatialResolution = ""
        self.length = 0
        self.spix = 0
        self.epix = -1
        self.npix = 0
        self.nscan = 0
        self.nbands = 0
        self.nbandsir = 0
        self.nlvl = 42  # fixed 42 GMAO FP-IT levels now
        self.n_refl_loc = 10  # reflectance location 3rd dim count
        self.n_cloud_phase = 2  # cloud phases we retrieve
        self.bindx = None  # index to closest seawifs band
        self.ndets = 1
        self.mode = READ
        self.l2prod = ""  # list of L2 products to be included
        self.def_l2prod = ""  # list of default L2 products
        self.sd_id = 0  # hdf file id for the opened output file
        self.tot_prod = 0  # total # of L2 products to be created
        self.l2_prod_names = ""
        self.prodptr = None  # array of product structures
        self.productInfos = None  # pointers to product info structures

        self.geofile = None
        self.gmpfile = None
        self.orbit_node_lon = -999.0
        self.orbit_number = 0
        self.node_crossing_time = 0
        self.flag_cnt = 0
        self.terrain_corrected = 0
        self.sv_with_moon = 0
        self.grp_id = np.full(8, -1)  # netCDF group IDs

        self.iwave = None
        self.fwave = None
        self.fwhm = None
        self.Fobar = None
        self.Fonom = None
        self.Tau_r = None
        self.k_oz = None
        self.k_no2 = None
        self.aw = None
        self.bbw = None

        self.private_data = None

        self.sensorInfo = None

    def open(self):
        print(self.nbands)
