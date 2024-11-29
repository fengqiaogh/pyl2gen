import numpy as np
from oel_util.libtimeutils.unix2yds import unix2yds
from l2_struc import L2str
from oel_hdf4.filetype.filetype import FileType
import h5netcdf
from version import VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH, GITSHA
from l12_parms import PROGRAM
from libl1.l1 import L1str
from libl1.goci import GOCIL1

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
        self.format = FileType.FT_L2NCDF
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
        self.l2_prod_names = []
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

        # netCDF group IDs
        self.grp_id = {i: None for i in range(8)}

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

    def readl1(self):
        # l1rec.l1file = self
        match self.format:
            case FileType.FT_GOCIL1B:
                gocil1 = GOCIL1()
                gocil1.open(self)
            case _:
                print(f"readl1 - Unknown L1 input file format specifier: {self.format}")

    def writel2(self, l2rec: L2str):
        year, day, sec = unix2yds(l2rec.l1rec.scantime)
        msec = sec * 1e3
        print(f"Opening: {self.name}")
        match self.format:
            case FileType.FT_L2NCDF:
                calibrationDataStr = "calibration_data"
                equatorCrossingLonStr = "equatorCrossingLongitude"
                historyStr = "history"
                inputFilesStr = "input_sources"
                maskNamesStr = "mask_names"
                orbitNumberStr = "orbit_number"
                processingVersionStr = "processing_version"
                productNameStr = "product_name"
                softwareNameStr = "software_name"
                softwareVersionStr = "software_version"
                titleStr = "title"
                numberOfScanLinesStr = "number_of_lines"
                pixelsPerScanLineStr = "pixels_per_line"
                bandsPerPixelStr = "bands_per_pixel"
                profLvlPerPixelStr = "profile_levels_per_pixel"
                n_refl_loc_str = "number_of_reflectance_location_values"
                n_cloud_phase_str = "number_of_cloud_phases"
                totalBandNumberStr = "number_of_bands"
                bandNumberStr = "number_of_reflective_bands"
                wavelength_3d_str = "wavelength_3d"
        self.l2_prod_names.append("l2_flags")
        print(f"The following products will be included in {self.name}.")
        for idx, prod in enumerate(self.l2_prod_names):
            print(idx, prod)
        # Write out some global attributes
        title = f"{sensorId2SensorName} Level-2 Data {input.suite}"
        soft_id = f"{VERSION_MAJOR}.{VERSION_MINOR}.{VERSION_PATCH}.{GITSHA}"
        with h5netcdf.File(self.name, "w") as f:
            sensor_band_parameters = f.create_group("sensor_band_parameters")
            sensor_band_parameters.create_dataset("wavelength", data=wavelength)
            sensor_band_parameters.create_dataset("vcal_gain", data=vcal_gain)
            sensor_band_parameters.create_dataset("vcal_offset", data=vcal_offset)
            sensor_band_parameters.create_dataset("F0", data=F0)
            sensor_band_parameters.create_dataset("aw", data=aw)
            sensor_band_parameters.create_dataset("bbw", data=bbw)
            sensor_band_parameters.create_dataset("k_oz", data=k_oz)
            sensor_band_parameters.create_dataset("k_no2", data=k_no2)
            sensor_band_parameters.create_dataset("Tau_r", data=Tau_r)

            scan_line_attributes = f.create_group("scan_line_attributes")
            scan_line_attributes.create_dataset("year", data=year)
            scan_line_attributes.create_dataset("day", data=day)
            scan_line_attributes.create_dataset("msec", data=msec)
            scan_line_attributes.create_dataset("time", data=time)
            scan_line_attributes.create_dataset("detnum", data=detnum)
            scan_line_attributes.create_dataset("mside", data=mside)
            scan_line_attributes.create_dataset("slon", data=slon)
            scan_line_attributes.create_dataset("clon", data=clon)
            scan_line_attributes.create_dataset("elon", data=elon)
            scan_line_attributes.create_dataset("slat", data=slat)
            scan_line_attributes.create_dataset("clat", data=clat)
            scan_line_attributes.create_dataset("elat", data=elat)
            scan_line_attributes.create_dataset("csol_z", data=csol_z)

            geophysical_data = f.create_group("geophysical_data")
            for idx, prod in enumerate(self.l2_prod_names):
                geophysical_data.create_dataset(prod, data=self.l2_data[idx])

            navigation_data = f.create_group("navigation_data")
            navigation_data.create_dataset("longitude", data=longitude)
            navigation_data.create_dataset("latitude", data=latitude)
            navigation_data.create_dataset("tilt", data=tilt)

            # Processing Control attibutes
            processing_control = f.create_group("processing_control")

            # TODO 这里应该是个属性
            processing_control.create_dataset(softwareNameStr, data=PROGRAM)
            processing_control.create_dataset(softwareVersionStr, data=soft_id)
            processing_control.create_dataset(inputFilesStr, data=l1_input.input_files)
        numPixels = self.npix
        numScans = self.nscan
        numBands = self.nbands
        numBandsIR = self.nbandsir
        numLvlProf = self.nlvl
        n_refl_loc = self.n_refl_loc
        n_cloud_phase = self.n_cloud_phase
        spix = 0
        cpix = numPixels / 2
        epix = numPixels - 1
        cscan = numScans / 2

        tot_prod = len(self.l2_prod_names)
        self.tot_prod = tot_prod
        dm = np.zeros(3, dtype=np.int32)
        dm[0] = numScans
        dm[1] = numPixels
