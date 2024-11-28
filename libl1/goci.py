import math
from datetime import datetime
import numpy as np
import h5py
from libl1.l1 import L1str
from oel_util.libgenutils.genutils_globals import want_verbose
from pyproj import CRS, Transformer

# from ..get_zenaz import get_zenaz
# from ..sunangs import sunangs


class GOCIL1:
    nbands = 8
    NAV_GRP = "HDFEOS/POINTS/Navigation for GOCI/Data"
    TABLE_NAME = "Navigation for GOCI"

    def __init__(self):
        self.slot_asg = None
        self.slot_rel_time = None
        self.sat_pos = np.zeros(3)

    def open(self, file):
        src_path = file.name
        fields = [
            "Band 1 Image Pixel Values",
            "Band 2 Image Pixel Values",
            "Band 3 Image Pixel Values",
            "Band 4 Image Pixel Values",
            "Band 5 Image Pixel Values",
            "Band 6 Image Pixel Values",
            "Band 7 Image Pixel Values",
            "Band 8 Image Pixel Values",
        ]
        dims = np.zeros(3, dtype=np.int32)
        with h5py.File(src_path) as f:

            dims[1] = f["HDFEOS/POINTS/Scene Header"].attrs["number of columns"][0]
            dims[0] = f["HDFEOS/POINTS/Scene Header"].attrs["number of rows"][0]
            self.npixels = dims[1]
            self.nscans = dims[0]

            cpos = f["HDFEOS/POINTS/Ephemeris"].attrs[
                "Satellite position XYZ (ECEF) at scene center time"
            ]
            self.sat_pos[2] = cpos[0] * np.sin(cpos[2]) / 1000
            radius_in_xy = cpos[0] * np.cos(cpos[2]) / 1000
            self.sat_pos[1] = radius_in_xy * np.sin(cpos[1])
            self.sat_pos[0] = radius_in_xy * np.cos(cpos[1])
            time_str = f["HDFEOS/POINTS/Ephemeris"].attrs["Scene Start time"]

        self.slot_asg, self.slot_rel_time = slot_init(src_path, dims)
        if want_verbose:
            print(f"GOCI Level-1B {file.name}")
        file.npix = self.npixels
        file.nscan = self.nscans
        file.bands = self.nbands
        file.spatialResolution = "500 m"

        self.get_datetime(time_str)
        if want_verbose:
            print(
                f"GOCI Scene Start time: {self.year:4d}-{self.month:02d}-{self.day:02d} {self.doy:03d} {self.hour:02d}:{self.minute:02d}:{self.second:02d} {self.base_msec}"
            )
            print(
                f"GOCI file has {file.nbands} bands, {file.npix} samples, {file.nscan} lines"
            )

        with h5py.File(src_path) as f:
            self.proj4_open(f)
        if want_verbose:
            print("GOCI using internal navigation")

    def proj4_open(self, f):
        # read Map Projection attributes
        map_proj = "/HDFEOS/POINTS/Map Projection"
        centralLat = f[map_proj].attrs["Central Latitude (parallel)"][0]
        centralLon = f[map_proj].attrs["Central Longitude (meridian)"][0]
        equitorialRadius = f[map_proj].attrs["Equitorial radius of Earth ellipsoid"][0]
        polarRadius = f[map_proj].attrs["Polar radius of Earth ellipsoid"][0]

        # read Scene Header attributes
        sce_hea = "/HDFEOS/POINTS/Scene Header"
        llLat = f[sce_hea].attrs["Scene lower-left latitude"]
        llLon = f[sce_hea].attrs["Scene lower-left longitude"]
        urLat = f[sce_hea].attrs["Scene upper-right latitude"]
        urLon = f[sce_hea].attrs["Scene upper-right longitude"]
        projStr = f"+proj=ortho +ellps=WGS84 +datum=WGS84 +lon_0={centralLon} +lat_0={centralLat} +a={equitorialRadius} +b={polarRadius}"
        target_proj = CRS.from_proj4(projStr)
        source_proj = CRS.from_proj4("+proj=longlat +ellps=WGS84 +datum=WGS84")
        transformer = Transformer.from_crs(source_proj, target_proj, always_xy=True)

        # Calculate start and delta for the grid
        ll_trans = transformer.transform(llLon, llLat)
        ur_trans = transformer.transform(urLon, urLat)

        self.startX = ll_trans[0]
        self.startY = ur_trans[1]
        self.deltaX = (ur_trans[0] - ll_trans[0]) / self.npixels
        self.deltaY = (ll_trans[1] - ur_trans[1]) / self.nscans

    def get_datetime(self, time_str):
        time_str = time_str.decode("utf-8")
        time_obj = datetime.strptime(time_str, "%d-%b-%Y %H:%M:%S.%f")
        self.year = time_obj.year
        self.month = time_obj.month
        self.day = time_obj.day
        self.hour = time_obj.hour
        self.minute = time_obj.minute
        self.second = time_obj.second
        self.doy = time_obj.timetuple().tm_yday
        self.base_msec = 1000 * (self.second + 60 * (self.minute + 60 * self.hour))


def slot_init(file_path, dims):
    nbnd, nslot = 8, 16
    trg_bnd = 7
    with h5py.File(file_path) as f:
        group = f["HDFEOS/POINTS/Navigation for GOCI/Data"]

        # 检查表的字段和记录
        table = group["Navigation for GOCI"]
        nfields = len(table.dtype.names)
        nrecords = len(table)

        print(f"# fields: {nfields}, # records: {nrecords}")
        # 读取表格数据
        navigation_table = table[()]

    slot_nav_data = []
    for row in navigation_table:
        slot = SlotNav()
        slot.band_num = row["Band number"]
        slot.slot_num = row["Slot number"]
        slot.rel_time = row["Relative time"]
        slot.sc_att = row["Spacecraft attitude"]
        slot.xo = row["XO"]
        slot.yo = row["YO"]
        slot.xs = row["XS"]
        slot.ys = row["YS"]
        slot.xpo = row["XPO"]
        slot.ypo = row["YPO"]
        slot.xps = row["XPS"]
        slot.yps = row["YPS"]
        slot.num_a_parm = row["Number of valid A parameters"]
        slot.a_parm = row["A parameters value"]
        slot.num_b_parm = row["Number of valid B parameters"]
        slot.b_parm = row["B parameters value"]
        slot.num_c_parm = row["Number of valid C parameters"]
        slot.c_parm = row["C parameters value"]
        slot.num_d_parm = row["Number of valid D parameters"]
        slot.d_parm = row["D parameters value"]
        slot.num_ap_parm = row["Number of valid A prime parameters"]
        slot.ap_parm = row["A prime parameters value"]
        slot.num_bp_parm = row["Number of valid B prime parameters"]
        slot.bp_parm = row["B prime parameters value"]
        slot.num_cp_parm = row["Number of valid C prime parameters"]
        slot.cp_parm = row["C prime parameters value"]
        slot.num_dp_parm = row["Number of valid D prime parameters"]
        slot.dp_parm = row["D prime parameters value"]
        slot_nav_data.append(slot)

    print("Begin GOCI slot assignment")

    nlin = dims[0]
    npix = dims[1]
    slot_asg = np.zeros((nlin, npix), dtype=np.int32)

    step = 18
    nsx = int(2 + npix // step)
    nsy = int(2 + nlin // step)

    slot_asg_sml = np.zeros((nsy, nsx), dtype=np.uint8)
    bnd_tile_lut = np.full((nbnd, nslot), 254, dtype=np.uint8)

    # Supergrid slot assignment
    for iy in range(nsy):
        ilin = iy * step
        for ix in range(nsx):
            ipix = ix * step
            minrad = 200
            for itile in range(nslot):
                crad = goci_slot_nav(
                    ipix, ilin, 7, itile, slot_nav_data, nbnd, nslot, bnd_tile_lut
                )
                if crad < minrad:
                    slot_asg_sml[iy, ix] = itile
                    minrad = crad

    # Full grid slot assignment
    for iy in range(nsy - 1):
        lin_st = iy * step
        lin_en = min((iy + 1) * step, nlin)
        for ix in range(nsx - 1):
            pix_st = ix * step
            pix_en = min((ix + 1) * step, npix)

            box_pts = [
                slot_asg_sml[iy, ix],
                slot_asg_sml[iy, ix + 1],
                slot_asg_sml[iy + 1, ix],
                slot_asg_sml[iy + 1, ix + 1],
            ]
            box_pts.sort()

            if box_pts[0] == box_pts[3]:
                slot_asg[lin_st:lin_en, pix_st:pix_en] = box_pts[0]
            else:
                for ilin in range(lin_st, lin_en):
                    for ipix in range(pix_st, pix_en):
                        minrad = 200
                        curtil = -1
                        for itile in box_pts:
                            if itile != curtil:
                                curtil = itile
                                crad = goci_slot_nav(
                                    ipix,
                                    ilin,
                                    7,
                                    curtil,
                                    slot_nav_data,
                                    nbnd,
                                    nslot,
                                    bnd_tile_lut,
                                )
                                if crad < minrad:
                                    slot_asg[ilin, ipix] = curtil
                                    minrad = crad

    # Compute time offsets per slot
    slot_rel_time = np.zeros(nslot, dtype=np.float32)
    for itile in range(nslot):
        min_t, max_t = 5000, -5000
        for ibnd in range(nbnd):
            ilut = bnd_tile_lut[ibnd, itile]
            rel_t = slot_nav_data[ilut].rel_time
            min_t = min(min_t, rel_t)
            max_t = max(max_t, rel_t)
        slot_rel_time[itile] = (min_t + max_t) / 2.0

    print("GOCI slot, time assignments completed.")
    return slot_asg, slot_rel_time

    def read(self, file, l1rec: L1str):
        npix = file.npix
        nbands = file.nbands
        msec = base_msec
        # 将字符串转换为 datetime 对象
        date_str = "20200404051643"
        date_obj = datetime.strptime(date_str, "%Y%m%d%H%M%S")
        # 计算当天零点的 datetime 对象
        midnight = datetime(date_obj.year, date_obj.month, date_obj.day)

        # 计算给定时间距离当天零点的秒数
        seconds_since_midnight = (date_obj - midnight).total_seconds()
        gmt = seconds_since_midnight / 3600.0

        solz, sola = sunangs(self.year, self.doy, gmt, self.lon, self.lat)
        sola[sola > 180] = sola[sola > 180] - 360
        senz, sena = get_zenaz(self.sat_pos, self.lon, self.lat)
        sena[sena > 180] = sena[sena > 180] - 360


def goci_slot_nav(ipix, ilin, bnd, itile, slot_nav, nbnd, nslot, bnd_tile_lut):
    """
    Transform a GOCI scene point into the point on a specific tile and return
    the largest linear normalized distance from the center.

    Parameters:
        ipix (int): Pixel to transform.
        ilin (int): Line to transform.
        bnd (int): Band number.
        itile (int): Tile or slot number of GOCI.
        slot_nav (list of dict): List of structures with transform coefficients and normalization values.
        nbnd (int): Number of bands.
        nslot (int): Number of slots.
        bnd_tile_lut (numpy array): Storage for a look-up for proper element in slot_nav.
        lindist (numpy array): The largest normalized linear distance from tile center (output).

    Returns:
        int: 0 if all is OK.
    """
    nlut = nbnd * nslot

    if bnd_tile_lut[0, 0] == 254:
        for ilut in range(nlut):
            band_num = slot_nav[ilut].band_num
            slot_num = slot_nav[ilut].slot_num
            bnd_tile_lut[band_num, slot_num] = ilut

    # get normalization coefficients for scene to normalized scene
    ilut = bnd_tile_lut[bnd, itile]
    xo = slot_nav[ilut].xo
    yo = slot_nav[ilut].yo
    xs = slot_nav[ilut].xs
    ys = slot_nav[ilut].ys

    # Extract parameters for the transform
    a_parm = slot_nav[ilut].a_parm
    b_parm = slot_nav[ilut].b_parm
    c_parm = slot_nav[ilut].c_parm
    d_parm = slot_nav[ilut].d_parm

    num_a_parm = slot_nav[ilut].num_a_parm
    num_b_parm = slot_nav[ilut].num_b_parm
    num_c_parm = slot_nav[ilut].num_c_parm
    num_d_parm = slot_nav[ilut].num_d_parm

    # Normalize scene location
    xn = (ipix - xo) / xs
    yn = (ilin - yo) / ys

    # Create the vector for polynomial terms
    vec = np.zeros(16)
    vec[0] = 1.0
    vec[1] = xn
    vec[2] = yn
    vec[3] = xn * yn
    vec[4] = xn**2
    vec[5] = yn**2
    vec[6] = (xn**2) * yn
    vec[7] = (yn**2) * xn
    vec[8] = (xn**2) * (yn**2)

    # X transform
    numer = np.dot(vec[:num_a_parm], a_parm[:num_a_parm])
    denom = 1.0 + np.dot(vec[1 : num_b_parm + 1], b_parm[:num_b_parm])
    xpn = numer / denom

    # Y transform
    numer = np.dot(vec[:num_c_parm], c_parm[:num_c_parm])
    denom = 1.0 + np.dot(vec[1 : num_d_parm + 1], d_parm[:num_d_parm])
    ypn = numer / denom

    lindist = 200.0

    # Find the largest linear normalized distance from the center
    dist = [1.0 - xpn, 1.0 + xpn, 1.0 - ypn, 1.0 + ypn]
    lindist = 1.0 - min(dist)

    return lindist


class SlotNav:
    def __init__(self):
        self.band_num = None
        self.slot_num = None
        self.rel_time = None
        self.sc_att = None
        self.xo = None
        self.yo = None
        self.xs = None
        self.ys = None
        self.xpo = None
        self.ypo = None
        self.xps = None
        self.yps = None
        self.num_a_parm = None
        self.a_parm = None
        self.num_b_parm = None
        self.b_parm = None
        self.num_c_parm = None
        self.c_parm = None
        self.num_d_parm = None
        self.d_parm = None
        self.num_ap_parm = None
        self.ap_parm = None
        self.num_bp_parm = None
        self.bp_parm = None
        self.num_cp_parm = None
        self.cp_parm = None
        self.num_dp_parm = None
        self.dp_parm = None
