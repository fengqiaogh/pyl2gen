from enum import Enum
from pathlib import Path

import h5py

from genutils_globals import want_verbose
from sensorDefs import GOCI


class FileType(Enum):
    FT_INVALID = -1
    FT_UNKNOWN = 0
    FT_AVIRIS = 1
    FT_CLASSAVHRR = 2
    FT_CZCSL1A = 3
    FT_GOCIL1B = 4
    FT_HICOL1B = 5
    FT_L1BNCDF = 6
    FT_L1HDF = 7
    FT_L1XCAL = 8
    FT_L2HDF = 9
    FT_L2NCDF = 10
    FT_L3BIN = 11
    FT_L3MAP = 12
    FT_MERISCC = 13
    FT_MERISL1B = 14
    FT_MERISL2 = 15
    FT_MERISL1BSAFE = 16
    FT_MODISGEO = 17
    FT_MODISL1B = 18
    FT_MOSL1B = 19
    FT_OCM2L1B = 20
    FT_OCML1B = 21
    FT_OCML1BDB = 22
    FT_OCTSL1A = 23
    FT_OCTSL1B = 24
    FT_OLCI = 25
    FT_OLCIGEO = 26
    FT_OLIL1B = 27
    FT_OCIA = 28
    FT_OCIL1B = 29
    FT_OCIS = 30
    FT_OSMIL1A = 31
    FT_PRISM = 32
    FT_SEAWIFSL1A = 33
    FT_VIIRSGEO = 34
    FT_VIIRSGEONC = 35
    FT_VIIRSL1A = 36
    FT_VIIRSL1B = 37
    FT_VIIRSL1BNC = 38
    FT_SGLI = 39
    FT_L5TML1B = 40
    FT_L7ETML1B = 41
    FT_MSIL1C = 42
    FT_HAWKEYEL1A = 43
    FT_MISR = 44
    FT_SEABASSRRS = 45
    FT_SPEXONE = 46
    FT_HARP2 = 47
    FT_HARP = 48
    FT_HKT = 49
    FT_L1C = 50


class FileFormat:
    def __init__(self):
        self.type = FileType.FT_INVALID
        self.sensor_id = -1
        self.subsensor_id = -1

    def get(self, filename: Path):
        if filename.exists():
            pass
        if filename.suffix == ".he5":
            # check for GOCI
            with h5py.File(filename) as f:
                group = f["HDFEOS/POINTS/Scene Header"]
                if "GOCI Level-1B Data" in str(group.attrs["Scene Title"]):
                    self.type = FileType.FT_GOCIL1B
                    self.sensor_id = GOCI
                    self.subsensor_id = -1
                    if want_verbose:
                        print(f"Input file {filename} is a GOCI L1B HDF 5 file.")
