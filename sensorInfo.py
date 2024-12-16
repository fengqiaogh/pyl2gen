import json
import os
import sys
from pathlib import Path

import numpy as np

from esdist import esdist
from oel_util.libgenutils.genutils_globals import want_verbose

SENSOR_NUM = 37
SENSOR_NAME = [
    "SeaWiFS",
    "MOS",
    "OCTS",
    "AVHRR",
    "OSMI",
    "CZCS",
    "MODIST",
    "MODISA",
    "OCM1",
    "OCM2",
    "MERIS",
    "VIIRSN",
    "OCRVC",
    "HICO",
    "GOCI",
    "OLIL8",
    "Aquarius",
    "OCIA",
    "AVIRIS",
    "PRISM",
    "OLCIS3A",
    "SGLI",
    "MSIS2A",
    "L5TM",
    "L7ETMP",
    "VIIRSJ1",
    "MSIS2B",
    "HAWKEYE",
    "MISR",
    "OLCIS3B",
    "OCI",
    "OCIS",
    "VIIRSJ2",
    "OLIL9",
    "SPEXONE",
    "HARP2",
    "HARP",
]
SENSOR_DIR = [
    "seawifs",
    "mos",
    "octs",
    "avhrr",
    "osmi",
    "czcs",
    "modis",
    "modis",
    "ocm1",
    "ocm2",
    "meris",
    "viirs",
    "ocrvc",
    "hico",
    "goci",
    "oli",
    "aquarius",
    "ocia",
    "aviris",
    "prism",
    "olci",
    "sgli",
    "msi",
    "l5tm",
    "l7etmp",
    "viirs",
    "msi",
    "hawkeye",
    "misr",
    "olci",
    "oci",
    "ocis",
    "viirs",
    "oli",
    "spexone",
    "harp2",
    "harp",
]
SUBSENSOR_DIR = [
    "gac",
    "lac",
    "terra",
    "aqua",
    "npp",
    "j1",
    "s2a",
    "s2b",
    "s3a",
    "s3b",
    "j2",
    "l8",
    "l9",
]


class SensorInfo:
    def __init__(self, sensorId, subsensorId=-1):
        self.sensorId = sensorId
        self.subsensorId = subsensorId
        self.sensorId2SensorName()
        self.sensorId2SensorDir()
        self.subsensorId2SubsensorDir()
        self.rdsensorinfo()

    def sensorId2SensorName(self):
        if self.sensorId < 0:
            self.sensorName = None
        elif self.sensorId >= SENSOR_NUM:
            self.sensorName = None
        else:
            self.sensorName = SENSOR_NAME[self.sensorId]

    def sensorId2SensorDir(self):
        if self.sensorId < 0:
            self.sensorDir = None
        elif self.sensorId >= SENSOR_NUM:
            self.sensorDir = None
        else:
            self.sensorDir = SENSOR_DIR[self.sensorId]

    def subsensorId2SubsensorDir(self):
        if self.subsensorId < 0:
            self.subsensorDir = None
        elif self.subsensorId >= SENSOR_NUM:
            self.subsensorDir = None
        else:
            self.subsensorDir = SUBSENSOR_DIR[self.subsensorId]

    def rdsensorinfo(self):
        self.OCDATAROOT = Path.cwd().joinpath("share")
        # ifile = Path(args.ifile)
        # tmpStr = dataRoot.joinpath("common", f"{l2genProgName}_defaults.par")

        # self.OCDATAROOT = os.environ.get("OCDATAROOT", None)
        if not self.OCDATAROOT.exists():
            print("OCDATAROOT environment variable is not defined.")
            sys.exit(1)

        if want_verbose:
            print(f"Loading characteristics for {self.sensorName}")

        # 打开 JSON 文件
        with open(
            self.OCDATAROOT.joinpath(self.sensorName.lower(), "msl12_sensor_info.json")
        ) as file:
            # 读取文件内容并将其转换为字典
            data = json.load(file)

            for key in data.keys():
                setattr(self, key, np.array(data[key]))
