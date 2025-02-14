import json
from pathlib import Path
from collections import ChainMap

from libl1.l1 import Level1
from oel_hdf4.filetype.filetype import FileFormat
from oel_util.libgenutils.genutils_globals import want_verbose
from sensor_info import SensorInfo


def msl12_input(level1: Level1, cli_args: dict):
    l2genProgName = "msl12"
    l1_defaults_prefix = "msl12"
    OCDATAROOT = Path("./share")

    defaults = OCDATAROOT.joinpath("common", f"{l2genProgName}_defaults.json")

    if want_verbose:
        print(f"Loading default parameters from {defaults}")
    with open(defaults, "r") as f:
        defaults_settings = json.load(f)
    ifile = Path(cli_args["ifile"])

    format = FileFormat()
    format.get(ifile)

    level1.name = ifile
    level1.format = format.type

    sensor_info = SensorInfo(sensorId=format.sensor_id)
    sensor_info.rdsensorinfo(OCDATAROOT)

    level1.sensorinfo = sensor_info
    level1.wave = sensor_info.Lambda

    sensor_defaults = OCDATAROOT.joinpath(
        sensor_info.sensorDir, f"{l1_defaults_prefix}_defaults.json"
    )
    if want_verbose:
        print(
            f"Loading default parameters for {sensor_info.sensorName} from {sensor_defaults}"
        )
    with open(sensor_defaults) as f:
        sensor_defaults_settings = json.load(f)
    sensor_defaults_settings["naermodels"] = 80
    sensor_defaults_settings["aermodels"] = [
        "r30f95v01",
        "r30f80v01",
        "r30f50v01",
        "r30f30v01",
        "r30f20v01",
        "r30f10v01",
        "r30f05v01",
        "r30f02v01",
        "r30f01v01",
        "r30f00v01",
        "r50f95v01",
        "r50f80v01",
        "r50f50v01",
        "r50f30v01",
        "r50f20v01",
        "r50f10v01",
        "r50f05v01",
        "r50f02v01",
        "r50f01v01",
        "r50f00v01",
        "r70f95v01",
        "r70f80v01",
        "r70f50v01",
        "r70f30v01",
        "r70f20v01",
        "r70f10v01",
        "r70f05v01",
        "r70f02v01",
        "r70f01v01",
        "r70f00v01",
        "r75f95v01",
        "r75f80v01",
        "r75f50v01",
        "r75f30v01",
        "r75f20v01",
        "r75f10v01",
        "r75f05v01",
        "r75f02v01",
        "r75f01v01",
        "r75f00v01",
        "r80f95v01",
        "r80f80v01",
        "r80f50v01",
        "r80f30v01",
        "r80f20v01",
        "r80f10v01",
        "r80f05v01",
        "r80f02v01",
        "r80f01v01",
        "r80f00v01",
        "r85f95v01",
        "r85f80v01",
        "r85f50v01",
        "r85f30v01",
        "r85f20v01",
        "r85f10v01",
        "r85f05v01",
        "r85f02v01",
        "r85f01v01",
        "r85f00v01",
        "r90f95v01",
        "r90f80v01",
        "r90f50v01",
        "r90f30v01",
        "r90f20v01",
        "r90f10v01",
        "r90f05v01",
        "r90f02v01",
        "r90f01v01",
        "r90f00v01",
        "r95f95v01",
        "r95f80v01",
        "r95f50v01",
        "r95f30v01",
        "r95f20v01",
        "r95f10v01",
        "r95f05v01",
        "r95f02v01",
        "r95f01v01",
        "r95f00v01",
    ]
    # load sensor suite file
    localSuite = "OC"
    sensor_defaults_oc = OCDATAROOT.joinpath(
        sensor_info.sensorDir, f"{l2genProgName}_defaults_{localSuite}.json"
    )
    if want_verbose:
        print(f"Loading parameters for suite {localSuite} from {sensor_defaults_oc}")
    with open(sensor_defaults_oc) as f:
        sensor_defaults_oc_settings = json.load(f)

    if want_verbose & defaults_settings["deflate"] > 0:
        print(
            f"Internal data compression requested at compression level: {defaults_settings["deflate"]}"
        )
    environ = {"OCDATAROOT": OCDATAROOT}
    settings = ChainMap(
        cli_args,
        environ,
        sensor_defaults_oc_settings,
        sensor_defaults_settings,
        defaults_settings,
    )
    return settings
