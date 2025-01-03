import json
from operator import le
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
    level1.sensor_id = format.sensor_id
    level1.subsensor_id = format.subsensor_id

    sensor_info = SensorInfo(sensorId=level1.sensor_id)
    sensor_info.rdsensorinfo(OCDATAROOT)

    # level1

    sensor_defaults = OCDATAROOT.joinpath(
        sensor_info.sensorDir, f"{l1_defaults_prefix}_defaults.json"
    )
    if want_verbose:
        print(
            f"Loading default parameters for {sensor_info.sensorName} from {sensor_defaults}"
        )
    with open(sensor_defaults) as f:
        sensor_defaults_settings = json.load(f)

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
