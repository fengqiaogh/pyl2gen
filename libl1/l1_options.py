import json
from pathlib import Path

from oel_hdf4.filetype.filetype import FileFormat
from oel_util.libgenutils.genutils_globals import want_verbose
from libl1.filehandle import Filehandle
from sensorInfo import SensorInfo


def l1_read_default_files(l1file: Filehandle, ifile: Path):
    l1_defaults_prefix = "msl12"
    l1file.name = ifile
    format = FileFormat()
    format.get(l1file.name)
    l1file.format = format.type
    l1file.sensorID = format.sensor_id
    l1file.subsensorID = format.subsensor_id
    l1file.sensorInfo = SensorInfo(sensorId=l1file.sensorID)
    l1file.nbands = l1file.sensorInfo.Nbands
    tmpStr = l1file.sensorInfo.OCDATAROOT.joinpath(
        l1file.sensorInfo.sensorDir, f"{l1_defaults_prefix}_defaults.json"
    )
    if want_verbose:
        print(
            f"Loading default parameters for {l1file.sensorInfo.sensorName} from {tmpStr}"
        )
    with open(tmpStr, "r") as f:
        d = json.load(f)
