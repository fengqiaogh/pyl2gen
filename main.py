import argparse
import os

import sys
import time
from ancillary import Ancillary
from convl12 import convl12
from l2_struc import Level2
from l12_parms import PROGRAM
from libl1.l1 import Level1
from msl12_input import msl12_input
from oel_util.libgenutils.genutils_globals import want_verbose
from version import GITSHA, VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH


def main():
    parser = argparse.ArgumentParser(
        description="multi-sensor level-1b to level-2 conversion"
    )
    parser.add_argument("ifile", type=str, help="input L1 file name")
    parser.add_argument("--ofile", type=str, default="output", help="output file name")
    parser.add_argument("--sline", default=1, type=int, help="start line number")
    parser.add_argument(
        "--eline", default=-1, type=int, help="end line number (-1=the last line)"
    )
    parser.add_argument("--spixl", default=1, type=int, help="start pixel number")
    parser.add_argument(
        "--epixl", default=-1, type=int, help="end pixel number (-1=the last pixel)"
    )
    parser.add_argument(
        "--l2prod", type=str, help="L2 products to be included in ofile"
    )
    parser.add_argument(
        "--cloud_wave",
        default=865,
        type=int,
        help="wavelength of cloud reflectance test",
    )
    parser.add_argument(
        "--extreme_glint",
        default=0.03,
        type=float,
        help="extreme sun glint threshold",
    )
    if want_verbose:
        print("Loading command line parameters")

    args = parser.parse_args()

    # 创建参数字典
    cli_args = {key: value for key, value in vars(args).items() if value}

    level1 = Level1()
    settings = msl12_input(level1, cli_args)
    if not os.access(settings["ifile"], os.F_OK) or not os.access(
        settings["ifile"], os.R_OK
    ):
        print(f"l2gen: Input file {settings["ifile"]} does not exist or cannot open.")
        sys.exit()

    level1.read(
        settings["sline"], settings["eline"], settings["spixl"], settings["epixl"]
    )

    print(
        f"Begin {PROGRAM} Version {VERSION_MAJOR}.{VERSION_MINOR}.{VERSION_PATCH}-{GITSHA} Processing",
    )
    print(f"Sensor is {level1.sensorinfo.sensorName}")
    print(f"Sensor ID is {level1.sensorinfo.sensorId}")
    print(f"Sensor has {level1.sensorinfo.Nbands} reflective bands")
    # print(f"Sensor has {level1.nbandsir} emissive bands")
    # print(f"Number of along-track detectors per band is {l1file.ndets}")
    print(f"Number of input pixels per scan is {level1.npix}")
    print(f"Processing pixels {settings["spixl"]} to {settings["epixl"]}")
    print(f"Processing scans {settings["sline"]} to {settings["eline"]}")

    start_time = time.time()
    print(f"Begin MSl12 processing at {time.ctime(start_time)}")

    # Load ancillary data
    anc = Ancillary()
    anc.set(level1.lat, level1.lon, level1.datetime, settings)

    level1.load(settings, anc)

    level2 = Level2()

    # Convert the L1B radiances to L2
    convl12(level1, level2, anc, settings)


if __name__ == "__main__":
    main()
