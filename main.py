import argparse
import os
import time

from convl12 import convl12
from input_struc import Instr
from l2_struc import L2str
from libl1.filehandle import Filehandle, WRITE
from libl1.l1 import FATAL_ERROR
from msl12_input import msl12_input
from libl1.l1_input import L1_input_t


def main():

    l2rec = L2str()
    l1file = Filehandle()
    tgfile = Filehandle()
    aefile = Filehandle()
    ofile = Filehandle()

    parser = argparse.ArgumentParser()
    parser.add_argument("ifile")
    parser.add_argument("--sline", type=int)
    parser.add_argument("--eline", type=int)
    parser.add_argument("--spixl", type=int)
    parser.add_argument("--epixl", type=int)
    parser.add_argument("--l2prod")
    args = parser.parse_args()

    # Parse input parameters
    input = Instr("")
    l1_input = L1_input_t()
    msl12_input(args, l1file, input, l1_input)
    if not os.access(input.ifile, os.F_OK) or not os.access(input.ifile, os.R_OK):
        print(f"-E- l2gen: Input file {input.ifile} does not exist or cannot open.")
        exit(FATAL_ERROR)
    l1file.readl1()

    npix = l1file.npix
    # Set the end pixel if it was not set by command argument
    if l1_input.epixl == -1 or l1_input.epixl > l1file.npix:
        l1_input.epixl = l1file.npix
    if l1_input.eline == -1 or l1_input.eline > l1file.nscan:
        l1_input.eline = l1file.nscan
    if l1_input.spixl < 1:
        l1_input.spixl = 1
    if l1_input.sline < 1:
        l1_input.sline = 1
    spix = max(l1_input.spixl - 1, 0)
    epix = min(l1_input.epixl - 1, l1file.npix - 1)
    dpix = max(l1_input.dpixl, 1)
    sscan = max(l1_input.sline - 1, 0)
    escan = min(l1_input.eline - 1, l1file.nscan - 1)
    dscan = max(l1_input.dline, 1)

    if sscan > escan or spix > epix:
        print(f"-E- l2gen: scan and pixel limits make no sense.")
        print(f" start scan  = {sscan + 1}")
        print(f" end   scan  = {escan + 1}")
        print(f" start pixel = {spix + 1}")
        print(f" end   pixel = {epix + 1}")
        raise ValueError("FATAL_ERROR: Invalid scan and pixel limits.")

    l1file.spix = spix
    l1file.epix = epix

    # Open output file
    ofile.name = input.ofile
    ofile.mode = WRITE
    ofile.sensorID = l1file.sensorID
    ofile.nbands = l1file.nbands
    ofile.nbandsir = l1file.nbandsir
    ofile.nlvl = l1file.nlvl
    ofile.bindx = l1file.bindx
    ofile.ndets = l1file.ndets
    ofile.subsensorID = l1file.subsensorID
    ofile.spatialResolution = l1file.spatialResolution
    ofile.spix = spix
    ofile.epix = epix
    ofile.npix = (epix - spix) / dpix + 1
    ofile.length = l2rec.length
    ofile.nscan = (escan - sscan) / dscan + 1
    ofile.l2prod, input.l2prod
    ofile.def_l2prod, input.def_l2prod

    ofile.orbit_node_lon = l1file.orbit_node_lon
    ofile.orbit_number = l1file.orbit_number
    ofile.node_crossing_time = l1file.node_crossing_time
    ofile.private_data = l1file.private_data

    print(
        "Begin %s Version %d.%d.%d-%s Processing\n",
        PROGRAM,
        VERSION_MAJOR,
        VERSION_MINOR,
        VERSION_PATCH,
        GITSHA,
    )
    print(f"Sensor is {l1file.sensorInfo.SensorName}")
    print(f"Sensor ID is {l1file.sensorID}")
    print(f"Sensor has {l1file.nbands} reflective bands")
    print(f"Sensor has {l1file.nbandsir} emissive bands")
    print(f"Number of along-track detectors per band is {l1file.ndets}")
    print(f"Number of input pixels per scan is {l1file.npix}")
    print(f"Processing pixels {spix + 1} to {epix + 1} by {dpix}")
    print(f"Processing scans {sscan + 1} to {escan + 1} by {dscan}")

    start_time = time.time()
    print(f"Begin MSl12 processing at {time.ctime(start_time)}")

    # Convert the L1B radiances to L2
    # convl12()


if __name__ == "__main__":
    main()
