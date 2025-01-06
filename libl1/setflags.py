import numpy as np
from libl1.l1 import ON


class Flag:
    def __init__(self, shape):
        # view angle limits
        self.senzmax = np.zeros(shape, dtype=np.int16)
        self.solzmax = np.zeros(shape, dtype=np.int16)

        # clouds
        self.cloud = np.zeros(shape, dtype=np.int16)

    def set(self, wave, settings):
        ib412 = np.abs(412 - wave).argmin()
        ib555 = np.abs(555 - wave).argmin()
        ib670 = np.abs(670 - wave).argmin()
        ib865 = np.abs(865 - wave).argmin()
        ibcloud = np.abs(settings["cloud_wave"] - wave).argmin()
        print(f"Using {wave[ibcloud]} nm channel for cloud flagging over water.")
        print(f"Using {wave[ib412]} nm channel for cloud flagging over land.")

        self.senzmax[senz > settings["satzen"]] = ON
        self.solzmax[solz > settings["sunzen"]] = ON

        # Check for clouds (daytime only)
        cloud_albedo = (
            rhos[:, :, ibcloud] - TLg[:, :, ibcloud] * np.pi / mu0 / Fo[ibcloud]
        )
        albedo = settings["albedo"]
        self.cloud[cloud_albedo > albedo] = ON
