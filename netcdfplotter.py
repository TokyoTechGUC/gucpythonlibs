import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from glob import glob
from wrf import getvar
import numpy as np
from datetime import datetime, timedelta
import os.path
from os import mkdir, getcwd
import subprocess
import ffmpeg


class ncdfproc:
    """
    Author: Alvin C.G. Varquez
    Updated: November 4, 2021
    Note: Currently supports only 3D variables of WRF and WRF processed with CDO.

    Usage example:

    test = ncdfproc(files='*',outdir='./',nctype=None,overwrite=False)

    Arguments:
        files: A required string (wildcards acceptable) of the
            netCDF files to be processed.
        outdir (optional): Path of the output directory. Default
            is the current working directory.
        nctype (optional): String corresponding to the type of the
            netCDF file. Setting this to None will make the program
            guess the type. Set to either 'wrf' (standard WRF
            outputs), 'cdowrf' (WRF outputs processed using CDO;
            requires XTIME,XLONG, and XLAT variables).
            Defaults to None.
        overwrite (optional): Option to overwrite the images
            to be created. Default is False.

    Functions:
        Create 2D png images of 3D variables with
        dimensions (time,latitude, longitude) with draw_map().

        test.draw_map(variable='T2',figsize=(8,6),
            coastwidth=0.0,dpi=100,title='',adjust_hr=0,**kwargs)

        Optional arguments:
            variable: Defaults to 'T2' for WRF.
            figsize: The size of the figure (plt.subplots option).
                Defaults to (8,6).
            coastwidth: The thickness of the coastlines to
                be drawn. Defaults to 0.0.
            title: The prefix added to the title of the image.
                Defaults to ''.
            adjust_hr: Shifts the time coordinate in the
                output image. For example, localtime=GMT+adjust_hr.
                Defaults to 0 or GMT.
            **kwargs: Optional keyword arguments of the pcolormesh
                of matplotlib. Example, you may adjust the range by
                adding the vmin and vmax.

        To create a timelapse video in mp4 of a sequence of images.

        test.create_video(string,framerate=4,out='out.mp4',overwrite_output=True)

        Arguments:
            string (required): A wildcard of the images.
                For example: '*.png'. Note, an error will happen
                if the image extension is not addded after the
                wildcard (e.g. '*' will result to an error).
            framerate (optional): The framerate of the video.
                Higher framerate will lead to faster video.
                Defaults to 4 (or 4 images per second).
            out (optional): The output mp4 filename.
                Defaults to 'out.mp4'.
            overwrite_output: Whether to overwrite an existing
                mp4 file entered in the out argument.
                Defaults to True.
    """

    def __init__(self, files, outdir="./", nctype=None, overwrite=False):
        def check_type(self):
            ftypes = {"HFX": "wrf", "UTCI": "wrf"}
            if self.nctype == None:
                print(
                    "No type specified. Determining type from either WRF, CESM (and sub-models), or CDO outputs."
                )
                df = nc.Dataset(self.files[0], "r")
                for key, value in ftypes.items():
                    try:
                        df.variables[key]
                        self.nctype = value
                        if len(df.variables["XLAT"].shape) == 2:
                            self.nctype = "cdowrf"
                        print(
                            f"{self.filecount} {value.upper()} file(s) detected (nctype={self.nctype})."
                        )
                        break
                    except:
                        pass
                if self.nctype == None:
                    raise KeyError("File type currently unsupported.")
            return

        self.files = glob(files)
        self.outdir = outdir
        if not os.path.exists(outdir):
            mkdir(self.outdir)
        self.filecount = len(self.files)
        self.nctype = nctype
        self.overwrite = overwrite
        check_type(self)

    def draw_map(
        self,
        variable=None,
        figsize=(8, 6),
        coastwidth=0.0,
        dpi=100,
        title="",
        adjust_hr=0,
        **kwargs,
    ):
        self.coastwidth = coastwidth
        self.variable = variable
        self.title = title
        self.fig, self.ax = plt.subplots(
            figsize=figsize, dpi=dpi, constrained_layout=False
        )
        self.ax.clear()

        def gen_time(time, source, startdate, adjust_hr):
            units = source.units
            if "minutes since simulation start" in source.description:
                time = datetime.strptime(startdate, "%Y-%m-%d_%H:%M:00") + timedelta(
                    minutes=time, hours=adjust_hr
                )
            else:
                time = datetime.strptime(
                    units, "minutes since %Y-%m-%d %H:%M:00"
                ) + timedelta(minutes=time, hours=adjust_hr)
            return time

        def construct_shell_cdowrf(self, **kwargs):
            df = nc.Dataset(self.files[0], "r")
            lats = df.variables["XLAT"][:, 0]
            lons = df.variables["XLONG"][0, :]
            ny = lats.shape[0]
            nx = lons.shape[0]
            m = Basemap(
                llcrnrlon=np.min(lons),
                llcrnrlat=np.min(lats),
                urcrnrlon=np.max(lons),
                urcrnrlat=np.max(lats),
                resolution="h",
                ax=self.ax,
            )
            lons1, lats1 = m.makegrid(nx, ny)
            self.x, self.y = m(lons1, lats1)
            return m

        def draw_cdowrf(self, adjust_hr, **kwargs):
            degree = r"$\degree$"
            if self.variable == None:
                print("Variable not specified. Defaulting to T2")
                self.variable = "T2"
            print(f"Creating maps for {self.variable}")
            m = construct_shell_cdowrf(self, **kwargs)
            for count, ifil in enumerate(self.files):
                df = nc.Dataset(ifil, "r")
                for it in range(df.dimensions["XTIME"].size):
                    fig = self.fig
                    ax = self.ax
                    time = df.variables["XTIME"]
                    time = gen_time(
                        int(time[it]), time, df.SIMULATION_START_DATE, adjust_hr
                    )
                    if self.title == "":
                        self.title = self.variable
                    outfile = f'{self.outdir}/{self.title.replace(" ","_")}_{time.strftime("%Y%m%d%H%M")}_UTCplus_{adjust_hr:02d}.png'
                    if os.path.exists(outfile) & (not self.overwrite):
                        print(
                            f"Skipping {outfile} since it exists. To overwrite, set overwrite=True"
                        )
                        continue
                    if self.variable == "tc2":
                        var = df.variables["T2"]
                        ax1 = ax.pcolormesh(
                            self.x,
                            self.y,
                            var[it, :, :] - 273.15,
                            shading="auto",
                            **kwargs,
                        )
                    else:
                        var = df.variables[self.variable]
                        ax1 = ax.pcolormesh(
                            self.x, self.y, var[it, :, :], shading="auto", **kwargs
                        )
                    m.drawcoastlines(linewidth=self.coastwidth)
                    ax.set_title(
                        f'{self.title} {time.strftime("%Y/%m/%d %H:%M")} (UTC+{adjust_hr})'
                    )
                    if (count == 0) & (it == 0):
                        cbar = fig.colorbar(ax1, ax=ax, orientation="horizontal")
                        cbar.ax.set_title(
                            f'{var.description.title()} ({degree+"C" if self.variable=="tc2" else var.units})'
                        )
                    else:
                        cbar.update_normal(ax1)
                    fig.savefig(outfile, bbox_inches="tight", transparent=False)
                    ax1 = None
                    del (fig, ax)

        def construct_shell_wrf(self, **kwargs):
            df = nc.Dataset(self.files[0], "r")
            lats = df.variables["XLAT"][0, :, 0]
            lons = df.variables["XLONG"][0, 0, :]
            ny = lats.shape[0]
            nx = lons.shape[0]
            m = Basemap(
                llcrnrlon=np.min(lons),
                llcrnrlat=np.min(lats),
                urcrnrlon=np.max(lons),
                urcrnrlat=np.max(lats),
                resolution="h",
                ax=self.ax,
            )
            lons1, lats1 = m.makegrid(nx, ny)
            self.x, self.y = m(lons1, lats1)
            return m

        def draw_wrf(self, adjust_hr, **kwargs):
            degree = r"$\degree$"
            if self.variable == None:
                print("Variable not specified. Defaulting to T2")
                self.variable = "T2"
            print(f"Creating maps for {self.variable}")
            m = construct_shell_wrf(self, **kwargs)
            for count, ifil in enumerate(self.files):
                df = nc.Dataset(ifil, "r")
                for it in range(df.dimensions["Time"].size):
                    fig = self.fig
                    ax = self.ax
                    time = getvar(df, "XTIME", timeidx=it)
                    time = gen_time(
                        int(time), time, df.SIMULATION_START_DATE, adjust_hr
                    )
                    if self.title == "":
                        self.title = self.variable
                    outfile = f'{self.outdir}/{self.title.replace(" ","_")}_{time.strftime("%Y%m%d%H%M")}_UTCplus_{adjust_hr:02d}.png'
                    if os.path.exists(outfile) & (not self.overwrite):
                        print(not self.overwrite)
                        print(
                            f"Skipping {outfile} since it exists. To overwrite, set overwrite=True"
                        )
                        continue
                    if self.variable == "tc2":
                        var = getvar(df, "T2", timeidx=it)
                        ax1 = ax.pcolormesh(
                            self.x, self.y, var - 273.15, shading="auto", **kwargs
                        )
                    else:
                        var = getvar(df, self.variable, timeidx=it)
                        ax1 = ax.pcolormesh(
                            self.x, self.y, var, shading="auto", **kwargs
                        )
                    m.drawcoastlines(linewidth=self.coastwidth)
                    ax.set_title(
                        f'{self.title} {time.strftime("%Y/%m/%d %H:%M")} (UTC+{adjust_hr})'
                    )
                    if (count == 0) & (it == 0):
                        cbar = fig.colorbar(ax1, ax=ax, orientation="horizontal")
                        cbar.ax.set_title(
                            f'{var.description.title()} ({degree+"C" if self.variable=="tc2" else var.units})'
                        )
                    else:
                        cbar.update_normal(ax1)
                    fig.savefig(outfile, bbox_inches="tight")
                    ax1 = None
                    del (fig, ax)

        if self.nctype == "wrf":
            draw_wrf(self, adjust_hr, **kwargs)
        if self.nctype == "cdowrf":
            draw_cdowrf(self, adjust_hr, **kwargs)
        self.ax.clear()

    def create_video(self, string, framerate=4, out="out.mp4", overwrite_output=True):
        output_options = {
            "crf": 20,
            "vf": "scale=trunc(iw/2)*2:trunc(ih/2)*2",
            "vcodec": "libx264",
            "pix_fmt": "yuv420p",
        }
        (
            ffmpeg.input(string, pattern_type="glob", framerate=framerate)
            .output(out, **output_options)
            .run(overwrite_output=overwrite_output)
        )
