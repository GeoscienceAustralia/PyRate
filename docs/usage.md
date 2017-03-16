# Usage

## PyRate workflow

`PyRate` installs an executable `pyrate`.

Use `help` for the different command line options:

    >> pyrate --help
    Usage: pyrate [OPTIONS] COMMAND [ARGS]...

    Options:
      -v, --verbosity [DEBUG|INFO|WARNING|ERROR]
                                      Level of logging
      --help                          Show this message and exit.

    Commands:
      linrate
      postprocess
      prepifg

As you can see `pyrate` has three command line options.

1. prepifg
1. linrate
1. postprocess

Below we discuss these options.

### Preparing interferograms: prepifg
The first step of PyRate is to convert the unwrapped interferograms into
geotiffs, followed by multilooking and cropping. Both of these operations are
performed by `pyrate prepifg` command:

    >> pyrate prepifg --help
    Usage: pyrate prepifg [OPTIONS] CONFIG_FILE

    Options:
      --help  Show this message and exit.

So one can use the `prepfig` command as the following:

    pyrate prepifg /path/to/config_file

The two major steps during the `prepifg` operation are described below.

#### Data formatting: convert to geotiff

The `prepifg` command will determine the input format from the value specified
at the *processor:* keyword in the config file (0: ROI\_PAC; 1: GAMMA)

A GAMMA translation requires a geographic DEM header file (\*.dem.par) and SLC parameter files (\*.slc.par) for both master and slave images to extract metadata required for the formatting. Therefore three header files are needed to format each geocoded unwrapped GAMMA interferogram <GAMMA_FILE>. The path and name of the DEM header file are specified in the config file under the *demHeaderFile:* keyword. The SLC parameter files should be in the same location as the interferogram file and are found automatically by date string pattern matching.

A ROI\_PAC translation requires a header/resource file (*.rsc* extension) for the geocoded unwrapped ROI_PAC interferogram (in the same directory) and either the geographic projection (e.g. 'WGS84') specified as an option or a header/resource file for the geographic DEM containing the geographic projection in the parameter DATUM:

#### Image transformations: multilooking and cropping
This `prepifg` command will also perform multi-looking (resampling) and
cropping the images.

Two examples of the config files are provided in `configs` directory,
with examples of the `roipac` and `gamma` prepifg configuration.
Both config files can be used with `prepifg`.

### Linear rate and time series analysis: linrate

    >> pyrate linrate --help
    Usage: pyrate linrate [OPTIONS] CONFIG_FILE

    Options:
      -r, --rows INTEGER  divide ifgs into this many rows
      -c, --cols INTEGER  divide ifgs into this many columns
      --help              Show this message and exit


This is the core of the processing tools, handled by the `linrate` command:

    pyrate linrate path/to/config_file -c 3 -r 4

This command will does the time series and linear rate analysis, but has the
options to break the interferograms into tiles of `r` rows and `c` columns.
So this above command will break the interferograms into 12 tiles and will
produce 12 linear rate and time series predictions corresponding to each tile.

The optional rows and columns help us create smaller `tiles` of the
interferograms that can be accommodated in the memory. The number of tiles
chosen should be as small as possible that fits in the system memory.

### Putting it back together: postprocess
The last step in `pyrate` is to put all the tiles back together from the
`linrate` part.

    >> pyrate postprocess --help
    Usage: pyrate postprocess [OPTIONS] CONFIG_FILE

    Options:
      -r, --rows INTEGER  divide ifgs into this many rows
      -c, --cols INTEGER  divide ifgs into this many columns
      --help              Show this message and exit.

Make sure to use the same number of rows and columns with `postprocess` as
with `linrate`:

    pyrate postprocess path/to/config_file -c 3 -r 4

## MPI Support
`PyRate` has been designed for supercomputers and instructions to use an HPC
 cluster can be found in the [pbs directory](pbs).

## Python multiprocessing support
In addition to the MPI support for HPC, `PyRate` can be use standard
multiprocessing simply by turning `parallel:  1` in the config file to take
advantage of multiple cores on a single PC.