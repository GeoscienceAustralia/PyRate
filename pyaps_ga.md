## How to use PyAPS with PyRate

### Create a `ECMWF` directory inside the `PyRate` main

    cd PyRate
    mkdir ECMWF

This is where the `.grb` will be saved by the `pyaps`.

### Clone/copy the `pyaps`

Clone the [pyaps package from this Geoscience Austalia github repo](https://github.com/GeoscienceAustralia/PyAPS.git) inside the `PyRate` folder. 


### Follow these instructions from 'GIAnT_doc.pdf' so you can download the .grb files

```ERA-Interim Re-Analysis products are downloaded from the ECMWF repository in Europe 3 . We use the variables Temperature, Geopotential Height and Relative Humidity (default) at each Pressure Levels. If you are working on a machine with a non-US IP address, you should use this option. You need to register on the ECMWF website and provide your password in the file GIAnT/pyaps/model.cfg (a template is provided). Follow these steps to set up your access```

    1. Agree to terms and conditions
    To get your personalized access key from ECMWF, by read and
    agree to this license http://data-portal.ecmwf.int/data/d/
    license/interim_full/.
    2. Register
    Register at https://apps.ecmwf.int/registration/ with the
    same email address used to agree to the terms and conditions.
    3. Login
    Login at https://apps.ecmwf.int/auth/login/.
    4. Sign the agreement
    Sign the agreement at http://apps.ecmwf.int/datasets/licences/
    general/.
    5. Copy your key to model.cfg Obtain your API key from https:
    //api.ecmwf.int/v1/key/ and copy it to your model.cfg file.
    PyAPS only downloads the fields of interest from the ECMWF archive
    - Humidity and Temperature as a function of pressure level. Each file
    (global for single time eopch) is about 28MB in size.

SB: Copy the `pyaps/model.cfg.template` as `model.cfg` and fill in with your credentials.    

### Run the `remove_aps_delay.py` from the 'PyRate' directory

    cd PyRate
    python remove_aps_delay.py