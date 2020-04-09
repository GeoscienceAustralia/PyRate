YEARS_PER_DAY = 1 / 365.25

PYRATE_DEFAULT_CONFIGRATION = {
    "obsdir": {
        "DataType": "path",
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "ifgfilelist": {
        "DataType": "path",
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "demfile": {
        "DataType": "path",
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "demHeaderFile": {
        "DataType": "path",
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "slcFileDir": {
        "DataType": "path",
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "slcfilelist": {
        "DataType": "path",
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "cohfiledir": {
        "DataType": "path",
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "cohfilelist": {
        "DataType": "path",
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "outdir": {
        "DataType": "path",
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "processor": {
        "DataType": int,
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [0, 1],
        "Required": False
    },
    "noDataAveragingThreshold": {
        "DataType": float,
        "DefaultValue": 0.0,
        "MinValue": 0,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "noDataValue": {
        "DataType": float,
        "DefaultValue": 0.0,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "nan_conversion": {
        "DataType": int,
        "DefaultValue": 0,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [0, 1],
        "Required": False
    },
    "parallel": {
        "DataType": int,
        "DefaultValue": 0,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [0, 1, 2],
        "Required": False
    },
    "processes": {
        "DataType": int,
        "DefaultValue": 8,
        "MinValue": 1,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "cohmask": {
        "DataType": int,
        "DefaultValue": 0,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [0, 1],
        "Required": False
    },
    "cohthresh": {
        "DataType": float,
        "DefaultValue": None,
        "MinValue": 0.0,
        "MaxValue": 1.0,
        "PossibleValues": None,
        "Required": False
    },
    "ifgcropopt": {
        "DataType": int,
        "DefaultValue": 1,
        "MinValue": 1,
        "MaxValue": 4,
        "PossibleValues": [1, 2, 3, 4],
        "Required": False
    },
    "ifglksx": {
        "DataType": int,
        "DefaultValue": 1,
        "MinValue": 1,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "ifglksy": {
        "DataType": int,
        "DefaultValue": 1,
        "MinValue": 1,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "ifgxfirst": {
        "DataType": float,
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "ifgxlast": {
        "DataType": float,
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "ifgyfirst": {
        "DataType": float,
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "ifgylast": {
        "DataType": float,
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "refx": {
        "DataType": float,
        "DefaultValue": -1,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "refy": {
        "DataType": float,
        "DefaultValue": -1,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "refnx": {
        "DataType": int,
        "DefaultValue": 10,
        "MinValue": 1,
        "MaxValue": 50,
        "PossibleValues": None,
        "Required": False
    },
    "refny": {
        "DataType": int,
        "DefaultValue": 10,
        "MinValue": 1,
        "MaxValue": 50,
        "PossibleValues": None,
        "Required": False
    },
    "refchipsize": {
        "DataType": int,
        "DefaultValue": 21,
        "MinValue": 1,
        "MaxValue": 101,
        "PossibleValues": None,
        "Note": "Must be an odd number.",
        "Required": False
    },
    "refminfrac": {
        "DataType": float,
        "DefaultValue": 0.5,
        "MinValue": 0.0,
        "MaxValue": 1.0,
        "PossibleValues": None,
        "Required": False
    },
    "refest": {
        "DataType": int,
        "DefaultValue": 1,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [1, 2],
        "Required": False
    },
    "orbfit": {
        "DataType": int,
        "DefaultValue": 0,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [0, 1],
        "Required": False
    },
    "orbfitmethod": {
        "DataType": int,
        "DefaultValue": 2,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [1, 2],
        "Required": False
    },
    "orbfitdegrees": {
        "DataType": int,
        "DefaultValue": 1,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [1, 2, 3],
        "Required": False
    },
    "orbfitlksx": {
        "DataType": int,
        "DefaultValue": 10,
        "MinValue": 1,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "orbfitlksy": {
        "DataType": int,
        "DefaultValue": 10,
        "MinValue": 1,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "apsest": {
        "DataType": int,
        "DefaultValue": 0,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [0, 1],
        "Required": False
    },
    "slpfmethod": {
        "DataType": int,
        "DefaultValue": 1,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [1, 2],
        "Required": False
    },
    "slpfcutoff": {
        "DataType": float,
        "DefaultValue": 1.0,
        "MinValue": 0.001,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "slpforder": {
        "DataType": int,
        "DefaultValue": 1,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [1, 2, 3],
        "Required": False
    },
    "slpnanfill": {
        "DataType": int,
        "DefaultValue": 0,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [0, 1],
        "Required": False
    },
    "slpnanfill_method": {
        "DataType": str,
        "DefaultValue": "cubic",
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": ["linear", "nearest", "cubic"],
        "Required": False
    },
    "tlpfmethod": {
        "DataType": int,
        "DefaultValue": 1,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [1, 2, 3],
        "Required": False
    },
    "tlpfcutoff": {
        "DataType": float,
        "DefaultValue": 1.0,
        "MinValue": YEARS_PER_DAY,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "tlpfpthr": {
        "DataType": int,
        "DefaultValue": 1,
        "MinValue": 1,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "tscal": {
        "DataType": int,
        "DefaultValue": 0,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [0, 1],
        "Required": False
    },
    "tsmethod": {
        "DataType": int,
        "DefaultValue": 2,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [1, 2],
        "Required": False
    },
    "smorder": {
        "DataType": int,
        "DefaultValue": None,
        "MinValue": None,
        "MaxValue": None,
        "PossibleValues": [1, 2],
        "Required": False
    },
    "smfactor": {
        "DataType": float,
        "DefaultValue": -1.0,
        "MinValue": -5.0,
        "MaxValue": 0,
        "PossibleValues": None,
        "Required": False
    },
    "ts_pthr": {
        "DataType": int,
        "DefaultValue": 3,
        "MinValue": 1,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "nsig": {
        "DataType": int,
        "DefaultValue": 2,
        "MinValue": 1,
        "MaxValue": 10,
        "PossibleValues": None,
        "Required": False
    },
    "pthr": {
        "DataType": int,
        "DefaultValue": 3,
        "MinValue": 1,
        "MaxValue": None,
        "PossibleValues": None,
        "Required": False
    },
    "maxsig": {
        "DataType": int,
        "DefaultValue": 10,
        "MinValue": 0,
        "MaxValue": 1000,
        "PossibleValues": None,
        "Required": False
    }
}
