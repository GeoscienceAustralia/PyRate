# -*- coding: utf-8 -*-

import os, sys
from pyrate.shared import Ifg
import pyrate.config as cfg
from pyrate.prepifg import (
    CUSTOM_CROP,
    extents_from_params,
    prepare_ifgs
)

def main():
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [OPTIONS] [config-file]')
    parser.add_option('-t', '--thresh', type=float,
                        help='Nodata averaging threshold')
    parser.add_option('-i', '--ifglist', type=str,
                        help='name of file containing list of interferograms')
    parser.add_option('-v', '--verbose', action='store_true', default=False,
                        help='Display more output')
    options, args = parser.parse_args()

    try:
        cfg_path = args[0] if args else 'pyrate.conf'
        pars = cfg.get_config_params(cfg_path)
    except IOError as err:
        errmsg = 'Config file error: %s "%s"'
        print(errmsg % (err.strerror, err.filename))
        sys.exit(err.errno)

    ifgListFile = options.ifglist or pars.get(cfg.IFG_FILE_LIST)
    if ifgListFile is None:
        print 'Interferogram file name not provided'
        sys.exit(err.errno)

    ifglist = cfg.parse_namelist(ifgListFile)
    obs_dir = pars.get(cfg.OBS_DIR) 
    if obs_dir is not None:
        ifg_paths = [os.path.join(obs_dir, p) for p in ifglist]
    else:
        ifg_paths = ifglist
    ifgs = [Ifg(p) for p in ifg_paths]

    crop = pars[cfg.IFG_CROP_OPT]
    exts = extents_from_params(pars) if crop == CUSTOM_CROP else None

    kwargs = { 'thresh': options.thresh if options.thresh else 0.5,
                'user_exts': exts, 'verbose': options.verbose, }

    prepare_ifgs(ifgs, crop_opt=crop, xlooks=pars[cfg.IFG_LKSX],
                ylooks=pars[cfg.IFG_LKSY], **kwargs)


if __name__ == '__main__':
    main()

