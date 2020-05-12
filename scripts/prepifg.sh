#!/usr/bin/env bash
# inputs are extents, thresh, coherence file, coherence_threshold
# coherence masking
# gdal average
# gdal average has two steps
# create a raster with the input raster nan_matrix, then compute the nan_frac of the input raster
# take average of the input raster, and set nan_frac > thresh = nan

inputdir=/home/sudipta/repos/PyRate/out/gamma/out
outdir_crop=/home/sudipta/repos/PyRate/out/gamma/out/cropped
outdir_coh=/home/sudipta/repos/PyRate/out/gamma/out/coherence

mkdir -p ${outdir_crop}
mkdir -p ${outdir_coh}

function crop_resample_average {
    f=$1
    echo ${f}
    gdalwarp -overwrite ${f}
}

function nan_matrix {
    f=$1
    echo ${f}
    python gdal_calc_local.py --overwrite ${f}
}

function resampled_average {
    f=$1
    echo ${f}
    python gdal_calc_local.py --overwrite ${f}
}

export -f crop_resample_average
#export -f mask
#export -f multilook
export -f nan_matrix
export -f resampled_average

#ls ${inputdir}/*.tif | parallel crop ${outdir_crop}

# while IFS= read -r lineA && IFS= read -r lineB <&3; do
#  echo "$lineA"; echo "$lineB"
# done <fileA 3<fileB


# outfile.txt is tab/space delimited file with `ifg - coh mask` pairs in each line

# 1. apply coh mask at the very beginning
#cat mask.txt | parallel mask

# 2. crop, average, resample
#cat crop.txt | parallel crop_resample_average

# 3. create nan fractions
#cat nan_fraction.txt | parallel nan_matrix

# 4. crop average resample nan fractions
# cat crop_average_nan_fraction.txt | parallel crop_resample_average

# 5. now resampled
# cat resampled_average.txt | parallel resampled_average

