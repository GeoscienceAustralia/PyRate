import pathlib


def delete_tsincr_files(params):

    out_dir = pathlib.Path(params["outdir"])
    for file_path in out_dir.iterdir():
        if "tsincr" in str(file_path):
            file_path.unlink()
