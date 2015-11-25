__author__ = 'sudipta'
"""
This script is used to delete files that are not used in any where in the package.

Inputs:

The directory to delete files from: dir_name_to_delete_files_from
The directory to find usage: dir_name_to_find_usage

All files of 'dir_name_to_delete_files_from' not used in 'dir_name_to_find_usage' are deleted.

The script accounts for file that are referenced inside the files in  'dir_name_to_delete_files_from'.

Usage:
just change the names of 'dir_name_to_delete_files_from' and `dir_name_to_find_usage` to the desired dirs and run the
script
"""

import os
import subprocess
import shlex


def list_of_all_files(dir):
    all_files = []
    for root, dirs, files in os.walk(dir):
        for file in files:
            all_files.append(os.path.join(root, file))

    return all_files


def delete_file(file_name):
    try:
        os.remove(file_name)
    except OSError:
        raise
        pass


def clean_files_from_this_directory_to_that(all_files, dir_name_to_find_usage):

    used_files_ = []
    potentially_deleted_files_ = []

    for f in all_files:
        command = shlex.split("grep -R -H {file_name} {dir}".format(
            file_name=os.path.split(f)[-1], dir=dir_name_to_find_usage))
        try:
            output = subprocess.check_output(command)
            used_files_.append(f)
        except:
            potentially_deleted_files_.append(f)

    print '++'*50
    return used_files_, potentially_deleted_files_


if __name__ == '__main__':
    dir_name_to_delete_files_from = "/home/sudipta/Dropbox/GA/PyRate/tests/"
    dir_name_to_find_usage = "/home/sudipta/Dropbox/GA/PyRate/pyrate/"
    all_files = list_of_all_files(dir_name_to_delete_files_from)

    used_files, _ = clean_files_from_this_directory_to_that(all_files, dir_name_to_find_usage)
    print "=="*50
    need_these_files_back, _ = clean_files_from_this_directory_to_that(all_files, dir_name_to_delete_files_from)

    used_files_filenames = map(lambda uf: os.path.split(uf)[-1], used_files)
    copied_files_filenames = map(lambda uf: os.path.split(uf)[-1], need_these_files_back)

    set_of_files_to_keep = set(used_files_filenames + copied_files_filenames)

    count = 0
    count_deleted = 0
    for f in all_files:
        if os.path.split(f)[-1] in set_of_files_to_keep:
            count += 1
        else:
            delete_file(f)
            count_deleted += 1

    print 'files deleted', count_deleted
    print 'files saved', count
