import os
requirment_file_path = "../requirements.txt"

print([line.strip() for line in open(requirment_file_path, "r").readlines() if len(line)>2])


