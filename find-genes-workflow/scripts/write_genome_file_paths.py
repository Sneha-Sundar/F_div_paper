import argparse
import sys

parser = argparse.ArgumentParser(prog='Write genome file paths to a single file for use in ANI calculation',
                    description='Write genome file paths to a single file for use in ANI calculation.')
parser.add_argument('-i','--input_file_paths',help="List of genome paths",required = True,nargs='+') 
parser.add_argument('-o','--output_file',help="Path to output file",required = True)
args = parser.parse_args()

file_paths = args.input_file_paths
output_file = args.output_file

with open(output_file, "w") as f:
    for path in file_paths:
        f.write(path + "\n")
