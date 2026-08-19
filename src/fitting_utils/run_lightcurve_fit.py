import sys
import os
import shutil
import glob
import Tiberius
from global_utils import parseInput

def make_folder(folder_name):
	if os.path.isdir(folder_name):
		pass
	else:
		os.system("mkdir %s"%folder_name)

def move_input_files_to_output_folder(input_dict):
    output_folder = input_dict['output_foldername']

    # Files listed in fitting_input.txt
    keys = [
        'prior_filename',
        'model_input_files',
        'GP_model_input_files',
    ]

    # Make sure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Also move fitting_input.txt itself
    files_to_move = ['fitting_input.txt']

    # Add all files from the relevant entries
    for key in keys:
        if key in input_dict:
            files = input_dict[key].split(';')
            files_to_move.extend(files)

    for filename in files_to_move:
        filename = filename.strip()

        if not filename:
            continue

        # If only a filename is given, look in the current directory
        source = filename

        # If it is already in the output folder, don't do anything
        destination = os.path.join(output_folder, os.path.basename(filename))

        # Already at the destination
        if os.path.abspath(source) == os.path.abspath(destination):
            # print(f"Already in output folder: {filename}")
            continue

        # File exists in current location
        if os.path.isfile(source):
            # print(f"Moving {filename} -> {destination}")
            shutil.copy(source, destination)

        # # File doesn't exist at source, but already exists in output folder
        # elif os.path.isfile(destination):
        #     print(f"Already in output folder: {os.path.basename(filename)}")

        # else:
        #     print(f"WARNING: File not found: {filename}")

Tiberius_path = "/".join(sys.argv[0].split("/")[:-1])
starting_bin = int(sys.argv[1])
stopping_bin = int(sys.argv[2])

input_dict = parseInput('fitting_input.txt')

make_folder("%s"%input_dict['output_foldername'])
make_folder("%s/tables"%input_dict['output_foldername'])
make_folder("%s/plots"%input_dict['output_foldername'])
make_folder("%s/pickled_objects"%input_dict['output_foldername'])

# Move required input files into output folder
move_input_files_to_output_folder(input_dict)

if bool(int(input_dict['generate_LDCs'])):
	os.system("python %s/generate_LDCs.py"%(Tiberius_path))

print('Running light_curve_fit.py..')
for i in range(starting_bin,stopping_bin):
	os.system("python %s/light_curve_fit.py %d"%(Tiberius_path,i))

print('Running plot_output.py to generate plots and tables...')
os.system("python %s/plot_output.py -s -st -cp"%Tiberius_path)

# os.system(f"cp *.txt {input_dict['output_foldername']}")
# # # # # os.system("python %s/model_table_generator.py"%Tiberius_path)

# exclude_keywords = ["fitting_input", "LD_coefficients", "model", "prior"]

# for f in glob.glob("*.txt"):

#     if not any(k in f for k in exclude_keywords):

#         dst = f"{input_dict['output_foldername']}/tables/{f}"

#         if os.path.exists(dst):
#             os.remove(dst)

#         shutil.move(f, dst)

# os.system("mv *.pickle %s/pickled_objects/"%input_dict['output_foldername'])
# os.system("mv *.png %s/plots/"%input_dict['output_foldername'])
# os.system("mv *.pdf %s/plots/"%input_dict['output_foldername'])
