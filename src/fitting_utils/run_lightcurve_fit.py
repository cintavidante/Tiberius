import sys
import os
import shutil
import glob
import Tiberius
from global_utils import parseInput

Tiberius_path = "/".join(sys.argv[0].split("/")[:-1])
starting_bin = int(sys.argv[1])
stopping_bin = int(sys.argv[2])

input_dict = parseInput('fitting_input.txt')

def make_folder(folder_name):
	if os.path.isdir(folder_name):
		pass
	else:
		os.system("mkdir %s"%folder_name)

make_folder("%s"%input_dict['output_foldername'])
make_folder("%s/tables"%input_dict['output_foldername'])
make_folder("%s/plots"%input_dict['output_foldername'])
make_folder("%s/pickled_objects"%input_dict['output_foldername'])

if bool(int(input_dict['generate_LDCs'])):
	os.system("python %s/generate_LDCs.py"%(Tiberius_path))

print('Running light_curve_fit.py..')
for i in range(starting_bin,stopping_bin):
	os.system("python %s/light_curve_fit.py %d"%(Tiberius_path,i))

# print('Running plot_output.py to generate plots and tables...')
# os.system("python %s/plot_output.py -s -st -cp"%Tiberius_path)

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
