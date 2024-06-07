import os
import shutil

# version control
version = "2.1.12";

# path
path_cur = os.getcwd();
path_dist = "_dist/";
path_dist_pkg = path_dist + "whatshow_phy_mod_otfs/";
# file
files = ["OTFS.py", "OTFSResGrid.py", "OTFSDetector.py"];
file_init = "__init__.py";
file_setup = "setup.py";
file_readme = "README.md";
# file content
file_init_content = "from .OTFS import OTFS\nfrom .OTFSResGrid import OTFSResGrid\nfrom .OTFSDetector import OTFSDetector";
file_setup_content = f'from setuptools import setup, find_packages\n\
\n\
# prepare the instruction\n\
description = "";\n\
with open("README.md", "r") as readme:\n\
    description = readme.read();\n\
setup(\n\
      name="whatshow_phy_mod_otfs",\n\
      version="{version}",\n\
      packages=find_packages(),\n\
      long_description = description,\n\
      long_description_content_type = \"text/markdown\"\n\
);';


# clear the distribution
if os.path.exists(path_dist):
    shutil.rmtree(path_dist);

# create the distribution folder
os.makedirs(path_dist_pkg);   

# copy files to the distribution
for file in files:
    if os.path.exists(file):
        shutil.copyfile(file, path_dist_pkg + file);

# add config files
# __init__.py
shutil.copyfile(file_init, path_dist_pkg + file_init);
# setup.py
with open(path_dist + file_setup, "w") as setup:
    setup.write(file_setup_content);

# README.md
description_pypi = "";
with open(file_readme, "r") as readme:
    description = readme.read(); 
    txt_1_start = description.find("## How to install");
    # txt_2_start = description.find("## How to use");
    # txt_3_start = description.find("## Samples");
    # description_pypi = description[:txt_1_start] + description[txt_2_start:txt_3_start];
    description_pypi = description[:txt_1_start];
with open(path_dist + file_readme, "w") as readme_pypi:
    readme_pypi.write(description_pypi);

# build & upload
os.chdir(path_dist);
os.system("python " + file_setup + " sdist bdist_wheel");
os.system("twine upload dist/*");
os.chdir(path_cur);
