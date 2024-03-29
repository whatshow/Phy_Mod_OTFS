import os
import shutil

# version control
version = "1.0.2";
install_requires = ['numpy>=1.23.5'];

# path
path_cur = os.getcwd();
path_dist = "_dist/";
path_dist_pkg = path_dist + "whatshow_phy_mod_otfs/";
# file
file = "OTFS.py"
file_init = "__init__.py";
file_setup = "setup.py";
file_readme = "README.md";
# file content
file_init_content = "from .OTFS import OTFS";
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
      install_requires={install_requires},\n\
      long_description = description,\n\
      long_description_content_type = \"text/markdown\"\n\
);';


# clear the distribution
if os.path.exists(path_dist):
    shutil.rmtree(path_dist);

# create the distribution folder
os.makedirs(path_dist_pkg);   

# copy files to the distribution
shutil.copyfile(file, path_dist_pkg + file);

# add config files
# __init__.py
with open(path_dist_pkg + file_init, "w") as init:
    init.write(file_init_content);
# setup.py
with open(path_dist + file_setup, "w") as setup:
    setup.write(file_setup_content);

# README.md
description_pypi = "";
with open(file_readme, "r") as readme:
    description = readme.read(); 
    txt_1_start = description.find("## How to install");
    txt_2_start = description.find("## How to use");
    description_pypi = description[:txt_1_start] + description[txt_2_start:];
with open(path_dist + file_readme, "w") as readme_pypi:
    readme_pypi.write(description_pypi);

# build & upload
os.chdir(path_dist);
os.system("python " + file_setup + " sdist bdist_wheel");
os.system("twine upload dist/*");
os.chdir(path_cur);
