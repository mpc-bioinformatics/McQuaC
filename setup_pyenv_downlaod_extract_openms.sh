#!/bin/bash

set -e

# Delete previously downloaded openms
rm -rf bin/openms/

# Download openms
wget -O openms.deb https://github.com/OpenMS/OpenMS/releases/download/Release2.8.0/OpenMS-2.8.0-Debian-Linux-x86_64.deb

# Extract only the data-part of the deb package
ar x openms.deb data.tar.gz

# Extract the data into the openms_folder
mkdir -p ./bin/openms

# Delete downloaded and extracted deb package
tar -xf data.tar.gz -C ./bin/openms
rm -rf openms.deb data.tar.gz

# Download Pia:
rm -rf bin/pia/
wget -O pia.zip https://github.com/mpc-bioinformatics/pia/releases/download/1.4.7/pia-1.4.7.zip
mkdir -p ./bin/pia
unzip -d ./bin/pia pia.zip
rm -rf pia.zip


# Download Dinosaur
wget  -O bin/Dinosaur.jar https://github.com/fickludd/dinosaur/releases/download/1.2.0/Dinosaur-1.2.0.free.jar

# Install environment
python -m pip install pyopenms
python -m pip install pandas
python -m pip install sqlalchemy
python -m pip install plotly
python -m pip install scikit-learn
python -m pip install numpy
python -m pip install fisher-py==1.0.21

