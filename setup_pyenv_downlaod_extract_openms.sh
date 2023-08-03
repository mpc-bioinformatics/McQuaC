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

# Download Dinosaur
wget  -O bin/Dinosaur.jar https://github.com/fickludd/dinosaur/releases/download/1.2.0/Dinosaur-1.2.0.free.jar

# Install environment
rm -rf Pipfile
pipenv --rm | true
pipenv install --python 3.9
pipenv run pip install pyopenms==2.7.0
pipenv run pip install pandas
pipenv run pip install sqlalchemy
pipenv run pip install plotly
pipenv run pip install scikit-learn
pipenv run pip install numpy
pipenv run pip install fisher-py==1.0.21

