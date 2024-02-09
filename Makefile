install-native:
	# Install dependencies
	apt-get install -y bash mono-complete openjdk-17-jre python3 python3-pip git curl make build-essential libssl-dev zlib1g-dev \
   		libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm python-is-python3 unzip \
   		libncursesw5-dev xz-utils tk-dev libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev \
   		libboost-date-time-dev libboost-iostreams-dev libboost-regex-dev libboost-math-dev \
   		libboost-random-dev zlib1g libbz2-dev libsvm3 libxerces-c-dev libglpk-dev libqt5network5 libqt5opengl5 libqt5svg5 libqt5webkit5 libqt5core5a \
   		libqt5sql5

	# OpenMS
	rm -rf ./bin/openms
	curl -L -o openms.deb https://github.com/OpenMS/OpenMS/releases/download/Release2.8.0/OpenMS-2.8.0-Debian-Linux-x86_64.deb
	ar x openms.deb data.tar.gz
	mkdir -p ./bin/openms
	tar -xf data.tar.gz -C ./bin/openms
	rm -rf openms.deb data.tar.gz

	# Comet
	rm -f ./bin/comet
	curl -L -o ./bin/comet https://github.com/UWPR/Comet/releases/download/v2022.01.2/comet.linux.exe 
	chmod +x ./bin/comet

	# Pia
	rm -rf bin/pia
	curl -L -o pia.zip https://github.com/medbioinf/pia/releases/download/1.4.8/pia-1.4.8.zip
	mkdir -p ./bin/pia
	unzip -d ./bin/pia pia.zip
	rm -rf pia.zip

	# thermorawfileparser
	rm -rf bin/thermorawfileparser
	curl -L -o thermorawfileparser.zip https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.3/ThermoRawFileParser1.4.3.zip
	mkdir -p ./bin/thermorawfileparser
	unzip -d ./bin/thermorawfileparser thermorawfileparser.zip
	rm -rf thermorawfileparser.zip

	# Python requirements
	pip install -r requirements.txt
	pip install -r nativ/requirements.txt

	# Copy wrapper files
	cp native/!(*.txt) ./bin


install-docker:
	docker build -t mpc/nextqcflow-python:latest -f docker/python/Dockerfile .
	docker build -t mpc/nextqcflow-comet:latest -f docker/comet/Dockerfile .
	docker pull quay.io/biocontainers/pia:1.4.8--hdfd78af_0
	docker pull quay.io/biocontainers/thermorawfileparser:1.4.3--ha8f3691_0
	docker pull mfreitas/tdf2mzml