install-docker:
	docker build -t mpc/nextqcflow-python:latest -f docker/python/Dockerfile .
	docker build -t mpc/nextqcflow-comet:latest -f docker/comet/Dockerfile .
	docker pull quay.io/biocontainers/pia:1.4.8--hdfd78af_0
	docker pull quay.io/biocontainers/thermorawfileparser:1.4.3--ha8f3691_0
	docker pull mfreitas/tdf2mzml