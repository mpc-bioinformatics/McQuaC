install-docker:
	docker build -t mpc/nextqcflow-python:latest -f docker/python/Dockerfile .
	docker pull quay.io/medbioinf/comet-ms:v2024.01.0
	docker pull quay.io/biocontainers/pia:1.5.0--hdfd78af_0
	docker pull quay.io/biocontainers/thermorawfileparser:1.4.3--ha8f3691_0
	docker pull mfreitas/tdf2mzml

comet-params:
	docker run --rm -it quay.io/medbioinf/comet-ms:v2024.01.0 comet -p > new.comet.params