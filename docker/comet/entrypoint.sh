#!/usr/bin/env bash

# Ensures that the user running the container has the binary folder in PATH
export PATH=/usr/local/bin:$PATH

set -e

exec "$@"
