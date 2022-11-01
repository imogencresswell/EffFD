#!/bin/bash

set -e
set -u
set -o pipefail

pycodestyle *.py
python3 test_utils.py
bash test_file_generation.sh
