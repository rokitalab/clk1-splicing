#!/bin/sh
#
# Script to remove results from all analysis modules
# Author: Patricia Sullivan
# Date: July 2025
#


set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Get base directory of project
cd ..

# Remove all contents inside analyses/*/plots/
for dir in analyses/*/plots/; do
  [ -d "$dir" ] && rm -rfv "$dir"/*
done

# Remove all contents inside analyses/*/results/
for dir in analyses/*/results/; do
  [ -d "$dir" ] && rm -rfv "$dir"/*
done
