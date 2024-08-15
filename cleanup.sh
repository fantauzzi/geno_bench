#!/usr/bash

set -e

rm -rf *.fna.*
rm *.fna
rm ncbi_dataset.zip
dir_name=$(find . -maxdepth 1 -type d ! -name "." -print -quit)
rm -rf "${dir_name##*/}"*

