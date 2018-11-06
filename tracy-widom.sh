#!/bin/bash
# Runs Tracy Widom test given a file of eigenvalues, one per line
# Syntax: bash tracy-widom <eigenvalues_file> <outfile>
# Example: bash tracy-widom eigenvalues.txt cts_tracy_widom.txt

ln -s /usr/local/Local_tools/EIG3.0/POPGEN/twtable
twstats -t twtable -i $1 -o $2
