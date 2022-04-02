#!/bin/bash

input=$1
output=$2

python -c "import sys; print('\n'.join('\t'.join(c) for c in zip(*(l.strip().split() for l in sys.stdin.readlines() if l.strip()))))" < "$input" > "$output"
