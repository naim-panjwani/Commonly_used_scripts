#!/bin/bash
module load ncdu/1.14
ncdu -1xo- fanwang/ |gzip >fanwang_ncdu_usage.gz

