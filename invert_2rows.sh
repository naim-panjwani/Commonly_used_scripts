#!/bin/bash

filename=$1
paste <(head -1 $filename |tr '\t' '\n') <(head -2 $filename |tail -1 |tr '\t' '\n')
