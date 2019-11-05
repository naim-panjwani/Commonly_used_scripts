#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 18:21:59 2019

@author: naim
"""

import string
import pandas as pd
import numpy as np
from tqdm import tqdm

import sqlalchemy as sa
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import Session
from sqlalchemy import create_engine, inspect, String, Integer
import pymysql
pymysql.install_as_MySQLdb()


def getSNPPositions(snp_list, build=38):
    try:
        if build == 19 or build == 37:
            buildstr = 'hg19'
        else:
            buildstr = 'hg38'
    except:
        ValueError('Please enter build in integer format')
    pos = []
    allowed_chars = 'rs' + string.digits
    engine = sa.create_engine(f'mysql://genome@genome-mysql.cse.ucsc.edu:3306/{buildstr}')
    for snp in tqdm(snp_list):
        querysnp = snp.split(';')[0]
        if any([char not in allowed_chars for char in querysnp]):
            pos.append(-1)
        else:
            pos.append(engine.execute(f"select chromEnd from snp151 where name='{querysnp}'").fetchall()[0][0])
    return pos


# Example
snp_list = ['rs662702', 'rs7512462']
getSNPPositions(snp_list, 38)

