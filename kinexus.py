#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 15:10:08 2019

@author: igocer
"""

import requests
import urllib.request
import time
from bs4 import BeautifulSoup
import random
from time import sleep
import numpy as np
import pandas as pd

url = 'http://www.phosphonet.ca/?search=Q9UBK2'
response = requests.get(url)
soup = BeautifulSoup(response.text, "html.parser")


phospho_sites = [ x.get_text() for x in soup.findAll("td", {"class": "pSiteNameCol"})]


site_base_url = "http://www.phosphonet.ca/kinasepredictor.aspx?uni=Q9UBK2&ps="
site_tables = []


for idx, site in enumerate(phospho_sites):
    site_url = site_base_url + site
    response = requests.get(site_url)
    soup = BeautifulSoup(response.text, "html.parser")
    site_tables.append(soup)
    print(site_url)
    sleep(random.uniform(2, 5))
    if idx%30 == 0:
        sleep(random.uniform(30, 40))

PGC_sites = pd.DataFrame()

for idx, site in enumerate(site_tables):
#test = site_tables[0]
    strings = list(site.html.stripped_strings)
    start_index = strings.index("Kinase 1:")
    end_index = 50*7+49
    string_array = np.array(strings[start_index:end_index])
    string_array = string_array.reshape((50,7), order='C')
    
    df = pd.DataFrame(string_array)
    df.drop([0,4,5], axis=1, inplace=True)
    df.reset_index(inplace=True)
    df.columns = ["kinase_no", "kinase_name", "uniprot_id", "kinexus_score", "kinexus_score_v2"]
    df.insert(0, "site", phospho_sites[idx][1:])
    df.insert(0, "aa", phospho_sites[idx][0])
    if(PGC_sites.empty):
       PGC_sites = df 
    else:
        PGC_sites = pd.concat([PGC_sites, df])

PGC_sites['site'] = PGC_sites['site'].astype(dtype='int32')
PGC_sites['kinexus_score'] = PGC_sites['kinexus_score'].astype(dtype='int32')        
PGC_sites['kinexus_score_v2'] = PGC_sites['kinexus_score_v2'].astype(dtype='int32')
PGC_sites.reset_index(drop=True)

PGC_sites['msite_shift'] = np.where(PGC_sites['site'] < 170, 2, 1)
PGC_sites['msite'] = PGC_sites['site'] - PGC_sites['msite_shift']


PGC_sites.query('kinexus_score>700')


PGC_sites.to_csv("PPARGC1A_kinexus.csv", sep='\t', index=False)
