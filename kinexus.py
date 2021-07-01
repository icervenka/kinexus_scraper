#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 15:10:08 2019

@author: icervenka
"""

import argparse
import requests
from bs4 import BeautifulSoup
import random
from time import sleep
import numpy as np
import pandas as pd

# -------------------------
# batch and sleep constants
# -------------------------
sleep_individual_low = 2
sleep_individual_high = 5

batch_size = 30
sleep_batch_low = 30
sleep_batch_high = 40

# -------------------------
# parser
# -------------------------
parser = argparse.ArgumentParser(description='Get protein phosphosites from phosphonet.ca')
parser.add_argument('ids', metavar='uniprot_id', type=str, nargs='+',
                    help='Uniprot accession numbers of protiens to retrieve')
parser.add_argument('-o', '--outdir', type=str, default='.',
                    help='path where to store output (default: current dir)')
args = parser.parse_args()

# -------------------------
# main loop for specified ids
# -------------------------
for uniprot_id in args.ids:
    url = 'http://www.phosphonet.ca/?search=' + uniprot_id
    response = requests.get(url)
    
    # only process if page is OK
    if response.status_code == 200:
        soup = BeautifulSoup(response.text, "html.parser")
        
        phospho_sites = [ x.get_text() for x in soup.findAll("td", {"class": "pSiteNameCol"})]
        
        site_base_url = 'http://www.phosphonet.ca/kinasepredictor.aspx?uni=' + uniprot_id + '&ps='
        site_tables = []
        
        for idx, site in enumerate(phospho_sites):
            site_url = site_base_url + site
            response = requests.get(site_url)
            soup = BeautifulSoup(response.text, "html.parser")
            site_tables.append(soup)
            print(site_url)
            sleep(random.uniform(sleep_individual_low, sleep_individual_high))
            if idx%batch_size == 0:
                sleep(random.uniform(sleep_batch_low, sleep_batch_high))
        
        phos_sites_df = pd.DataFrame()
        
        for idx, site in enumerate(site_tables):
            strings = list(site.html.stripped_strings)
            start_index = strings.index("Kinase 1:")
            end_index = 50*7+49 # what is this
            string_array = np.array(strings[start_index:end_index])
            string_array = string_array.reshape((50,7), order='C')
            
            df = pd.DataFrame(string_array)
            df.drop([0,4,5], axis=1, inplace=True)
            df.reset_index(inplace=True)
            df.columns = ["kinase_no", "kinase_name", "uniprot_id", "kinexus_score", "kinexus_score_v2"]
            df.insert(0, "site", phospho_sites[idx][1:])
            df.insert(0, "aa", phospho_sites[idx][0])
            if(phos_sites_df.empty):
               phos_sites_df = df 
            else:
                phos_sites_df = pd.concat([phos_sites_df, df])
        
        # final reorganization of phospho site dataframe
        phos_sites_df['site'] = phos_sites_df['site'].astype(dtype='int32')
        phos_sites_df['kinexus_score'] = phos_sites_df['kinexus_score'].astype(dtype='int32')        
        phos_sites_df['kinexus_score_v2'] = phos_sites_df['kinexus_score_v2'].astype(dtype='int32')
        phos_sites_df.reset_index(drop=True)
        phos_sites_df['msite_shift'] = np.where(phos_sites_df['site'] < 170, 2, 1)
        phos_sites_df['msite'] = phos_sites_df['site'] - phos_sites_df['msite_shift']
        
        # TODO optionally include user supplied score filtering 
        # phos_sites_df.query('kinexus_score>700')
        
        # save to csv
        phos_sites_df.to_csv(args.outdir + "/" + uniprot_id + "_phos_kinexus.csv", 
                             sep='\t', 
                             index=False)
