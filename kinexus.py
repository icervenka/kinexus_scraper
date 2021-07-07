#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 15:10:08 2019

@author: icervenka
"""

import argparse
import requests
import random
from bs4 import BeautifulSoup
from time import sleep
import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------
# constants
# ----------------------------------------------------------------------------
# urls and class names needed to extract phospho sites
phosphonet_base_url = 'http://www.phosphonet.ca/?search='
phosphonet_kinase_url = 'http://www.phosphonet.ca/kinasepredictor.aspx?uni='
phos_site_class = "pSiteNameCol"

# number of kinases reported per phospho site (50 as of now)
# TODO maybe move this to parser so the user can specify number of top kinases to return 
num_kinases_per_phos = 50
# each phospho site has currently 7 descriptor values
phos_column_values = 7

# ----------------------------------------------------------------------------
# functions
# ----------------------------------------------------------------------------

def get_phospho_sites(uniprot_id, base_url, tag_class):  
    url = base_url + uniprot_id
    response = requests.get(url)
    
    # only process if page is OK
    if response.status_code == 200:
        soup = BeautifulSoup(response.text, "html.parser")
        
        # extract all phospho sites based on css class
        phospho_sites = [ x.get_text() for x in soup.findAll("td", {"class": tag_class}) ]
        
    return phospho_sites
     

def get_kinases(uniprot_id, phospho_site, base_url):
    # Needs to reflect website changes
    kinase_begin_string = "Kinase 1:" 
    string_array = []
    url = base_url + uniprot_id + '&ps=' + phospho_site
    response = requests.get(url)
    
    # only process if page is OK
    if response.status_code == 200:
        print("querying kinases for: " + url)
        soup = BeautifulSoup(response.text, "html.parser")
        
        strings = list(soup.html.stripped_strings)
            
        # find where the kinase data begins based on string
        start_index = strings.index(kinase_begin_string)
        
        # end index is number of kinases timess columns plus the offset
        # of the first kinase
        end_index = num_kinases_per_phos * phos_column_values + start_index
        string_array = np.array(strings[start_index:end_index])
        string_array = string_array.reshape((num_kinases_per_phos,
                                             phos_column_values), 
                                            order='C')
    
    return string_array

def kinase_array_to_df(kinase_array):
    # convert reshaped array to df
    df = pd.DataFrame(kinase_array)
    # drop unused/duplicated columns
    df.drop([0,4,5], axis=1, inplace=True)
    df.reset_index(inplace=True)
    # rename columns to something useful
    # TODO maybe move to constants section or let user choose
    df.columns = ["kinase_no", "kinase_name", "uniprot_id", "kinexus_score", "kinexus_score_v2"]
    # insert the site position and amino acid at the beginning of data frame
    df.insert(0, "site", phospho_sites[idx][1:])
    df.insert(0, "aa", phospho_sites[idx][0])
    
    return(df)

def prettify_phos_df(phos_df):
    # final reorganization of phospho site dataframe
    # change to proper types
    phos_df['site'] = phos_df['site'].astype(dtype='int32')
    phos_df['kinexus_score'] = phos_df['kinexus_score'].astype(dtype='int32')        
    phos_df['kinexus_score_v2'] = phos_df['kinexus_score_v2'].astype(dtype='int32')
    phos_df.reset_index(drop=True)
    return phos_df


# ----------------------------------------------------------------------------
# commandline parser
# ----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Get protein phosphosites from phosphonet.ca')
parser.add_argument('ids', metavar='uniprot_id', type=str, nargs='+',
                    help='Human uniprot accession numbers of protiens to retrieve')
parser.add_argument('-o', '--outdir', type=str, default='.',
                    help='path where to store output, path must exist (default: current dir)')
parser.add_argument('--sil', type=int, default='2',
                    help='min sleep value between phosphosite queries (default: 2)')
parser.add_argument('--sih', type=int, default='5',
                    help='max sleep value between phosphosite queries (default: 5)')
parser.add_argument('--bs', type=int, default='30',
                    help='number of phosphosites per batch to query between sleep periods (default: 30)')
parser.add_argument('--sbl', type=int, default='30',
                    help='min sleep value between batch queries in seconds (default: 30)')
parser.add_argument('--sbh', type=int, default='40',
                    help='max sleep value between batch queries in seconds (default: 40)')
args = parser.parse_args()

# ----------------------------------------------------------------------------
# main loop for specified ids
# ----------------------------------------------------------------------------
for uniprot_id in args.ids:
    phospho_sites = get_phospho_sites(uniprot_id, phosphonet_base_url, phos_site_class)
    phospho_site_df = pd.DataFrame()
    
    for idx, site in enumerate(phospho_sites):
        kinases = get_kinases(uniprot_id, site, phosphonet_kinase_url)
        kinases = kinase_array_to_df(kinases)
        
        # if first run of the loop, initialize the final df or append otherwise
        if(phospho_site_df.empty):
            phospho_site_df = kinases
        else:
            phospho_site_df = pd.concat([phospho_site_df, kinases])
            
        # sleep periods are introduced to reduce the burden on the server
        # to avoid being kicked out
        # default sleep values make processing not feasible for large number
        # of proteins, but have been tested to work
        # you can speed things up, but might get disconnected
        sleep(random.uniform(args.sil, args.sih))
        if (idx+1)%args.bs == 0:
            print("wating between batches...")
            sleep(random.uniform(args.sbl, args.sbh))
    
    phospho_site_df = prettify_phos_df(phospho_site_df)
    
    # TODO optionally include user supplied score filtering 
    # phos_sites_df.query('kinexus_score>700')
    
    # save to csv
    phospho_site_df.to_csv(args.outdir + "/" + uniprot_id + "_phos_kinexus.csv", 
                           sep='\t', 
                           index=False)
    