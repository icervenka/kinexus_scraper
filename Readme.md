# Kinexus Kinase Scraper
Simple web scraper that downloads lists of kinases that phosphorylate specific
protein residues based on http://www.phosphonet.ca. Takes uniprot protein id as
an input (only human proteins are supported). Several sleep timers are available
for customization to decrease the load on the webpage or avoid disconnect.
This scraper can be used in conjunction with [pamgene](https://github.com/icervenka/pamgene_analysis) analysis script to download
the needed phosphorylation data for kinase prediction.

With default settings, the runtime is around 1 min for 10 phosphorylatable amino acids.

# Example usage
Downloads top predicted kinases for all residues of PPARGC1A protein

`
kinexus.py Q9UBK2
`

# Result
| aa | site | kinase_no | kinase_name          | uniprot_id | kinexus_score | kinexus_score_v2 |
|----|------|-----------|----------------------|------------|---------------|------------------|
| S  | 15   | 1         | CK2a1 (CSNK2A1)      | P68400     | 390           | 250              |
| S  | 15   | 2         | CK2a2 (CSNK2A2)      | P19784     | 382           | 245              |
| S  | 15   | 3         | BARK2 (ADRBK2, GRK3) | P35626     | 315           | 202              |
| S  | 15   | 4         | PFTAIRE2 (ALS2CR7)   | Q96Q40     | 309           | 199              |
| S  | 15   | 5         | BARK1 (ADRBK1, GRK2) | P25098     | 309           | 198              |
| S  | 15   | 6         | NEK8 (AI086865)      | Q86SG6     | 299           | 192              |


# Commandline arguments
Sleep time is chosen as random between minimum and maximum arguments

- `--sil` - min sleep time between individual queries
- `--sih` - max sleep time between individual queries
- `--bs` - batch size of queries
- `--sbl` - min sleep time between query batches
- `--sbh` - max sleep time between query batches
