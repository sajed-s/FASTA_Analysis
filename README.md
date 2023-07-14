# Gene Data Processing

This repository contains code for processing gene data obtained from GenBank. The code extracts gene IDs, dates, organisms, titles, journals, and PubMed IDs from GenBank records. It also downloads the DNA sequences associated with the gene IDs using computer vision techniques.

## Dependencies

The following dependencies are required to run the code:

- DendroPy
- pandas
- BioPython

Install the dependencies using the following command:

```bash
pip install dendropy pandas biopython
```
#Importing the libreries

```bash
import dendropy
from dendropy.interop import genbank
import pandas as pd
from Bio import SeqIO
```
#Importing the text file of gene IDs
In this part you should have a text file consist of your NCBI gene IDs

```bash
with open("E:\python\FMD\GI_A.txt", "r") as my_file:
    data = my_file.read()
    GIlist = data.replace('\n', ' ').split(" ")

data_part1 = GIlist[0:200]
data_part2 = GIlist[200:]
```
#Extracting the data
In this part I would like to extract some information for each of my samples including:
- Gene_id
- Colloction date
- organism name
- Research_title
- Journal_name
- Pubmed_Id
* I had to segmented the data into two sections because the number of the samples were to much

```bash
gene_id = []
date = []
organism = []
title = []
journal = []
pubmed_id = []

gb_dna = genbank.GenBankDna(data_part1)
for gb_rec in gb_dna:
    gene_id.append(gb_rec.gi)
    date.append(gb_rec.create_date)
    organism.append(gb_rec.organism)
    for ref in gb_rec.references:
        pubmed_id.append(ref.pubmed_id)
        journal.append(ref.journal)
        title.append(ref.title)

gb_dna = genbank.GenBankDna(data_part2)
for gb_rec in gb_dna:
    gene_id.append(gb_rec.gi)
    date.append(gb_rec.create_date)
    organism.append(gb_rec.organism)
    for ref in gb_rec.references:
        pubmed_id.append(ref.pubmed_id)
        journal.append(ref.journal)
        title.append(ref.title)
```
#Saving the data to CSV files

```bash
data_id_org_date = pd.DataFrame({'gene_id': gene_id, 'organism': organism, 'date': date})
data_pubmed_journal_title = pd.DataFrame({'journal': journal, 'pubmed_id': pubmed_id, 'title': title})
data_id_org_date.to_csv('data_id_org_date_a.csv')
data_pubmed_journal_title.to_csv('data_pubmed_journal_a.csv')
```

#Downloading the strains
Now lets download the FASTA sequences of each samples
```bash
def get_fasta_a(ids):
    for i in ids:
        yield 'strain_a', genbank.GenBankDna(ids=i)

sampled1 = open('sample_a1.fasta', 'w')
for species, recs in get_fasta_a([data_part1]):
    tn = dendropy.TaxonNamespace()
    char_mat = recs.generate_char_matrix(taxon_namespace=tn, gb_to_taxon_fn=lambda gb: tn.require_taxon(label='%s_%s' % (species, gb.accession)))
    char_mat.write_to_stream(sampled1, 'fasta')
sampled1.close()

sampled2 = open('sample_a2.fasta', 'w')
for species, recs in get_fasta_a([data_part2]):
    tn = dendropy.TaxonNamespace()
    char_mat = recs.generate_char_matrix(taxon_namespace=tn, gb_to_taxon_fn=lambda gb: tn.require_taxon(label='%s_%s' % (species, gb.accession)))
    char_mat.write_to_stream(sampled2, 'fasta')
sampled2.close()
```
#Merging the fasta files
As mentioned before, I had to segment the data into two different parts. Now let's merge them all.

```bash
input_files = ["sample_a1.fasta", "sample_a2.fasta"]
output_file = "total_a.fasta"

with open(output_file, "w") as outfile:
    for input_file in input_files:
        for record in SeqIO.parse(input_file, "fasta"):
            SeqIO.write(record, outfile, "fasta")

print("Merged fasta files written to", output_file)

```


