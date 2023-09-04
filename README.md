# Gene Data Processing

This repository contains code for processing gene data obtained from GenBank. The code extracts gene IDs, dates, organisms, titles, journals, and PubMed IDs from GenBank records. It also downloads the DNA sequences associated with the gene IDs.

<div style="text-align:center;">
    ![](https://github.com/sajed-s/FASTA_Analysis/blob/main/istockphoto-1055935724-612x612.jpg)
</div>

This repository contains code for processing gene data obtained from GenBank. The code extracts gene IDs, dates, organisms, titles, journals, and PubMed IDs from GenBank records. It also downloads the DNA sequences associated with the gene IDs.
## Dependencies

The following dependencies are required to run the code:

- DendroPy
- pandas
- BioPython

Install the dependencies using the following command:

```bash
pip install dendropy pandas biopython
```
# Importing the libreries

```bash
import dendropy
from dendropy.interop import genbank
import pandas as pd
from Bio import SeqIO
```
# Importing the text file of gene IDs
In this part you should have a text file consist of your NCBI gene IDs

```bash
with open("E:\python\FMD\GI_A.txt", "r") as my_file:
    data = my_file.read()
    GIlist = data.replace('\n', ' ').split(" ")

data_part1 = GIlist[0:200]
data_part2 = GIlist[200:]
```
# Extracting the data
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
# Saving the data to CSV files

```bash
data_id_org_date = pd.DataFrame({'gene_id': gene_id, 'organism': organism, 'date': date})
data_pubmed_journal_title = pd.DataFrame({'journal': journal, 'pubmed_id': pubmed_id, 'title': title})
data_id_org_date.to_csv('data_id_org_date_a.csv')
data_pubmed_journal_title.to_csv('data_pubmed_journal_a.csv')
```

# Downloading the strains
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
# Detecting the mutated amino acids
To detect mutations in amino acid sequences, it is essential to align the sequences with each other. Several tools and libraries are available to perform sequence alignment, and the choice of the most suitable one depends on your specific goals. It is important to consider the lengths of your sequences and the number of samples when selecting an alignment method.

For optimal results, conducting the analysis on a Linux platform is recommended since many alignment libraries are not implemented on Windows. Now, let's delve into the code after the sequence alignment step:

```bash
from Bio import SeqIO
import pandas as pd
import numpy as np

# Read the aligned FASTA file
with open('ready_for_ditect.fas') as fasta_file:
    identifiers = []
    lengths = []
    sequence =[]
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):
        identifiers.append(seq_record.id)
        lengths.append(len(seq_record.seq))
        sequence.append(seq_record.seq)
```
## now lets make a pandas data frame to preform the analysis:
```bash
# Create a DataFrame from the aligned sequences
data = pd.DataFrame(sequence)

# Select the desired section of the aligned sequences
edited_data = data.iloc[:, 0:2336]

```
## Identify unique and mutated amino acid positions
```bash
non_mutated = []
mutated = []
def is_unique(df):
    for name in df:
        s = df[name]
        if len(np.unique(s)) != 1:
            mutated.append(name)
        else:
            non_mutated.append(name)

# Call the function to identify unique and mutated amino acids
is_unique(edited_data)
```

## Create DataFrames for non-mutated and mutated amino acids
```bash
result_nonmutated = pd.DataFrame({"non_mutated_AA": non_mutated})
result_mutated = pd.DataFrame({"mutated_AA": mutated})

# Convert the DataFrames to lists
nonmutated = result_nonmutated.iloc[:, 0].values.tolist()
mutated = result_mutated.iloc[:, 0].values.tolist()

mutated_aa = []
type_aa = []
quantity = []
```
## Obtain the types of amino acids for each mutated position
```bash
for AA in mutated:
    type_aa.append(edited_data[AA].unique())

# Pad the arrays to ensure they have the same length
max_length = max(len(arr) for arr in type_aa)
padded_type_aa = [np.pad(arr, (0, max_length - len(arr)), mode='constant', constant_values=np.nan) for arr in type_aa]

# Create a DataFrame to indicate the type of each mutated amino acid
df = {'mutated': mutated, 'type': type_aa}
pddf = pd.DataFrame(df)
```
## Determine the product being coded based on the position of mutated amino acids
```bash
conditions = [
    pddf['mutated'] < 201,
    (pddf['mutated'] >= 202) & (pddf['mutated'] <= 286),
    (pddf['mutated'] >= 287) & (pddf['mutated'] <= 504),
    (pddf['mutated'] >= 505) & (pddf['mutated'] <= 725),
    (pddf['mutated'] >= 726) & (pddf['mutated'] <= 935),
    (pddf['mutated'] >= 736) & (pddf['mutated'] <= 953),
    (pddf['mutated'] >= 954) & (pddf['mutated'] <= 1107),
    (pddf['mutated'] >= 1108) & (pddf['mutated'] <= 1425),
    (pddf['mutated'] >= 1426) & (pddf['mutated'] <= 1578),
    (pddf['mutated'] >= 1579) & (pddf['mutated'] <= 1601),
    (pddf['mutated'] >= 1602) & (pddf['mutated'] <= 1625),
    (pddf['mutated'] >= 1626) & (pddf['mutated'] <= 1649),
    (pddf['mutated'] >= 1650) & (pddf['mutated'] <= 1862)
]

choices = [
    'Leader protease', 'VP4', 'VP2', 'VP3', 'VP1', 'Protein 2A',
    'Protein 2B', 'Protein 2C', 'Protein 3A', 'Protein 3B-1',
    'Protein 3B-2', 'Protein 3b-3', 'Protein 3C'
]

pddf['product'] = np.select(conditions, choices, default='unknown')
```

## Save the results to a CSV file
```bash
pddf.to_csv('result_AA.csv')
```
By following this code, you can detect mutations in amino acid sequences that have been aligned. Remember to provide the aligned sequences in the appropriate FASTA format.







