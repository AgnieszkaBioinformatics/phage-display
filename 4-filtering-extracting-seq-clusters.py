import glob
from Bio import SeqIO
import re
import os

# filtering and extracting sequences from clusters (that have at least 10 sequences) for the assembly
# firstly we prepare txt file with reads per each cluster for seqtk to extract reads from original fastq files.

# we have a folder for each barcode where there are stored clusters | E.g. barcode 02 for 400 nt seq length


cwd02400=f"{cwd}/usearch-02-400-filtered/*"

count02400_lista=[]
for file in glob.glob(cwd02400):
    with open(file, 'r') as f:
        count02400=0
        for record in SeqIO.parse(f, 'fasta'):
            count02400+=1
        # it has to be at least 10 sequences per cluster
        if count02400 > 10:
            count02400_lista.append(file)



# we need to delete _poczatek/_koniec name description for seqtk to work
pattern='_poczatek|_koniec'

for file in count02400_lista:
    with open(file, 'r') as original, open(f'{file}_corrected', 'w') as corrected:
        for record in SeqIO.parse(original, 'fasta'):
            desc=record.description
            if desc.endswith('_poczatek'):
                new_desc=re.split(pattern, desc)[0]
                record.description=new_desc
                record.id=record.description
            elif desc.endswith('_koniec'):
                new_desc=re.split(pattern, desc)[0]
                record.description=new_desc
                record.id=record.description
            SeqIO.write(record, corrected, 'fasta')                


# saving sequences for each cluster
corrected_files=glob.glob(f"{cwd}/usearch-02-400-filtered/*_corrected")
names=[]
for klaster in corrected_files:
    with open(klaster, 'r') as k:
        klastry_names=[]
        for record in SeqIO.parse(k, 'fasta'):
            desc=record.description
            klastry_names.append(desc)

        with open(f'{cwd}/02-400-do-assembly/{os.path.basename(klaster)}', 'w') as kn:
            kn.write(str(klastry_names)) 
        names.append(klastry_names)