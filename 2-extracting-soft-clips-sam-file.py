import PhageDisplayHelpers
import polars
import re
'''
extracting soft clips from sam file (reads mapped to the phage seq assembly) that were not mapped to the phage assembly 
to get sequences of potential antibody chains
'''

# loading sam file to the dataframe
map=polars.read_csv('all-mapped-reads.sam', has_header=False, separator="\t", quote_char=None, skip_rows=3, new_columns=["qname", "flag", "rname",
                                                                                                                          "pos", "mapq", "cigar", 
                                                                                                                          "rnext", "pnext", "tlen", 
                                                                                                                          "seq", "qscore", "mapqscore", 
                                                                                                                          "identity"])
map=map.select(polars.col("qname", "flag", "cigar", "seq", "qscore", "mapqscore", "identity"))

# adding positions of soft clips to the dataframe
q_s=[]
q_e=[]


for match in map["cigar"]:
    ps,pe = PhageDisplayHelpers.decode_cigar(match)
    q_s.append(ps)
    q_e.append(pe)
    
ps=polars.Series("start",q_s)
pe=polars.Series("end", q_e)
map=map.hstack([ps, pe])

# adding information about read length
slen_lista=[]

for row in map.iter_rows():
    slen=len(row[3])
    slen_lista.append(slen)

l=polars.Series('len-seq', slen_lista)
map=map.hstack([l])


# filtering out the unmapped reads
new_df=map.filter(map['seq']!="*")

# finding and extracting soft clips information in the cigar string
pattern=re.compile('\d+S')

sc_list=[]

for row in new_df.iter_rows(named=True):   
    match=pattern.findall(row['cigar'])
    sc_list.append(match)

sc_seria=polars.Series('soft-clips', sc_list)
new_df=new_df.hstack([sc_seria])

df = new_df.with_columns(new_df['soft-clips'].apply(lambda lst: len(lst) > 0).alias('non_empty'))

# Filtering rows where 'soft-clips' column is not empty
df_sc = df.filter(df['non_empty'])


# getting coordinates of soft clips positions
sc_coord_d_lista=[]

for row in df_sc.iter_rows(named=True):
    sc_coord=[]
    for i in row['soft-clips']:
        it=i.split('S')[0]
        sc_coord.append(int(it))
    sc_coord_d_lista.append(sc_coord)

int_seria=polars.Series('sc-coord', sc_coord_d_lista)
df_sc=df_sc.hstack([int_seria])


# extracting soft clips sequences
uciete_sekw_d_lista = []

for row in df_sc.iter_rows(named=True):
    row_seq = []
    # checking if the sequence has soft clips on both ends
    if len(row['sc-coord']) > 1:
        seq = row['seq'][:row['sc-coord'][0]]  
        seqd = row['seq'][-row['sc-coord'][1]:]
        # flag "16" indicates that sequences are reversed
        if row['flag']==16:
            seq=PhageDisplayHelpers.revcomp(seq)
            seqd=PhageDisplayHelpers.revcomp(seqd)
        row_seq.append(seq)
        row_seq.append(seqd)
    # if the sequence has soft clip only on one side, we need to check on which
    elif len(row['sc-coord']) == 1:
        # checking if it is in the beginning of the sequence
        if row['cigar'].startswith(row['soft-clips'][0]):
            seq = 'P'+row['seq'][:row['sc-coord'][0]]    # with 'P' we indicate that it is in the begining
            if row['flag']==16:
                seq=PhageDisplayHelpers.revcomp(seq)
            row_seq.append(seq)
        else:
            seq = 'K'+row['seq'][-row['sc-coord'][0]:]  # with 'K' we indicate that it is in the ending of the seq
            if row['flag']==16:
                seq=PhageDisplayHelpers.revcomp(seq)
                row_seq.append(seq)
    uciete_sekw_d_lista.append(row_seq)

    # saving soft clips sequences in the dataframe
    sc_sekwencjes=polars.Series('sekwencje', uciete_sekw_d_lista)
    df_sc=df_sc.hstack([sc_sekwencjes])

# first we do this for single chain of the antibody
# we store the potential antibody sequences with distinction if the soft clips were located in the begining or in the ending of the sequence
sc1_lista400=[]
sc2_lista400=[]

for row in df_sc.iter_rows(named=True):
    if len(row['sekwencje'])==1: 
        sc1400=row['sekwencje'][0]
        name1=row['qname']
        if sc1400.startswith('P') and len(sc1400) in range(400, 550):
            sc1_lista400.append(('>', name1+'_poczatek', '\n', sc1400[1:], '\n'))
        elif sc1400.startswith('K') and len(sc1400) in range(400, 550):
            sc2_lista400.append(('>', name1+'_koniec', '\n', sc1400[1:], '\n'))
    
    if len(row['sekwencje']) > 1:
        sc11400=row['sekwencje'][0]
        sc2400=row['sekwencje'][1]
        name2=row['qname']
        if len(sc11400) in range(400,550):
            sc1_lista400.append(('>', name2+'_poczatek', '\n', sc11400, '\n'))
        if len(sc2400) in range(400,550):
            sc2_lista400.append(('>', name2+'_koniec', '\n', sc2400, '\n'))

# now for double chain of the antibody
sc1_lista700=[]
sc2_lista700=[]

for row in df_sc.iter_rows(named=True):
    if len(row['sekwencje'])==1: 
        sc1700=row['sekwencje'][0]
        name1=row['qname']
        if sc1700.startswith('P') and len(sc1700) in range(700, 850):
            sc1_lista700.append(('>', name1+'_poczatek', '\n', sc1700, '\n'))
        elif sc1700.startswith('K') and len(sc1700) in range(700, 850):
            sc2_lista700.append(('>', name1+'_koniec', '\n', sc1700, '\n'))
    
    if len(row['sekwencje']) > 1:
        sc11700=row['sekwencje'][0]
        sc2700=row['sekwencje'][1]
        name2=row['qname']
        if len(sc11700) in range(700,850):
            sc1_lista700.append(('>', name2+'_poczatek', '\n', sc11700, '\n'))
        if len(sc2700) in range(700,850):
            sc2_lista700.append(('>', name2+'_koniec', '\n', sc2700, '\n'))

all_lista400=sc1_lista400+sc2_lista400
all_lista700=sc1_lista700+sc2_lista700

# the sequences were saved under the name:

with open('all_sc_02_400.fasta', 'w') as f:
    for i in all_lista400:
        for t in i:
            f.write(str(t))


with open('all_sc_02_700.fasta', 'w') as f:
    for i in all_lista700:
        for t in i:
            f.write(str(t))