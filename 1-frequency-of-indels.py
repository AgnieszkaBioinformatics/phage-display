import PhageDisplayHelpers
import polars
import matplotlib.pyplot as plt

# loading sam file to the dataframe
map=polars.read_csv('mapowanie.sam', has_header=False, separator="\t", quote_char=None, skip_rows=5, new_columns=["qname", "flag", "rname", 
                                                                                                                  "pos", "mapq", "cigar", 
                                                                                                                  "rnext", "pnext", "tlen", 
                                                                                                                  "seq", "qscore", "mapqscore", "identity"])
map=map.select(polars.col("qname", "flag", "cigar", "seq", "qscore", "mapqscore", "identity"))

# adding inforamtion about inserions and deletions to the dataframe
dfs=[]
deletions=[]
insertions=[]

for match in map["cigar"]:
    ps,pe = PhageDisplayHelpers.decode_cigar(match)
    deletions.append(ps)
    insertions.append(pe)
    
ps=polars.Series("deletions",deletions)
pe=polars.Series("insertions", insertions)
map=map.hstack([ps, pe])

# adding information about length of the read
slen_lista=[]

for row in map.iter_rows():
    slen=len(row[3])
    slen_lista.append(slen)

l=polars.Series('len-seq', slen_lista)
map=map.hstack([l])

# calculating percentage of deletions and insertions
proc_ins_lista=[]
proc_dels_lista=[]

for row in map.iter_rows(named=True):
    proc_ins=(int(row['insertions']))*100/(int(row['len-seq']))
    proc_ins_lista.append(proc_ins)
    proc_dels=(int(row['deletions']))*100/(int(row['len-seq']))
    proc_dels_lista.append(proc_dels)

ins=polars.Series('proc-ins', proc_ins_lista)
dels=polars.Series('proc-dels', proc_dels_lista)
map=map.hstack([ins, dels])

# plotting results
plt.hist(map['proc-ins'], range =(0,1))
plt.hist(map['proc-dels'],  range =(0,1))