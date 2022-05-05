import sys
import numpy as np
import pandas as pd
#bowtie = pd.read_csv("./bowtie_compare/sequence2.mapq", delimiter="\t",header=None)
local_count=pd.read_csv(sys.argv[1], delimiter=" ",header=None)
global_count=pd.read_csv(sys.argv[2], delimiter=" ",header=None)
	

with open(sys.argv[3]) as f:
    firstline = f.readline().rstrip()

kmer=31
chrom,source,type,start,end,score,strand,phase,attributes=firstline.split()
name=chrom+" "+start+" "+end
frames = [local_count[1], global_count[1]]
df=pd.concat(frames,axis=1)
df['ratio']=df.iloc[:, 0] /df.iloc[:,1]

df['position']=df.index
df['position'] = df['position'].apply(lambda x: x+int(start))
df['chromosome']=chrom
df['source']=source
df['type']=type
df['start']=df['position']
df['end']=df['position']
df['stand']=strand
df['phase']=phase
df['attributes']=df['ratio']
df['ratio_mean_base']=df['ratio'].rolling(kmer).mean()
df['ratio_median_base']=df['ratio'].rolling(kmer).median()
df[0:kmer-1]=df[0:kmer-1].assign(ratio_mean_base=df['ratio'][0:kmer-1].expanding().mean(),ratio_median_base=df['ratio'][0:kmer-1].expanding().median())
df['ratio_mean_kmer']=df['ratio_median_base'].rolling(kmer).mean()
df['ratio_median_kmer']=df['ratio_median_base'].rolling(kmer).median()
df['ratio_max_kmer']=df['ratio_median_base'].rolling(kmer).max()

df[0:kmer-1]=df[0:kmer-1].assign(ratio_mean_kmer=df['ratio_median_base'][0:kmer-1].expanding().mean(),ratio_median_kmer=df['ratio_median_base'][0:kmer-1].expanding().median(),ratio_max_kmer=df['ratio_median_base'][0:kmer-1].expanding().max())
#df['bt_per_base']=(df['bt'].rolling(kmer).median())
#df[0:kmer-1]=df[0:kmer-1].assign(bt_per_base=df['bt'][0:kmer-1].expanding().median())
#df['bt_median']=df['bt'].rolling(kmer).median()
import matplotlib.pyplot as plt
import pandas as pd
  
  
# Initialize the lists for X and Y
  
plt.figure(figsize=(20, 3))  # width:20, height:3
 
X = list(df['position'])
Y = list(df['ratio_max_kmer'])

# Plot the data using bar() method
plt.bar(X, Y, color='g',align='center',width=1,alpha=1)
#plt.bar(X, Z, color='y',align='center',width=1,alpha=1)

plt.title(name)
plt.xlabel("Position")
plt.ylabel("Ratio")
png=name+'_newscore.png'
plt.savefig(png)
df['start_bed']=df['start']-1
df['end_bed']=df['start_bed']+kmer-1
bedfile=pd.DataFrame()
bedfile['chromosome']=df['chromosome']
bedfile['start']=df['start_bed']
bedfile['end']=df['end_bed']
bedfile['score']=df['ratio_max_kmer']
bedname=name+'.bed'
bedfile.to_csv(bedname,index=False)
