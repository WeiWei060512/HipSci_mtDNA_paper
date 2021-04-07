#!/usr/bin/env python
import sys
import glob
import numpy as np

def RemovedRegion_mt(all_pos):
        RemovedRegion=[66,67,68,69,70,71,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,513,514,515,516,517,518,519,520,521,522,523,524,525,3106,3107,12418,12419,12420,12421,12422,12423,12424,12425,16180,16181,16182,16183,16184,16185,16186,16187,16188,16189,16190,16191,16192,16193,16194]
        keep_pos = all_pos[~all_pos['POS'].isin(RemovedRegion)]
        return keep_nonDloop

def remove_indels(Vari):
        mask = (Vari['REF'].str.len() == 1) & (Vari['ALT'].str.len() == 1)
        SNPonly = Vari.loc[mask]
        return SNPonly

def filter_DP_HF(df):

        df2 = pd.DataFrame(df.ix[:, -1].str.split(':').tolist(),columns = ['GT','GQ','SDP','DP','RD','AD','FREQ','PVAL','RBQ','ABQ','RDF','RDR','ADF','ADR'])
        df3 = pd.concat([df, df2], axis=1)
        df4 = df3[df3['DP'].astype(int)>200]
        df4['FREQ'] = df4['FREQ'].str.rstrip('%').astype('float') / 100.0
        print (df4['FREQ'])
        df5 = df4[df4['FREQ'].astype(float)>0.02]
        df6 = df5[(df5['ADF'].astype(int)>1) & (df5['ADR'].astype(int)>1)]
        return df5

vcf_files = sys.argv[1]
for sample in vcf_files:
        with open(sample, 'rU') as input1:
                cellName = sample.replace('.vcf','')
                df = pd.read_csv(input1,skiprows=23, sep="\t", engine='python')
                print (sample)
                df1 = filter_DP_HF(df)
                df2 = RemovedRegion_mt(df1)
                df3 = remove_indels(df2)
                df3['cellID']=cellName
                df3.to_csv(cellName + '.filtered.tsv', sep = '\t',header=True, index=False)
