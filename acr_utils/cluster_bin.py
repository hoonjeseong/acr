import os,sys
import collections
import shutil
import kmeans1d
import pandas as pd
import numpy as np
from subprocess import call
from sklearn.cluster import KMeans
from Bio.SeqIO.FastaIO import SimpleFastaParser

BAC120_MARKERS = set(["PF00380.14", "PF00410.14", "PF00466.15", "PF01025.14", "PF02576.12", "PF03726.9","TIGR00006", "TIGR00019", "TIGR00020", "TIGR00029", "TIGR00043", "TIGR00054", "TIGR00059", "TIGR00061", "TIGR00064", "TIGR00065", "TIGR00082", "TIGR00083", "TIGR00084", "TIGR00086", "TIGR00088", "TIGR00090", "TIGR00092", "TIGR00095", "TIGR00115", "TIGR00116", "TIGR00138", "TIGR00158", "TIGR00166", "TIGR00168", "TIGR00186", "TIGR00194", "TIGR00250", "TIGR00337", "TIGR00344", "TIGR00362", "TIGR00382", "TIGR00392", "TIGR00396", "TIGR00398", "TIGR00414", "TIGR00416", "TIGR00420", "TIGR00431", "TIGR00435", "TIGR00436", "TIGR00442", "TIGR00445", "TIGR00456", "TIGR00459", "TIGR00460", "TIGR00468", "TIGR00472", "TIGR00487", "TIGR00496", "TIGR00539", "TIGR00580", "TIGR00593", "TIGR00615", "TIGR00631", "TIGR00634", "TIGR00635", "TIGR00643", "TIGR00663", "TIGR00717", "TIGR00755", "TIGR00810", "TIGR00922", "TIGR00928", "TIGR00959", "TIGR00963", "TIGR00964", "TIGR00967", "TIGR01009", "TIGR01011", "TIGR01017", "TIGR01021", "TIGR01029", "TIGR01032", "TIGR01039", "TIGR01044", "TIGR01059", "TIGR01063", "TIGR01066", "TIGR01071", "TIGR01079", "TIGR01082", "TIGR01087", "TIGR01128", "TIGR01146", "TIGR01164", "TIGR01169", "TIGR01171", "TIGR01302", "TIGR01391", "TIGR01393", "TIGR01394", "TIGR01510", "TIGR01632", "TIGR01951", "TIGR01953", "TIGR02012", "TIGR02013", "TIGR02027", "TIGR02075", "TIGR02191", "TIGR02273", "TIGR02350", "TIGR02386", "TIGR02397", "TIGR02432", "TIGR02729", "TIGR03263", "TIGR03594", "TIGR03625", "TIGR03632", "TIGR03654", "TIGR03723", "TIGR03725", "TIGR03953"])

AR122_MARKERS = set(["PF01868.11", "PF01282.14", "PF01655.13", "PF01092.14", "PF01000.21", "PF00368.13", "PF00827.12", "PF01269.12", "PF00466.15", "PF01015.13", "PF13685.1", "PF02978.14", "PF04919.7", "PF01984.15", "PF04104.9", "PF00410.14", "PF01798.13", "PF01864.12", "PF01990.12", "PF07541.7", "PF04019.7", "PF00900.15", "PF01090.14", "PF02006.11", "PF01157.13", "PF01191.14", "PF01866.12", "PF01198.14", "PF01496.14", "PF00687.16", "PF03874.11", "PF01194.12", "PF01200.13", "PF13656.1", "PF01280.15", "TIGR00468", "TIGR01060", "TIGR03627", "TIGR01020", "TIGR02258", "TIGR00293", "TIGR00389", "TIGR01012", "TIGR00490", "TIGR03677", "TIGR03636", "TIGR03722", "TIGR00458", "TIGR00291", "TIGR00670", "TIGR00064", "TIGR03629", "TIGR00021", "TIGR03672", "TIGR00111", "TIGR03684", "TIGR01077", "TIGR01213", "TIGR01080", "TIGR00501", "TIGR00729", "TIGR01038", "TIGR00270", "TIGR03628", "TIGR01028", "TIGR00521", "TIGR03671", "TIGR00240", "TIGR02390", "TIGR02338", "TIGR00037", "TIGR02076", "TIGR00335", "TIGR01025", "TIGR00471", "TIGR00336", "TIGR00522", "TIGR02153", "TIGR02651", "TIGR03674", "TIGR00323", "TIGR00134", "TIGR02236", "TIGR03683", "TIGR00491", "TIGR00658", "TIGR03680", "TIGR00392", "TIGR00422", "TIGR00279", "TIGR01052", "TIGR00442", "TIGR00308", "TIGR00398", "TIGR00456", "TIGR00549", "TIGR00408", "TIGR00432", "TIGR00264", "TIGR00982", "TIGR00324", "TIGR01952", "TIGR03626", "TIGR03670", "TIGR00337", "TIGR01046", "TIGR01018", "TIGR00936", "TIGR00463", "TIGR01309", "TIGR03653", "TIGR00042", "TIGR02389", "TIGR00307", "TIGR03673", "TIGR00373", "TIGR01008", "TIGR00283", "TIGR00425", "TIGR00405", "TIGR03665", "TIGR00448"])

def Write_Fa(Fasta,NewF,Contigs):
    with open(Fasta,'r') as in_handle, open(NewF,'w') as out_handle:
        for title, seq in SimpleFastaParser(in_handle):
            if title.split()[0] in Contigs:
                out_handle.write('>'+title+'\n'+seq+'\n')

def Cal_size(Fasta):
    Size=0
    with open(Fasta,'r') as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            Size+=len(seq)
    return Size

def Clustering(Fasta,Coverage,Size,cpu):
    Bin=[];Fa_Size={}
    with open(Fasta,'r') as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            Bin.append(title.split()[0])
            Fa_Size[title.split()[0]]=len(seq)
    df=pd.read_table(Coverage, delimiter='\t',index_col=0)
    df=df.loc[Bin,:]
    if int(Size) > len(df.index):
        return {},0
    df1 = df.values
    df_name=df.index
    # Number of clusters
    kmeans = KMeans(n_clusters=int(Size),n_init=30,n_jobs=cpu)
    # Fitting the input data
    kmeans = kmeans.fit(df1)
    # Getting the cluster labels
    labels = kmeans.predict(df1)
    # Centroid values
    centroids = kmeans.cluster_centers_
    Sep_fa={}
    for i,j in enumerate(df_name):
        Sep_fa.setdefault(str(labels[i]),[]).append(j)
    return Sep_fa,Fa_Size

def Checking_info(Km,SCGs,Out_P,Fa_Size,Marker,Log,ID,tag,Path,Ex_contig,sub,Prefix):
    Maximum={}
    for c in Km:
        Check=[];Size=0
        for i in Km[c]:
            Size+=Fa_Size[i]
            if i in SCGs:
                Check+=SCGs[i]
        BinStat=collections.Counter(Check)
        Completeness=len(list(filter(lambda x: BinStat[x]>=1,BinStat.keys())))*100/float(len(Marker))
        Redundancy=sum([BinStat[k]-1 for k in BinStat.keys() if BinStat[k]>1])*100/float(len(Marker))
        Score=float(Completeness)-float(5*Redundancy)
        print ('Completeness:',Completeness,'Redundancy:',Redundancy,'Size:',Size,'Score:',Score,'ID:',ID.split('.fa')[0], sub)
        if Score >= 50 and not bool(set(Km[c]) & set(Ex_contig)):
            NewID=Out_P+'/'+Prefix+'.'+ID.split('.fa')[0]+'-'+str(sub)+'.'+tag
            Write_Fa(Path+'/'+ID,NewID,Km[c])
            Log[str(ID.split('.fa')[0]+'-'+str(sub)+'.'+tag)]={'Completeness':Completeness,'Redundancy':Redundancy,'Size':Size,'Score':Score}
            sub+=1
            Ex_contig+=Km[c]
    return Log,Ex_contig,sub 

def make_Marker(l,Marker,Check,SCGs):
    if l[3] in Marker:
        Check.append(l[3])
        if l[0].count('_')>1:   # <- in case : contig name includes '_' already
            SCGs.setdefault('_'.join(l[0].split('_')[:-1]),[]).append(l[3])
        else:
            SCGs.setdefault(l[0].split('_')[0],[]).append(l[3])
    return Check,SCGs

def check_Marker(Path,Ex,Out_P,i,Cov_P,BinStat,SCGs,Marker,tag,Prefix,cpu):
    Log={}
    Completeness=len(list(filter(lambda x: BinStat[x]>=1,BinStat.keys())))*100/float(len(Marker))
    Redundancy=sum([BinStat[k]-1 for k in BinStat.keys() if BinStat[k]>1])*100/float(len(Marker))
    if float(Completeness) - 5*float(Redundancy) >= 50:
        shutil.copy(Path+'/'+i,Out_P+'/'+Prefix+'.'+i.split('.fa')[0]+'.'+tag)
        Size=Cal_size(Path+'/'+i)
        Log[str(i.split(Ex)[0]+'.'+tag)]={'Completeness':Completeness,'Redundancy':Redundancy,'Size':Size,'Score':Completeness- 5*Redundancy}
    elif Completeness >= 50:
        ex_contig=[];sub=0
        for n in range(15):
            Km,Fa_Size=Clustering(Path+'/'+i,Cov_P,n+1,cpu)
            Log,ex_contig,sub=Checking_info(Km,SCGs,Out_P,Fa_Size,Marker,Log,i,tag,Path,ex_contig,sub,Prefix)
            print ('times:',n)
    return Log
