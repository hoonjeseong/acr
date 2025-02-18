import os,sys
import collections
import pathlib
import shutil
import kmeans1d
import pandas as pd
import numpy as np
from glob import glob
from subprocess import call
from sklearn.cluster import KMeans
from Bio.SeqIO.FastaIO import SimpleFastaParser

# get markers
BAC120_MARKERS = set(["PF00380.14", "PF00410.14", "PF00466.15", "PF01025.14", "PF02576.12", "PF03726.9","TIGR00006", "TIGR00019", "TIGR00020", "TIGR00029", "TIGR00043", "TIGR00054", "TIGR00059", "TIGR00061", "TIGR00064", "TIGR00065", "TIGR00082", "TIGR00083", "TIGR00084", "TIGR00086", "TIGR00088", "TIGR00090", "TIGR00092", "TIGR00095", "TIGR00115", "TIGR00116", "TIGR00138", "TIGR00158", "TIGR00166", "TIGR00168", "TIGR00186", "TIGR00194", "TIGR00250", "TIGR00337", "TIGR00344", "TIGR00362", "TIGR00382", "TIGR00392", "TIGR00396", "TIGR00398", "TIGR00414", "TIGR00416", "TIGR00420", "TIGR00431", "TIGR00435", "TIGR00436", "TIGR00442", "TIGR00445", "TIGR00456", "TIGR00459", "TIGR00460", "TIGR00468", "TIGR00472", "TIGR00487", "TIGR00496", "TIGR00539", "TIGR00580", "TIGR00593", "TIGR00615", "TIGR00631", "TIGR00634", "TIGR00635", "TIGR00643", "TIGR00663", "TIGR00717", "TIGR00755", "TIGR00810", "TIGR00922", "TIGR00928", "TIGR00959", "TIGR00963", "TIGR00964", "TIGR00967", "TIGR01009", "TIGR01011", "TIGR01017", "TIGR01021", "TIGR01029", "TIGR01032", "TIGR01039", "TIGR01044", "TIGR01059", "TIGR01063", "TIGR01066", "TIGR01071", "TIGR01079", "TIGR01082", "TIGR01087", "TIGR01128", "TIGR01146", "TIGR01164", "TIGR01169", "TIGR01171", "TIGR01302", "TIGR01391", "TIGR01393", "TIGR01394", "TIGR01510", "TIGR01632", "TIGR01951", "TIGR01953", "TIGR02012", "TIGR02013", "TIGR02027", "TIGR02075", "TIGR02191", "TIGR02273", "TIGR02350", "TIGR02386", "TIGR02397", "TIGR02432", "TIGR02729", "TIGR03263", "TIGR03594", "TIGR03625", "TIGR03632", "TIGR03654", "TIGR03723", "TIGR03725", "TIGR03953"])

AR122_MARKERS = set(["PF01868.11", "PF01282.14", "PF01655.13", "PF01092.14", "PF01000.21", "PF00368.13", "PF00827.12", "PF01269.12", "PF00466.15", "PF01015.13", "PF13685.1", "PF02978.14", "PF04919.7", "PF01984.15", "PF04104.9", "PF00410.14", "PF01798.13", "PF01864.12", "PF01990.12", "PF07541.7", "PF04019.7", "PF00900.15", "PF01090.14", "PF02006.11", "PF01157.13", "PF01191.14", "PF01866.12", "PF01198.14", "PF01496.14", "PF00687.16", "PF03874.11", "PF01194.12", "PF01200.13", "PF13656.1", "PF01280.15", "TIGR00468", "TIGR01060", "TIGR03627", "TIGR01020", "TIGR02258", "TIGR00293", "TIGR00389", "TIGR01012", "TIGR00490", "TIGR03677", "TIGR03636", "TIGR03722", "TIGR00458", "TIGR00291", "TIGR00670", "TIGR00064", "TIGR03629", "TIGR00021", "TIGR03672", "TIGR00111", "TIGR03684", "TIGR01077", "TIGR01213", "TIGR01080", "TIGR00501", "TIGR00729", "TIGR01038", "TIGR00270", "TIGR03628", "TIGR01028", "TIGR00521", "TIGR03671", "TIGR00240", "TIGR02390", "TIGR02338", "TIGR00037", "TIGR02076", "TIGR00335", "TIGR01025", "TIGR00471", "TIGR00336", "TIGR00522", "TIGR02153", "TIGR02651", "TIGR03674", "TIGR00323", "TIGR00134", "TIGR02236", "TIGR03683", "TIGR00491", "TIGR00658", "TIGR03680", "TIGR00392", "TIGR00422", "TIGR00279", "TIGR01052", "TIGR00442", "TIGR00308", "TIGR00398", "TIGR00456", "TIGR00549", "TIGR00408", "TIGR00432", "TIGR00264", "TIGR00982", "TIGR00324", "TIGR01952", "TIGR03626", "TIGR03670", "TIGR00337", "TIGR01046", "TIGR01018", "TIGR00936", "TIGR00463", "TIGR01309", "TIGR03653", "TIGR00042", "TIGR02389", "TIGR00307", "TIGR03673", "TIGR00373", "TIGR01008", "TIGR00283", "TIGR00425", "TIGR00405", "TIGR03665", "TIGR00448"])

# functions 
def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()

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

def Clustering(Fasta,Coverage,Size,process,jgi):
    np.random.seed(123456789)
    Bin=[];Fa_Size={}
    with open(Fasta,'r') as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            Bin.append(title.split()[0])
            Fa_Size[title.split()[0]]=len(seq)
    df=pd.read_table(Coverage, delimiter='\t',index_col=0)
    if jgi!='N':
        df=df[df.columns[2:][::2]]
    df=df.loc[Bin,:]
    if int(Size) > len(df.index):
        return {},0
    # start clustering
    if len(df.columns)==1: # 1-d kmeans
        df1 =df.iloc[:,0].to_list()
        df_name=df.index
        labels, centroids = kmeans1d.cluster(df1, int(Size))
    else: # kmeans
        df1 = df.values
        df_name=df.index
        # Number of clusters
        kmeans = KMeans(n_clusters=int(Size),n_init=30)
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

def Checking_info(Km,SCGs,Out_P,Fa_Size,Marker,Log,ID,comp,cont,tag,Path,Ex_contig,sub,Prefix,log_P,N,Ex):
    Maximum={}
    Log_F=open(log_P+'/acr.kmean.out','a') #open log
    ########################################### skip lower completeness 
    Check=[]
    for c in Km:
        for i in Km[c]:
            if i in set(Ex_contig):
                continue
            if i in SCGs:
                Check+=SCGs[i]
    BinStat=collections.Counter(Check)
    Completeness=len(list(filter(lambda x: BinStat[x]>=1,BinStat)))*100/float(len(Marker))
    if Completeness<comp:
        Log_F.write('# Break the loop because of lower completeness, Domain: '+tag+', ID: '+ID.split(Ex)[0]+', number of K:'+str(N)+', Completeness: '+str(Completeness)+'\n')
        return False,False,False
    ############################################
    for c in Km:
        Check=[];Size=0
        for i in Km[c]:
            Size+=Fa_Size[i]
            if i in SCGs:
                Check+=SCGs[i]
        BinStat=collections.Counter(Check)
        Completeness=len(list(filter(lambda x: BinStat[x]>=1,BinStat)))*100/float(len(Marker))
        Redundancy=sum([BinStat[k]-1 for k in BinStat.keys() if BinStat[k]>1])*100/float(len(Marker))
        Score=float(Completeness)-float(5*Redundancy)
        Log_F.write('Domain: '+tag+', ID: '+ID.split(Ex)[0]+', number of K:'+str(N)+', Completeness: '+str(Completeness)+', Redundancy: '+str(Redundancy)+', Size: '+str(Size)+', Score: '+str(Score)+', Group: '+ str(sub-1)+'\n')
        #if Score >= f_Score and not bool(set(Km[c]) & set(Ex_contig)): #revision 0731
        if Completeness >= comp and Redundancy < cont and not bool(set(Km[c]) & set(Ex_contig)): #revised 0731
            NewID=Out_P+'/'+Prefix+'.'+ID.split(Ex)[0]+'-'+str(sub-1)+'.'+tag+'.fa'
            Write_Fa(Path+'/'+ID,NewID,Km[c])
            Log[str(ID.split(Ex)[0]+'-'+str(sub-1)+'.'+tag)]={'Domain':tag,'Completeness':Completeness,'Redundancy':Redundancy,'Size':Size,'Score':Score}
            sub+=1
            Ex_contig+=Km[c]
    Log_F.close()
    return Log,Ex_contig,sub

def Euk_marker_score(Annot,Marker):
    for e in Marker:
        tmp_SCGs=[i for i in Annot if i in Marker[e]]
        Completness=float(len(tmp_SCGs)*100/float(len(Marker[e])))
        tmp_SCGs=collections.Counter(tmp_SCGs)
        Redundancy=sum([tmp_SCGs[i]-1 for i in tmp_SCGs if tmp_SCGs[i]>1])*100/float(len(Marker[e]))
        Score=Completness-5*Redundancy
        yield e,Score

def make_Marker(l,Marker,Check,SCGs):
    gID=l[3]
    if l[3]=='-': #for EukCC 
        gID=l[2].split('.')[0]
    if gID in Marker:
        Check.append(gID)
        if l[0].count('_')>1:   # <- in case : contig name includes '_' already
            SCGs.setdefault('_'.join(l[0].split('_')[:-1]),[]).append(gID)
        else:
            SCGs.setdefault(l[0].split('_')[0],[]).append(gID)
    return Check,SCGs

def select_E_marker(Marker):
    euk_m=str(pathlib.Path(__file__).parent.absolute())+'/../data/eukcc/sets/'
    EUK_MARKER=set()
    with open(euk_m+Marker+'.set','r') as M:
        for l in M:
            EUK_MARKER.add(l.rstrip('\n'))
    return EUK_MARKER

def check_Marker(Path,Ex,Out_P,i,Cov_P,BinStat,SCGs,Marker,comp,cont,tag,Prefix,process,jgi,log_P):
    Log={}
    Completeness=len(list(filter(lambda x: BinStat[x]>=1,BinStat)))*100/float(len(Marker))
    Redundancy=sum([BinStat[k]-1 for k in BinStat if BinStat[k]>1])*100/float(len(Marker))
    #if float(Completeness) - 5*float(Redundancy) >= f_Score: #revision comment 0731
    if float(Completeness) >= comp and float(Redundancy) < cont: #revised 0731
        shutil.copy(Path+'/'+i,Out_P+'/'+Prefix+'.'+i.split(Ex)[0]+'.'+tag+'.fa')
        Size=Cal_size(Path+'/'+i)
        Log[str(i.split(Ex)[0]+'.'+tag)]={'Domain':tag,'Completeness':Completeness,'Redundancy':Redundancy,'Size':Size,'Score':Completeness- 5*Redundancy}
    elif Completeness >= comp:
        ex_contig=[];sub=1
        for n in range(15): # k-test to 50
            Km,Fa_Size=Clustering(Path+'/'+i,Cov_P,n+1,process,jgi)
            Log,ex_contig,sub=Checking_info(Km,SCGs,Out_P,Fa_Size,Marker,Log,i,comp,cont,tag,Path,ex_contig,sub,Prefix,log_P,n+1,Ex)
            if Log==False:
                break
    return Log
