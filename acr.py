import sys,os
import shutil
import pathlib
import datetime
import optparse
import logging
import subprocess
import collections
import pandas as pd
import acr_utils.cluster_bin as cb

def run_cmd(cmd,log_dir,shell=False,program=False):
    # make log
    log_P=os.path.abspath(log_dir)
    Out = os.path.join(log_P+'/acr.log.out')
    Err = os.path.join(log_P+'/acr.log.err')

    with open(Out, 'a') as L:
        now=str(datetime.datetime.now().date()) + '_' + \
            str(datetime.datetime.now().time()).replace(':', '.')
        if program:
            now='#### '+program+' running : '+now+' ####'
        if shell:
            L.write('time: '+now+'\n'+str(cmd)+ '\n')
        else:
            L.write('time: '+now+'\n'+' '.join(cmd)+'\n')
    # run 
    if shell:
        with open(Out,'a') as sto, open(Err,'a') as ste:
            subprocess.call(cmd,stdout=sto,stderr=ste,shell=True)
    else:
        with open(Out,'a') as sto, open(Err,'a') as ste:
            proc=subprocess.Popen(cmd,stdout=sto, stderr=ste)
            proc.communicate()

def find_p(name):
    # check the program path
    loc = shutil.which(name)
    if loc != None:
        try:
            o=subprocess.check_output([loc,'-h'],stderr=subprocess.STDOUT,shell=False,universal_newlines=True)
            return loc
        except:
            logging.error("Please set the {0} path !\n".format(loc)) 
            sys.exit()

def score_sort(blast_out,blast_filter):
    temp=set()
    evalue={}
    with open(blast_out,'r') as B, open(blast_filter,'w') as N:
        while True:
            line=B.readline()
            if not line: break
            if not line.startswith('#'):
                spl=line.rstrip('\n').split()
                temp.add(spl[0])
                evalue.setdefault(spl[0],{})
                evalue[spl[0]].setdefault(float(spl[5]),spl[3])
        for i in evalue:
            evalue[i]=evalue[i][max(set(evalue[i].keys()))]
        B.seek(0)
        while True:
            line=B.readline()
            if not line: break
            spl=line.rstrip('\n').split()
            if spl[0] in evalue and spl[3] == evalue[spl[0]]:
                N.write(line)

def run_HMM_to_Marker(prodigal,hmmsearch,bin_F,bin_P,name,Pfam,Tigrfam,log_P):
    cmd=[prodigal,'-i',bin_F,'-a',bin_P+'/'+name+'.faa','-d',bin_P+'/'+name+'.fna','-f','gff','-p','meta','-o',bin_P+'/'+name+'.gff']
    run_cmd(cmd,program='prodigal',log_dir=log_P)
    statinfo = os.stat(bin_P+'/'+name+'.faa')
    if float(statinfo.st_size)==0:
        return 1
    cmd=[hmmsearch,'--tblout',bin_P+'/gene.P.hit','--cut_ga',Pfam,bin_P+'/'+name+'.faa','>',bin_P+'/gene.Pfam']
    run_cmd(' '.join(cmd),program='hmmsearch',log_dir=log_P,shell=True)
    score_sort(bin_P+'/gene.P.hit',bin_P+'/gene.f.P.hit')

    cmd=[hmmsearch,'--tblout',bin_P+'/gene.T.hit','--cut_ga',Tigrfam,bin_P+'/'+name+'.faa','>',bin_P+'/gene.Tfam']
    run_cmd(' '.join(cmd),program='hmmsearch',log_dir=log_P,shell=True)
    score_sort(bin_P+'/gene.T.hit',bin_P+'/gene.f.T.hit')

    cmd=['cat',bin_P+'/gene.f.P.hit',bin_P+'/gene.f.T.hit']#,'>',bin_P+'/gene.Total.hit'] 
    with open(bin_P+'/gene.Total.hit','w') as sto:
        proc=subprocess.Popen(cmd,stdout=sto)#,stderr=subprocess.PIPE)
        proc.communicate()
    return 0

def check_cluster(bin_P,ex,gF_P,oF_P,i,cov,preF,cpu):
    BSCGs={};ASCGs={}
    BCheck=[];ACheck=[]
    with open(bin_P+'/gene.Total.hit','r') as B:
        for l in B:
            if not l.startswith('#'):
                l=l.split()
                BCheck,BSCGs=cb.make_Marker(l,cb.BAC120_MARKERS,BCheck,BSCGs)
                ACheck,ASCGs=cb.make_Marker(l,cb.AR122_MARKERS,ACheck,ASCGs)
    BBinStat=collections.Counter(BCheck)
    ABinStat=collections.Counter(ACheck)
    Sb=cb.check_Marker(gF_P,ex,oF_P,i,cov,BBinStat,BSCGs,cb.BAC120_MARKERS,'Bac',preF,cpu)
    Sa=cb.check_Marker(gF_P,ex,oF_P,i,cov,ABinStat,ASCGs,cb.AR122_MARKERS,'Arc',preF,cpu) 
    return Sb,Sa

def main(gF,cov,oF,ex,preF,sizeG,cpu,bypass):
    # check DB and program path
    info_P=str(pathlib.Path(__file__).parent.absolute())+'/programs.txt'
    Pfam=str(pathlib.Path(__file__).parent.absolute())+'/data/gtdbtk86/Pfam-A.hmm'
    Tigrfam=str(pathlib.Path(__file__).parent.absolute())+'/data/gtdbtk86/tigrfam.hmm'
    if not os.path.isfile(info_P):
        logging.error("Please check program.txt file !\n")
        sys.exit()
    with open(info_P,'r') as I:
        for i in I:
            i=i.rstrip('\n').split(':')
            if i[0]=='prodigal':
                prodigal=find_p(i[1])
            if i[0]=='hmmsearch':
                hmmsearch=find_p(i[1])
    if not os.path.isfile(Pfam) or not os.path.isfile(Tigrfam):
        logging.error("Please check the Pfam, Tigrfam DB path !\n")
        sys.exit()
    
    # running the program
    gF_P=os.path.abspath(gF)
    oF_P=os.path.abspath(oF)
    pathlib.Path(oF_P).mkdir(parents=True, exist_ok=True)

    bins=os.listdir(gF_P) 
    bins=filter(lambda x:x.endswith(ex),bins)
    Stat={}

    for i in bins:
        print (i)
        flag=0
        name=i.split('.'+ex)[0]
        bin_F=os.path.join(gF_P,i)
        bin_P=os.path.join(gF_P,name)
        #filtering MAGs size 500k
        statinfo = os.stat(bin_F)
        if float(statinfo.st_size)<sizeG:
            continue
       
        pathlib.Path(bin_P).mkdir(parents=True, exist_ok=True) 
        if os.path.isfile(bin_P+'/gene.Total.hit'):
            Sb,Sa=check_cluster(bin_P,ex,gF_P,oF_P,i,cov,preF,cpu)
        elif bypass == "N":
            flag=run_HMM_to_Marker(prodigal,hmmsearch,bin_F,bin_P,name,Pfam,Tigrfam,log_P=oF_P)
            if flag==1:
                continue
            Sb,Sa=check_cluster(bin_P,ex,gF_P,oF_P,i,cov,preF,cpu)
        try:
            Stat.update(Sb)
        except:
            continue
        try:
            Stat.update(Sa)
        except:
            continue

    df=pd.DataFrame(Stat).fillna(0)
    df.T.to_csv(oF_P+'/Bin_Stat.csv',index_label='Bin')

if __name__=="__main__":
    usage = """refine_bins.py -g [bin folder] -c [coverage file] -o [output]"""
    parser = optparse.OptionParser(usage)
    parser.add_option("-g","--genome",dest="genome",
        help="genome file path",type="string")
    parser.add_option("-o","--output",dest="output",
        help="output folder",type="string")
    parser.add_option("-c","--converage",dest="coverage",
        help="contig coverage file",type="string")
    parser.add_option("-e","--extention",dest="extension",
        help="genome file extension |default: fa",type="string")
    parser.add_option("-p","--prefix",dest="prefix",
        help="prefix |default: refine",type="string")
    parser.add_option("-s","--size",dest="size",
        help="MAG size filter |default: 500k",type="string")
    parser.add_option("-t","--thread",dest="thread",
        help="number of workers | default = 1",type="int")
    parser.add_option("-b","--bypass",dest="bypass",
        help="bypass prodigal - hmmsearch | default = N",type="string")
    (opts,args)=parser.parse_args()

    if opts.genome is None or opts.output is None or opts.coverage is None:
        parser.print_help()
    else:
        if opts.extension is None:
            opts.extension='fa'
        if opts.prefix is None:
            opts.prefix='refine'
        if opts.size is None:
            opts.size=500000
        if opts.thread is None:
            opts.thread=4
        if opts.bypass is None:
            opts.bypass="N"
        main(opts.genome,opts.coverage,opts.output,opts.extension,opts.prefix,opts.size,opts.thread,opts.bypass)
