import sys,os
import shutil
import pathlib
import datetime
import optparse
import logging
import subprocess
import collections
import json
import pandas as pd
import acr_utils.cluster_bin as cb
import acr_utils.NRmarker as NR
from Bio.SeqIO.FastaIO import SimpleFastaParser as SFP
from glob import glob

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

def run_cmd2(cmd,log_dir,shell=False,program=False,check=True):
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
    try:
        with open(Out,'a') as sto, open(Err,'a') as ste:
            subprocess.run(cmd, check=check, shell=shell, stdout=sto, stderr=ste)
    except subprocess.CalledProcessError:
        if program:
            logging.error("an error occured while executing {}".format(program))
        else:
            logging.error("an error occured while executing")
        exit(1)

def write_log(log_dir,message):
    log_P=os.path.abspath(log_dir)
    Out = os.path.join(log_P+'/acr.log.out')
    with open(Out, 'a') as L:
        L.write(message+'\n')

def find_p(name):
    # check the program path
    loc = shutil.which(name)
    if loc != None:
        try:
            help_m='-h'
            if 'pplacer' in name or 'guppy' in name:
                help_m='--help'
            if 'guppy' in name:
                o=subprocess.check_output([loc,'tog',help_m],stderr=subprocess.STDOUT,shell=False,universal_newlines=True,encoding='utf-8')
            else:
                o=subprocess.check_output([loc,help_m],stderr=subprocess.STDOUT,shell=False,universal_newlines=True,encoding='utf-8')
            return loc
        except:
            logging.error("Please check the {0} !\n".format(loc)) 
            sys.exit()
    else:
        logging.error("Please set the {0} path !\n".format(loc))
        sys.exit()

def filter_overlapped(evalue,strand,cds_f):
    tmp={k:v for k,v in evalue.items() if k in strand}
    cds_d=collections.defaultdict(list)
    for i in tmp:
        ctg='_'.join(i.split('_')[:-1])
        order=int(i.split('_')[-1])
        cds_d[ctg].append(order)
    for i in cds_d:
        flag=0
        new_order=0
        cds_group=collections.defaultdict(set)
        cds_d[i].sort()
        r_cds=len(cds_d[i])
        for j in range(r_cds):
            if j==r_cds-1:
                break
            if cds_d[i][j+1]-cds_d[i][j] <= cds_f:
                if flag==0:
                    new_order+=1
                cds_group[new_order].add(cds_d[i][j])
                cds_group[new_order].add(cds_d[i][j+1])
                flag=1
            else:
                flag=0
        for j in cds_group:
            if len(cds_group[j])<2:
                continue
            bit_max=0
            for k in cds_group[j]:
                bit_max2=max(set(tmp[i+'_'+str(k)].keys()))
                if bit_max2>bit_max:
                    bit_max=bit_max2
                else:
                    del(evalue[i+'_'+str(k)])
    return evalue

def score_sort(blast_out,blast_filter,evalue_f=False,cds_f=False):
    temp=set()
    evalue={}
    strand=collections.defaultdict(set)
    with open(blast_out,'r') as B, open(blast_filter,'w') as N:
        while True:
            line=B.readline()
            if not line: break
            if not line.startswith('#'):
                spl=line.rstrip('\n').split()
                temp.add(spl[0])
                gID=spl[3]
                if spl[3]=='-':
                    gID=spl[2]
                if evalue_f:
                    if float(spl[4]) < evalue_f:
                        evalue.setdefault(spl[0],{})
                        evalue[spl[0]].setdefault(float(spl[5]),gID)
                        if cds_f:
                            strand[spl[23]].add(spl[0])
                else:
                    evalue.setdefault(spl[0],{})
                    evalue[spl[0]].setdefault(float(spl[5]),gID)
                    if cds_f:
                        strand[spl[23]].add(spl[0])
        # near CDS filtering (same strand with in 1 gene distance)
        if cds_f:
            evlaue=filter_overlapped(evalue,strand['1'],cds_f) # +strand filtering
            evlaue=filter_overlapped(evalue,strand['-1'],cds_f) # -strand filtering

        # get max bit score
        for i in evalue:
            evalue[i]=evalue[i][max(set(evalue[i].keys()))]
        B.seek(0)
        while True:
            line=B.readline()
            if not line: break
            if line.startswith('#'): continue
            spl=line.rstrip('\n').split()
            gID=spl[3]
            if spl[3]=='-':
                gID=spl[2]
            if spl[0] in evalue and gID == evalue[spl[0]]:
                N.write(line)

def run_prodigal(prodigal,bin_F,bin_P,name,log_P):
    cmd=[prodigal,'-i',bin_F,'-a',bin_P+'/'+name+'.faa','-d',bin_P+'/'+name+'.fna','-f','gff','-o',bin_P+'/'+name+'.gff']
    run_cmd2(cmd,program='prodigal',log_dir=log_P)
    statinfo = os.stat(bin_P+'/'+name+'.faa')
    if float(statinfo.st_size)==0:
        write_log(log_P,"# faa file is not exist")
        return 1
    return 0

def run_HMM_to_Marker(hmmsearch,process,bin_P,name,Pfam,Tigrfam,Eukcc_DB,log_P):
    #cmd=[hmmsearch,'--cpu',str(process),'--tblout',bin_P+'/gene.P.hit','--cut_ga',Pfam,bin_P+'/'+name+'.faa','>',bin_P+'/gene.Pfam']
    #run_cmd(' '.join(cmd),program='hmmsearch',log_dir=log_P,shell=True)
    cmd=[hmmsearch,'--cpu',str(process),'-o',bin_P+'/gene.Pfam','--tblout',bin_P+'/gene.P.hit','--cut_ga',Pfam,bin_P+'/'+name+'.faa']
    run_cmd2(cmd,program='hmmsearch',log_dir=log_P,shell=False)
    score_sort(bin_P+'/gene.P.hit',bin_P+'/gene.f.P.hit')

    #cmd=[hmmsearch,'--cpu',str(process),'--tblout',bin_P+'/gene.T.hit','--cut_ga',Tigrfam,bin_P+'/'+name+'.faa','>',bin_P+'/gene.Tfam']
    #run_cmd(' '.join(cmd),program='hmmsearch',log_dir=log_P,shell=True)
    cmd=[hmmsearch,'--cpu',str(process),'-o',bin_P+'/gene.Tfam','--tblout',bin_P+'/gene.T.hit','--cut_ga',Tigrfam,bin_P+'/'+name+'.faa']
    run_cmd2(cmd,program='hmmsearch',log_dir=log_P,shell=False)
    score_sort(bin_P+'/gene.T.hit',bin_P+'/gene.f.T.hit')

    if Eukcc_DB:
        #cmd=[hmmsearch,'--cpu',str(process),'--tblout',bin_P+'/gene.E.hit','--cut_ga',Eukcc_DB,bin_P+'/'+name+'.faa','>',bin_P+'/gene.Efam']
        #run_cmd(' '.join(cmd),program='hmmsearch',log_dir=log_P,shell=True)
        cmd=[hmmsearch,'--cpu',str(process),'-o',bin_P+'/gene.Efam','--tblout',bin_P+'/gene.E.hit','--cut_ga',Eukcc_DB,bin_P+'/'+name+'.faa']
        run_cmd2(cmd,program='hmmsearch',log_dir=log_P,shell=False)
        score_sort(bin_P+'/gene.E.hit',bin_P+'/gene.f.E.hit',evalue_f=float(1e-5),cds_f=1)
        
        cmd=['cat',bin_P+'/gene.f.P.hit',bin_P+'/gene.f.T.hit',bin_P+'/gene.f.E.hit'] 
        with open(bin_P+'/gene.Total.hit','w') as sto:
            proc=subprocess.Popen(cmd,stdout=sto)#,stderr=subprocess.PIPE)
            proc.communicate() 
    else:
        cmd=['cat',bin_P+'/gene.f.P.hit',bin_P+'/gene.f.T.hit']#,'>',bin_P+'/gene.Total.hit'] 
        with open(bin_P+'/gene.Total.hit','w') as sto:
            proc=subprocess.Popen(cmd,stdout=sto)#,stderr=subprocess.PIPE)
            proc.communicate()

def reduceJplace(jplace, jplaceselection, placementCutoff=0.4):
    """
    reduces a jplace file to placements with at least placementCutoff post_prob
     from EukCC with default Cutoff
    """
    with open(jplace) as json_file:
        j = json.load(json_file)
        fields = j["fields"]
        newplacements = []
        for placement in j["placements"]:
            ps = placement["p"]
            kp = []
            for p in ps:
                d = {k: v for k, v in zip(fields, p)}
                if d["post_prob"] >= placementCutoff:
                    kp.append(p)

            if len(kp) > 0:
                placementn = placement.copy()
                placementn["p"] = kp
                newplacements.append(placementn)
        jn = j.copy()
        jn["placements"] = newplacements

    with open(jplaceselection, "w") as f:
        json.dump(jn, f)
    return jplaceselection

def find_euk_Marker(hmmsearch,hmmalign,pplacer,guppy,process,bin_P,name,DB,refpkg,Eukcc_core_list,Eukcc_core_tree,Eukcc_marker_setinfo,log_P):
    cmd=[hmmsearch,'--cpu',str(process),'-o',bin_P+'/gene.Ecfam','--tblout',bin_P+'/gene.Ec.hit','--cut_ga',DB,bin_P+'/'+name+'.faa']
    #run_cmd(' '.join(cmd),program='hmmsearch',log_dir=log_P,shell=True)
    run_cmd2(cmd,program='hmmsearch',log_dir=log_P,shell=False)
    score_sort(bin_P+'/gene.Ec.hit',bin_P+'/gene.f.Ec.hit')
    top_hit={}
    # get order of SCGs
    write_log(log_P,"# find euk core genes")
    cOrder=[]
    with open(Eukcc_core_list,'r') as C:
        for i in C:
            cOrder.append(i.rstrip('\n'))
    gOrder=[]
    # make tmp folder
    tmp_P=bin_P+'/tmp/'
    pathlib.Path(tmp_P).mkdir(parents=True, exist_ok=True)
    with open(bin_P+'/gene.f.Ec.hit','r') as H,\
        open(bin_P+'/'+name+'.faa','r') as F:
        for i in H:
            if i.startswith('#'):
                continue
            i=i.split()
            gID=i[2].split('.')[0]
            if not gID in gOrder:
                top_hit[i[0]]=gID
                gOrder.append(gID)
        for t,s in SFP(F):
            if t.split()[0] in top_hit:
                with open(tmp_P+'/'+top_hit[t.split()[0]]+'.faa','w') as N:
                    N.write('>'+name+'\n'+s)
    # align and concatanate
    with open(tmp_P+'/concat_aln.fasta','w') as A:
        cLength={}
        gSeq=collections.defaultdict(dict)
        for g in gOrder:
            cmd=[hmmalign,'-o',tmp_P+'/'+g+'.aln','--outformat','afa','--mapali',refpkg+'/'+g+'.refpkg/'+g+'.fasta','--amino',refpkg+'/'+g+'.refpkg/'+g+'.hmm',tmp_P+'/'+g+'.faa']
            #run_cmd(' '.join(cmd),program='hmmalign',log_dir=log_P,shell=True)
            run_cmd2(cmd,program='hmmalign',log_dir=log_P,shell=False)
            with open(tmp_P+'/'+g+'.aln','r') as T:
                for t,s in SFP(T):
                    gSeq[t][g]=s.replace(".","-").replace("*","-")
                    if not g in cLength:
                        cLength[g]=len(s)
        # get length of alignment of core genes
        for i in glob(refpkg+'/P*refpkg/*.fasta'):
            gID=i.split('/')[-1].split('.fasta')[0]
            with open(i,'r') as F:
                for t,s in SFP(F):
                    if not gID in cLength:
                        cLength[gID]=len(s)
                    break
        # write concat align fasta
        for g in gSeq:
            oSeq=''
            for c in cOrder:
                if c in gSeq[g]:
                    oSeq+=gSeq[g][c]
                else:
                    oSeq+='-'*cLength[c]
            A.write('>'+g+'\n'+oSeq+'\n')
    # remove core hmm results
    for c in glob(tmp_P+'/PTHR*'):
        os.remove(c)
    # run pplacer
    cmd=[pplacer,'-p','--keep-at-most','5','-m','LG','-j',str(process),'-o',tmp_P+'/concat.pp','-c',refpkg+'/concat.refpkg',tmp_P+'/concat_aln.fasta']
    #run_cmd(' '.join(cmd),program='pplacer',log_dir=log_P,shell=True)
    run_cmd2(cmd,program='pplacer',log_dir=log_P,shell=False,check=False)
    if not os.path.isfile(tmp_P+'/concat.pp'):
        write_log(log_P,"could not find euk genome signal...")
        return None 
    # to placements with al least placementCutoff (0.4
    reduceJplace(tmp_P+'/concat.pp',tmp_P+'/concat.pp.jplace',0.4)
    # run guppy
    cmd=[guppy,'tog',tmp_P+'/concat.pp.jplace','-o',tmp_P+'/placement.tree']
    #run_cmd(' '.join(cmd),program='guppy',log_dir=log_P,shell=True)
    run_cmd2(cmd,program='guppy',log_dir=log_P,shell=False)
    if not os.path.isfile(tmp_P+'/placement.tree'):
        write_log(log_P,"could not make euk tree...") 
        return None
    # get the nearest markerset
    markers=NR.find_N_marker(Eukcc_core_tree,tmp_P+'/placement.tree',Eukcc_marker_setinfo)
    if markers:
        if name==markers[0]: 
            marker=markers[1]
            sisters=markers[2]
            print ("find the best marker for eukaryotic MAG: ",marker,sisters)
            write_log(log_P,"marker of Eukcc: "+marker)
            return marker
    

def check_cluster(bin_P,ex,gF_P,oF_P,i,cov,preF,process,jgi,Eukcc_DB,log_P):
    BSCGs={};ASCGs={};ESCGs={}
    BCheck=[];ACheck=[];ECheck=[]
    with open(bin_P+'/gene.Total.hit','r') as B:
        for l in B:
            if not l.startswith('#'):
                l=l.split()
                BCheck,BSCGs=cb.make_Marker(l,cb.BAC120_MARKERS,BCheck,BSCGs)
                ACheck,ASCGs=cb.make_Marker(l,cb.AR122_MARKERS,ACheck,ASCGs)
    
    BBinStat=collections.Counter(BCheck)
    ABinStat=collections.Counter(ACheck)
    Sb=cb.check_Marker(gF_P,ex,oF_P,i,cov,BBinStat,BSCGs,cb.BAC120_MARKERS,50,'Bac',preF,process,jgi,log_P)
    Sa=cb.check_Marker(gF_P,ex,oF_P,i,cov,ABinStat,ASCGs,cb.AR122_MARKERS,50,'Arc',preF,process,jgi,log_P) 
    if Eukcc_DB:
        EUK_MARKERS=cb.select_E_marker(Eukcc_DB.split('/')[-1].split('.')[0])
        with open(bin_P+'/gene.Total.hit','r') as B:
            for l in B:
                if not l.startswith('#'):
                    l=l.split()
                    ECheck,ESCGs=cb.make_Marker(l,EUK_MARKERS,ECheck,ESCGs)
        EBinStat=collections.Counter(ECheck)
        Se=cb.check_Marker(gF_P,ex,oF_P,i,cov,EBinStat,ESCGs,EUK_MARKERS,50,'Euk',preF,process,jgi,log_P) #filter_Q>50
        return Sb,Sa,Se
    else:
        return Sb,Sa

def main(gF,cov,oF,ex,preF,sizeG,process,bypass,jgi):
    # check DB and program path
    info_P=str(pathlib.Path(__file__).parent.absolute())+'/programs.txt'
    Pfam=str(pathlib.Path(__file__).parent.absolute())+'/data/gtdbtk86/Pfam-A.hmm'
    Tigrfam=str(pathlib.Path(__file__).parent.absolute())+'/data/gtdbtk86/tigrfam.hmm'
    Eukcc_core=str(pathlib.Path(__file__).parent.absolute())+'/data/eukcc_db_20191023/concat.hmm'
    Eukcc_core_tree=str(pathlib.Path(__file__).parent.absolute())+'/data/eukcc_db_20191023/refpkg/concat.refpkg/concat.tree'
    Eukcc_marker_setinfo=str(pathlib.Path(__file__).parent.absolute())+'/data/eukcc_db_20191023/sets/setinfo.csv'
    Eukcc_marker_set=str(pathlib.Path(__file__).parent.absolute())+'/data/eukcc_db_20191023/hmm_sets/sets/'
    Eukcc_core_info=str(pathlib.Path(__file__).parent.absolute())+'/data/eukcc_db_20191023/refpkg/'
    Eukcc_core_list=str(pathlib.Path(__file__).parent.absolute())+'/data/eukcc_db_20191023/profile.list'

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
            if i[0]=='EukRep':
                EukRep=find_p(i[1])
            if i[0]=='hmmalign':
                hmmalign=find_p(i[1])
            if i[0]=='pplacer':
                pplacer=find_p(i[1])
            if i[0]=='guppy':
                guppy=find_p(i[1])
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
        Eukcc_DB=False
        flag=0
        name=i.split(ex)[0]
        bin_F=os.path.join(gF_P,i)
        bin_P=os.path.join(gF_P,name)
        #filtering MAGs size 500k
        statinfo = os.stat(bin_F)
        if float(statinfo.st_size)<sizeG:
            continue
        #check euk MAGs?
        statinfo = os.stat(bin_F)
        if float(statinfo.st_size)>=5000000:
            cmd=[EukRep,"-i",bin_F,"-o",gF_P+"/eukrep.tmp","-ff"]
            #proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding='utf8')
            #run_cmd(' '.join(cmd),program='eukrep',log_dir=oF_P,shell=True)
            run_cmd2(cmd,program='eukrep',log_dir=oF_P,shell=False)
            if os.stat(gF_P+"/eukrep.tmp").st_size >=500000:
                write_log(oF_P,"# Pass Eukrep for a MAG")
                Eukcc_DB=Eukcc_marker_set
            os.remove(gF_P+"/eukrep.tmp")

        #if not Eukcc_DB: #tmp skip by hoonje
        #    continue

        # make output path
        pathlib.Path(bin_P).mkdir(parents=True, exist_ok=True)
        # run prodigal
        if run_prodigal(prodigal,bin_F,bin_P,name,log_P=oF_P)==1:
            continue
        # find nearest the euk marker
        if Eukcc_DB:
            Eukcc_DB=find_euk_Marker(hmmsearch,hmmalign,pplacer,guppy,process,bin_P,name,Eukcc_core,Eukcc_core_info,Eukcc_core_list,Eukcc_core_tree,Eukcc_marker_setinfo,log_P=oF_P)
            if Eukcc_DB:
                Eukcc_DB=os.path.join(Eukcc_marker_set,Eukcc_DB+'.hmm')
            else:
                Eukcc_DB=False
        # bypass
        if os.path.isfile(bin_P+'/gene.Total.hit'):
           bypass='Y'
        # run hmmsearch to marker set
        if bypass == "N":
            run_HMM_to_Marker(hmmsearch,process,bin_P,name,Pfam,Tigrfam,Eukcc_DB,log_P=oF_P)
        # bypass or run
        if not os.path.isfile(bin_P+'/gene.Total.hit'):
            write_log(oF_P,"# There are no matches for marker sets")
            continue
        # get SCGs data
        if Eukcc_DB:
            write_log(oF_P,"# start refining bac/arc/euk MAGs")
            Sb,Sa,Se=check_cluster(bin_P,ex,gF_P,oF_P,i,cov,preF,process,jgi,Eukcc_DB,log_P=oF_P)
        else:
            write_log(oF_P,"# start refining bac/arc MAGs")
            Sb,Sa=check_cluster(bin_P,ex,gF_P,oF_P,i,cov,preF,process,jgi,Eukcc_DB,log_P=oF_P)
        try:
            Stat.update(Sb)
        except:
            continue
        try:
            Stat.update(Sa)
        except:
            continue
        try:
            Stat.update(Se)
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
    parser.add_option("-c","--coverage",dest="coverage",
        help="contig coverage file",type="string")
    parser.add_option("-e","--extention",dest="extension",
        help="genome file extension |default: fa",type="string")
    parser.add_option("-p","--prefix",dest="prefix",
        help="prefix |default: refine",type="string")
    parser.add_option("-s","--size",dest="size",
        help="MAG size filter |default: 500k",type="string")
    parser.add_option("-t","--process",dest="process",
        help="number of workers | default = 1",type="int")
    parser.add_option("-b","--bypass",dest="bypass",
        help="bypass prodigal - hmmsearch | default = N",type="string")
    parser.add_option("-j","--from_jgi_cov",dest="jgi",
        help="please insert the number of samples (using jgi coverage file) | default = N",type="string")
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
        if opts.process is None:
            opts.process=4
        if opts.bypass is None:
            opts.bypass="N"
        if opts.jgi is None:
            opts.jgi="N"
        main(opts.genome,opts.coverage,opts.output,'.'+opts.extension,opts.prefix,opts.size,opts.process,opts.bypass,opts.jgi)
