import re
import sys
from collections import Counter
import itertools
import subprocess
import os
from time import sleep
import pandas as pd
import math
import numpy as np

def compute_MHC(hla_exome, hla_rna, epiopes):

    a=open(HLA_location)
    b=open(MHCI_list)
    c=open(MHCII_list)
    # Extract HLAs
    HLA_A=[]
    HLA_B=[]
    HLA_C=[]
    HLA_DQA1=[]
    HLA_DQB1=[]
    HLA_DRB1=[]


    for line in a:
        columns=line.split('\t')
        if columns[0] == mrn and line not in used_lines:
            # print('FOUND ONE!!')
            HLA_list=columns[3].split('_')
            if (HLA_list[1]=='A')and not re.search(r'locus', line)\
            and len(columns)>= 6:
                hla1=columns[4].split(':')
                hla2=columns[5].split(':')
                HLA1='HLA-'+str(':'.join(hla1[0:2]))
                HLA2='HLA-'+str(':'.join(hla2[0:2]))
                HLA_A.append(HLA1)
                HLA_A.append(HLA2)
            elif (HLA_list[1]=='B')and not re.search(r'locus', line)\
            and len(columns)>= 6:
                hla1=columns[4].split(':')
                hla2=columns[5].split(':')
                HLA1='HLA-'+str(':'.join(hla1[0:2]))
                HLA2='HLA-'+str(':'.join(hla2[0:2]))
                HLA_B.append(HLA1)
                HLA_B.append(HLA2)
            elif (HLA_list[1]=='C')and not re.search(r'locus', line)\
            and len(columns)>= 6:
                hla1=columns[4].split(':')
                hla2=columns[5].split(':')
                HLA1='HLA-'+str(':'.join(hla1[0:2]))
                HLA2='HLA-'+str(':'.join(hla2[0:2]))
                HLA_C.append(HLA1)
                HLA_C.append(HLA2)
            elif (HLA_list[1]=='DQA1')and not re.search(r'locus', line)\
            and len(columns)>= 6:
                hla1=columns[4].replace('*','').split(':')
                hla2=columns[5].replace('*','').split(':')
                HLA1=str(''.join(hla1[0:2]))
                HLA2=str(''.join(hla2[0:2]))
                HLA_DQA1.append(HLA1)
                HLA_DQA1.append(HLA2)
            elif (HLA_list[1]=='DQB1')and not re.search(r'locus', line)\
            and len(columns)>= 6:
                hla1=columns[4].replace('*','').split(':')
                hla2=columns[5].replace('*','').split(':')
                HLA1=str(''.join(hla1[0:2]))
                HLA2=str(''.join(hla2[0:2]))
                HLA_DQB1.append(HLA1)
                HLA_DQB1.append(HLA2)
            elif (HLA_list[1]=='DRB1')and not re.search(r'locus', line)\
            and len(columns)>= 6:
                hla1=columns[4].replace('*',':').split(':')
                hla2=columns[5].replace('*',':').split(':')
                HLA1=str(hla1[0])+'_'+str(hla1[1])+str(hla1[2])
                HLA2=str(hla2[0])+'_'+str(hla2[1])+str(hla2[2])
                HLA_DRB1.append(HLA1)
                HLA_DRB1.append(HLA2)
            used_lines.add(line)
    HLAs=[]
    HLADR=[]
    HLADQA1=[]
    HLADQB1=[]
    f = Counter(HLA_A)
    items=f.most_common()
    if len(items) ==1:
        cutoff=items[0][1]
    if len(items) >1 :
        cutoff=items[1][1]
    for x in items:
        if x[1]>= cutoff:
            HLAs.append(x[0])
    f = Counter(HLA_B)
    items=f.most_common()
    if len(items) ==1 :
        cutoff=items[0][1]
    if len(items) >1 :
        cutoff=items[1][1]
    for x in items:
        if x[1]>= cutoff:
            HLAs.append(x[0])
    f = Counter(HLA_C)
    items=f.most_common()
    if len(items) ==1 :
        cutoff=items[0][1]
    if len(items) >1 :
        cutoff=items[1][1]
    for x in items:
        if x[1]>= cutoff:
            HLAs.append(x[0])
    f = Counter(HLA_DQA1)
    items=f.most_common()
    if len(items) ==1 :
        cutoff=items[0][1]
    if len(items) >1 :
        cutoff=items[1][1]
    for x in items:
        if x[1]>= cutoff:
            HLADQA1.append(x[0])
    f = Counter(HLA_DQB1)
    items=f.most_common()
    if len(items) ==1 :
        cutoff=items[0][1]
    if len(items) >1 :
        cutoff=items[1][1]
    for x in items:
        if x[1]>= cutoff:
            HLADQB1.append(x[0])
    f = Counter(HLA_DRB1)
    items=f.most_common()
    if len(items) ==1 :
        cutoff=items[0][1]
    if len(items) >1 :
        cutoff=items[1][1]
    for x in items:
        if x[1]>= cutoff:
            HLADR.append(x[0])
    a.close()
    #########
    ##Check##
    #########
    if len(HLAs) <1:
        sys.exit('NO HLAS FOR {}'.format(mrn))
    #####################
    ##Create set of MHC##
    ##capable of being ##
    ##run              ##
    #####################
    MHC_I_list=set()
    ##create list of netmhcpan HLAs
    for line in b:
        MHC=line.strip()
        MHC_I_list.add(MHC)
    b.close()
    ##add to list from netmhcpanII HLAs
    MHC_II_list=set()
    for line in c:
        MHC=line.strip()
        MHC_II_list.add(MHC)
    c.close()
    ######################
    ##CREATE FASTA FILES##
    ##FOR NETMHC*       ##
    ######################
    if os.path.exists(loc+str(mrn)+'_overlap_final.txt'):
        print('Using '+str(mrn)+'_overlap_final.txt to generate FASTAs for netMHCpan.')
        df=pd.read_table(loc+str(mrn)+'_overlap_final.txt').fillna('-')
        epitope_fasta=open(loc+mrn+'_EpitopesI.fasta','w')
        epitope_fasta.close()
        epitope_fasta=open(loc+mrn+'_EpitopesII_15.fasta','w')
        epitope_fasta.close()
        epitope_fasta=open(loc+mrn+'_EpitopesII_14.fasta','w')
        epitope_fasta.close()
        epitope_fasta=open(loc+mrn+'_EpitopesII_13.fasta','w')
        epitope_fasta.close()
        epitope_fasta=open(loc+mrn+'_EpitopesII_12.fasta','w')
        epitope_fasta.close()
        mer_14=0
        mer_13=0
        mer_12=0
        for i,r in df.iterrows():
            index=i
            mer=r['Mut Epitope']
            Cosmic=r['Cosmic Info']
            if Cosmic==0:
                Cosmic='-'
            # print(r)
            if (len(mer) >= 8 and mer is not None and mer != '' and not re.search(r'\?', mer) and not re.search(r'X',mer))\
            and (r['Number of Exomes passing filters'] > 0 or re.search(r'ID', Cosmic)):
                epitope_fasta=open(loc+mrn+'_EpitopesI.fasta','a')
                epitope_fasta.write('>'+str(index)+'\n'+mer+'\n')
                epitope_fasta.close()
            if (len(mer) >= 15 and mer is not None and mer != '' and not re.search(r'\?', mer) and not re.search(r'X',mer))\
            and (r['Number of Exomes passing filters'] > 0 or re.search(r'ID', Cosmic)):
                epitope_fasta=open(loc+mrn+'_EpitopesII_15.fasta','a')
                epitope_fasta.write('>'+str(index)+'\n'+mer+'\n')
                epitope_fasta.close()
            elif (len(mer) == 14 and mer is not None and mer != '' and not re.search(r'\?', mer) and not re.search(r'X',mer))\
            and (r['Number of Exomes passing filters'] > 0 or re.search(r'ID', Cosmic)):
                epitope_fasta=open(loc+mrn+'_EpitopesII_14.fasta','a')
                epitope_fasta.write('>'+str(index)+'\n'+mer+'\n')
                epitope_fasta.close()
                mer_14=mer_14+1
            elif (len(mer) == 13 and mer is not None and mer != '' and not re.search(r'\?', mer) and not re.search(r'X',mer))\
            and (r['Number of Exomes passing filters'] > 0 or re.search(r'ID', Cosmic)):
                epitope_fasta=open(loc+mrn+'_EpitopesII_13.fasta','a')
                epitope_fasta.write('>'+str(index)+'\n'+mer+'\n')
                epitope_fasta.close()
                mer_13=mer_13+1
            elif (len(mer) == 12 and mer is not None and mer != '' and not re.search(r'\?', mer) and not re.search(r'X',mer))\
            and (r['Number of Exomes passing filters'] > 0 or re.search(r'ID', Cosmic)):
                epitope_fasta=open(loc+mrn+'_EpitopesII_12.fasta','a')
                epitope_fasta.write('>'+str(index)+'\n'+mer+'\n')
                epitope_fasta.close()
                mer_12=mer_12+1
        print('The number of 12mers is {}\nThe number of 13mers is {}\nThe number of 14mers is {}'.format(mer_12, mer_13,mer_14))
        df.to_csv(loc+str(mrn)+'_index2fasta.txt', sep='\t', index=True, index_label='Index')
    ##########################
    ##netmhc not recognizing##
    ##last line in fasta had##
    ##to do this            ##
    ##########################
    if os.path.exists(loc+mrn+'_EpitopesI.fasta'):
        epitope_fasta=open(loc+mrn+'_EpitopesI.fasta','a')
        epitope_fasta.write('\n')
        epitope_fasta.close()
    if os.path.exists(loc+mrn+'_EpitopesII_15.fasta'):
        epitope_fasta=open(loc+mrn+'_EpitopesII_15.fasta','a')
        epitope_fasta.write('\n')
        epitope_fasta.close()
    if os.path.exists(loc+mrn+'_EpitopesII_14.fasta'):
        epitope_fasta=open(loc+mrn+'_EpitopesII_14.fasta','a')
        epitope_fasta.write('\n')
        epitope_fasta.close()
    if os.path.exists(loc+mrn+'_EpitopesII_13.fasta'):
        epitope_fasta=open(loc+mrn+'_EpitopesII_13.fasta','a')
        epitope_fasta.write('\n')
        epitope_fasta.close()
    if os.path.exists(loc+mrn+'_EpitopesII_13.fasta'):
        epitope_fasta=open(loc+mrn+'_EpitopesII_13.fasta','a')
        epitope_fasta.write('\n')
        epitope_fasta.close()
    ##########################
    ##Going torun each HLA  ##
    ##prediction as part of ##
    ##a swarm file on our   ##
    ##cluster to save time  ##
    ##########################
    waitfile=open(loc+mrn+'_wait_file.txt','a')
    waitfile.write('This is a file it\'s presences signifies that the netmhc predictions haven\'t finished\n')
    waitfile.close()
    swarmfile=open(loc+mrn+'.netMHC.swarm', 'a')
    swarmfile.write('##Swarmfile to run a fasta file through netMHC for multiple HLAs\n')
    swarmfile.close()
    i=0
    failed_HLAs=set()
    for a in HLAs:
        y=open(loc+mrn+'_HLA_list.txt', 'a')
        y.write(a+'\n')
        y.close()
        ##need to format HLAs to see if they are in the list of MHC for netmhc4.0 or Pan
        HLA_pan=a.replace('*','')
        HLA=HLA_pan.replace(':','')
        # HLA=HLA_pan.replace(':','')
        print('searching for '+HLA_pan+' in netMHCpan list of HLAs.')
        if HLA_pan in MHC_I_list:
            print('found '+HLA_pan+' in netMHCpan list of HLAs.')
            swarmfile=open(loc+mrn+'.netMHC.swarm', 'a')
            swarmfile.write(str(netmhcpan)+' '+loc+mrn+'_EpitopesI.fasta -a '+HLA_pan+' -l 8,9,10,11,12 -s -xls  -xlsfile '+loc+mrn+'_netMHCpan_'+HLA+'.txt\n')
            swarmfile.close()
            i=i+1
        else:
            print(a+' not found in netMHC list of HLAs.')
            failed_HLAs.add(HLA_pan)
    #################################################
    ## need to make possible iterations of all DQ  ##
    ## combinations                                ##
    #################################################
    HLA_DQII=set()
    for x in itertools.product(HLADQA1,HLADQB1):
        HLA_pan='HLA-'+x[0]+'-'+x[1]
        HLA_DQII.add(HLA_pan)
    for HLA_pan in HLA_DQII:
        print('looking for '+HLA_pan+'in netHCIIpan list')
        if HLA_pan in MHC_II_list:
            print('found '+HLA_pan+' in netMHCIIpan list of HLAs.')
            swarmfile=open(loc+mrn+'.netMHC.swarm', 'a')
            swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_15.fasta -a '+HLA_pan+' -length 15 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_15_15_'+HLA_pan+'.txt\n')
            swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_15.fasta -a '+HLA_pan+' -length 14 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_15_14_'+HLA_pan+'.txt\n')
            swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_15.fasta -a '+HLA_pan+' -length 13 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_15_13_'+HLA_pan+'.txt\n')
            swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_15.fasta -a '+HLA_pan+' -length 12 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_15_12_'+HLA_pan+'.txt\n')
            swarmfile.close()
            i=i+4
            if mer_14 >0 :
                swarmfile=open(loc+mrn+'.netMHC.swarm', 'a')
                swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_14.fasta -a '+HLA_pan+' -length 14 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_14_14_'+HLA_pan+'.txt\n')
                swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_14.fasta -a '+HLA_pan+' -length 13 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_14_13_'+HLA_pan+'.txt\n')
                swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_14.fasta -a '+HLA_pan+' -length 12 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_14_12_'+HLA_pan+'.txt\n')
                swarmfile.close()
                i=i+3
            if mer_13 >0 :
                swarmfile=open(loc+mrn+'.netMHC.swarm', 'a')
                swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_13.fasta -a '+HLA_pan+' -length 13 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_13_13_'+HLA_pan+'.txt\n')
                swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_13.fasta -a '+HLA_pan+' -length 12 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_13_12_'+HLA_pan+'.txt\n')
                swarmfile.close()
                i=i+2
            if mer_12 >0 :
                swarmfile=open(loc+mrn+'.netMHC.swarm', 'a')
                swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_12.fasta -a '+HLA_pan+' -length 12 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_12_12_'+HLA_pan+'.txt\n')
                swarmfile.close()
                i=i+1
        else:
            print(a+' not found in netMHCIIpan list of HLAs.')
            failed_HLAs.add(HLA_pan)
    ############################################
    ## need to add the DR samples and use the ##
    ## netmhcII predictions                   ##
    ############################################
    for x in HLADR:
        HLA_pan=x
        print('looking for '+HLA_pan+'in netHCIIpan list')
        if HLA_pan in MHC_II_list:
            print('found '+HLA_pan+' in netMHCIIpan list of HLAs.')
            swarmfile=open(loc+mrn+'.netMHC.swarm', 'a')
            swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_15.fasta -a '+HLA_pan+' -length 15 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_15_'+HLA_pan.replace('_','-')+'.txt\n')
            swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_15.fasta -a '+HLA_pan+' -length 14 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_15_14_'+HLA_pan.replace('_','-')+'.txt\n')
            swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_15.fasta -a '+HLA_pan+' -length 13 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_15_13_'+HLA_pan.replace('_','-')+'.txt\n')
            swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_15.fasta -a '+HLA_pan+' -length 12 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_15_12_'+HLA_pan.replace('_','-')+'.txt\n')
            swarmfile.close()
            i=i+4
            if mer_14 >0 :
                swarmfile=open(loc+mrn+'.netMHC.swarm', 'a')
                swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_14.fasta -a '+HLA_pan+' -length 14 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_14_14_'+HLA_pan.replace('_','-')+'.txt\n')
                swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_14.fasta -a '+HLA_pan+' -length 13 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_14_13_'+HLA_pan.replace('_','-')+'.txt\n')
                swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_14.fasta -a '+HLA_pan+' -length 12 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_14_12_'+HLA_pan.replace('_','-')+'.txt\n')
                swarmfile.close()
                i=i+3
            if mer_13 >0 :
                swarmfile=open(loc+mrn+'.netMHC.swarm', 'a')
                swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_13.fasta -a '+HLA_pan+' -length 13 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_13_13_'+HLA_pan.replace('_','-')+'.txt\n')
                swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_13.fasta -a '+HLA_pan+' -length 12 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_13_12_'+HLA_pan.replace('_','-')+'.txt\n')
                swarmfile.close()
                i=i+2
            if mer_12 >0 :
                swarmfile=open(loc+mrn+'.netMHC.swarm', 'a')
                swarmfile.write(str(netmhcIIpan)+' -f '+loc+mrn+'_EpitopesII_12.fasta -a '+HLA_pan+' -length 12 -s -xls  -xlsfile '+loc+mrn+'_netMHCIIpan31_12_12_'+HLA_pan.replace('_','-')+'.txt\n')
                swarmfile.close()
                i=i+1
        else:
            print(a+' not found in netMHCIIpan list of HLAs.')
            failed_HLAs.add(HLA_pan)