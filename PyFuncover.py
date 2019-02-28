#! /usr/bin/ python
# -*- coding: utf8 -*-

###########
# IMPORTS #
###########

from __future__ import division
import argparse
from Bio import SearchIO, SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
import ftplib
from ftplib import FTP
import gzip
import json
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
import os
import pandas as pd
import re
import requests
import shutil
import subprocess
import tarfile
from threading import Thread
import urllib
import textwrap
import sys

#get system
if sys.platform in ['win32']:
    FLAG = True
else:
    FLAG = False

#disable BioPython Warning about future deprecated Search.IO
import warnings
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO

#############
# FUNCTIONS #
#############

def init():
    """
    Download genomes references and taxonomic database
    """
    
    print("Connection to the FTP NCBI Serveur...")
    ftp = ftplib.FTP()
    ftp.connect('ftp.ncbi.nlm.nih.gov')
    ftp.login(user='anonymous',passwd='your.email@gmail.com')
    
    print('Download the reference genomes list...')
    with open(os.getcwd()+"/assembly_summary_refseq.txt", "w") as f:
        ftp.retrbinary('RETR %s' % 'genomes/refseq/assembly_summary_refseq.txt', f.write)
        
    print('Download the genomes assembly list for eukaryotes organisms...')
    with open(os.getcwd()+"/eukaryotes.txt", "w") as f:
        ftp.retrbinary('RETR %s' % 'genomes/GENOME_REPORTS/eukaryotes.txt', f.write)
        
    print('Download the genomes assembly list for prokaryotes organisms...')
    with open(os.getcwd()+"/prokaryotes.txt", "w") as f:
        ftp.retrbinary('RETR %s' % 'genomes/GENOME_REPORTS/prokaryotes.txt', f.write)
    
    ftp.close()
    #check if folder dont already exist
    if os.path.isdir(os.getcwd()+"/TAXONOMY") == False:
        os.mkdir(os.getcwd()+"/TAXONOMY")
    print('Download the Taxonomic database...')
    with open(os.getcwd()+"/TAXONOMY/taxdump.tar.gz", "w") as f:
        urllib.urlretrieve("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", os.getcwd()+"/TAXONOMY/taxdump.tar.gz")
    urllib.urlcleanup()

    print('Extracting Database...')
    tfile = tarfile.open(os.getcwd()+"/TAXONOMY/taxdump.tar.gz", 'r:gz')
    tfile.extract("nodes.dmp",os.getcwd()+"/TAXONOMY/.")
    tfile.extract("names.dmp",os.getcwd()+"/TAXONOMY/.")
    tfile.close()
    os.remove(os.getcwd()+"/TAXONOMY/taxdump.tar.gz")
    print('Finish ! Exit')

def dic_Name():
    """
    Create the dicName and dicNameInverse
    
    dicName :
        key : Taxid : value : specie (scientific name)
    dicNameInverse :
        Key : Species (scientific name) : Value : Taxid

    return : dicName, dicNameInverse
    """
    #open file    
    names = open(os.getcwd()+"/TAXONOMY/names.dmp","r")
    #create dict with names
    dicName = dict()
    dicNameInverse = dict()
    for ligne in names:
        ligne = ligne.replace("\t|\n","").split("\t|\t")
        fin = ligne[len(ligne)-1]
        if fin == 'scientific name':
            dicName[ligne[0]]=ligne[1]
            dicNameInverse[ligne[1]]=ligne[0]
    return dicName, dicNameInverse

def dic_Node():
    """
    return dicNode
    """
    #open file
    nodes = open(os.getcwd()+"/TAXONOMY/nodes.dmp","r")
    #create a dict with nodes
    dicNode = dict()
    for ligne in nodes:
        ligne = ligne.replace("\t|\n","").split("\t|\t")
        dicNode[ligne[0]] = [ligne[1],ligne[2]]
    return dicNode

def getLineage(specie,dic_Node,dic_Name):
    #create dicNode and dicName
    dicNode = dic_Node
    dicName = dic_Name[0]
    dicNameInverse = dic_Name[1]
    sp = specie
    for j in "[]":
        sp = sp.replace(j,"")
    i = dicNameInverse[sp]
    genre = [specie]
    id_type = dicNode[i][1] #controle la boucle while si "root"
    id_parent = dicNode[i][0]
    while genre[0] != "root":
        parent = [dicName[id_parent]]
        genre = parent + genre
        id_type = dicNode[id_parent][1] #met a jour le controle
        id_parent = dicNode[id_parent][0]
    return genre

def lineage(dicName,dicNode):
    """
    Create a dictionnary with all Lineage from self.nodes and self.name
    """
    dicLineage = dict()
    for i in dicNode.keys():
        genre = [dicName[i]]
        id_type = dicNode[i][1] #controle la boucle while si "root"
        id_parent = dicNode[i][0]
        while genre[0] != "root":
            parent = [dicName[id_parent]]
            genre = parent + genre
            id_type = dicNode[id_parent][1] #met a jour le controle
            id_parent = dicNode[id_parent][0]
        dicLineage[genre[-1]]=genre
    return dicLineage

def groupedTaxId(dicLineage, dicNameinverse):
    """
    Create a dictionnary with all sub taxid for a given taxid
    """
    dicGroupedTaxid = dict()
    for species in dicLineage.keys():
        lineage_liste = dicLineage[species]
        taxid = dicNameinverse[species]
        dicGroupedTaxid[taxid] = set([taxid])
        for sub_rank in lineage_liste:
            taxid_rank = dicNameinverse[sub_rank]
            if taxid_rank not in dicGroupedTaxid:
                dicGroupedTaxid[taxid_rank] = set()
                dicGroupedTaxid[taxid_rank].add(taxid)
            else:
                dicGroupedTaxid[taxid_rank].add(taxid)
    return dicGroupedTaxid

def chunks(l, n):
    for i in xrange(0, len(l), n):
        yield l[i:i + n]

def uniprot2ncbi_accession(l):
    """
    l : liste of UNIPROT Accessions
    """
    dico_data = dict()
    
    headers = {
        'Connection': 'keep-alive',
        'Pragma': 'no-cache',
        'Cache-Control': 'no-cache',
        'Upgrade-Insecure-Requests': '1',
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.110 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8',
        'Accept-Encoding': 'gzip, deflate, br',
        'Accept-Language': 'fr-FR,fr;q=0.9,en-US;q=0.8,en;q=0.7',
    }

    params = {
        'method': 'db2db',
        'input': 'UniProt Accession',
        'inputValues': '',
        'outputs': 'GenBank Protein Accession',
        'taxonId': '',
        'format': 'row',
    }


    lol_accession = list(chunks(list(l),250))
    cpt = 1
    for l in lol_accession:
        print("/".join([str(cpt),str(len(lol_accession))])+" (250 items)")
        cpt += 1
        json_data = bioDBnet_request(l, headers, params)
        for i in json_data:
            for key,data in i.items():
                uniprot = i['InputValue']
                genbank = i['GenBank Protein Accession'].split('//')
                dico_data[uniprot] = genbank
    return dico_data

def get_biodbnet_data(l,db):
    """
    l : liste of NCBI ref Seq Accessions
    db : list of output db
    """
    dico_data = dict()
    
    headers = {
        'Connection': 'keep-alive',
        'Pragma': 'no-cache',
        'Cache-Control': 'no-cache',
        'Upgrade-Insecure-Requests': '1',
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.110 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8',
        'Accept-Encoding': 'gzip, deflate, br',
        'Accept-Language': 'fr-FR,fr;q=0.9,en-US;q=0.8,en;q=0.7',
    }

    params = {
        'method': 'db2db',
        'input': 'RefSeq Protein Accession',
        'inputValues': '',
        'outputs': ','.join(db),
        'taxonId': '',
        'format': 'row',
    }


    lol_accession = list(chunks(list(l),250))
    cpt = 1
    for l in lol_accession:
        print("/".join([str(cpt),str(len(lol_accession))])+" (250 items)")
        cpt += 1
        json_data = bioDBnet_request(l, headers, params)
        for i in json_data:
            for key,data in i.items():
                if key not in dico_data:
                    dico_data[key] = [data]
                else:
                    dico_data[key] += [data]
    return dico_data

def bioDBnet_database(id_list):
    """
    Return the user choosen DB
    Defaut, return GO-TERM DB
    """
    data = {
        1 : 'Affy ID',
        2 : 'Agilent ID',
        3 : 'Allergome Code',
        4 : 'ApiDB_CryptoDB ID',
        5 : 'Biocarta Pathway Name',
        6 : 'BioCyc ID',
        7 : 'CCDS ID',
        8 : 'Chromosomal Location',
        9 : 'CleanEx ID',
        10 : 'CodeLink ID',
        11 : 'COSMIC ID',
        12 : 'CPDB Protein Interactor',
        13 : 'CTD Disease Info',
        14 : 'CTD Disease Name',
        15 : 'CYGD ID',
        16 : 'dbSNP ID',
        17 : 'dictyBase ID',
        18 : 'DIP ID',
        19 : 'DisProt ID',
        20 : 'DrugBank Drug ID',
        21 : 'DrugBank Drug Info',
        22 : 'DrugBank Drug Name',
        23 : 'EC Number',
        24 : 'EchoBASE ID',
        25 : 'EcoGene ID',
        26 : 'Ensembl Biotype',
        27 : 'Ensembl Gene ID',
        28 : 'Ensembl Gene Info',
        29 : 'Ensembl Protein ID',
        30 : 'Ensembl Transcript ID',
        31 : 'FlyBase Gene ID',
        32 : 'FlyBase Protein ID',
        33 : 'FlyBase Transcript ID',
        34 : 'GAD Disease Info',
        35 : 'GAD Disease Name',
        36 : 'GenBank Nucleotide Accession',
        37 : 'GenBank Nucleotide GI',
        38 : 'GenBank Protein Accession',
        39 : 'GenBank Protein GI',
        40 : 'Gene ID',
        41 : 'Gene Info',
        42 : 'Gene Symbol',
        43 : 'Gene Symbol and Synonyms',
        44 : 'Gene Symbol ORF',
        45 : 'Gene Synonyms',
        46 : 'GeneFarm ID',
        47 : 'GO - Biological Process',
        48 : 'GO - Cellular Component',
        49 : 'GO - Molecular Function',
        50 : 'GO ID',
        51 : 'GSEA Standard Name',
        52 : 'H-Inv Locus ID',
        53 : 'HAMAP ID',
        54 : 'HGNC ID',
        55 : 'HMDB Metabolite',
        56 : 'Homolog - All Ens Gene ID',
        57 : 'Homolog - All Ens Protein ID',
        58 : 'Homolog - All Gene ID',
        59 : 'Homolog - Human Ens Gene ID',
        60 : 'Homolog - Human Ens Protein ID',
        61 : 'Homolog - Human Gene ID',
        62 : 'Homolog - Mouse Ens Gene ID',
        63 : 'Homolog - Mouse Ens Protein ID',
        64 : 'Homolog - Mouse Gene ID',
        65 : 'Homolog - Rat Ens Gene ID',
        66 : 'Homolog - Rat Ens Protein ID',
        67 : 'Homolog - Rat Gene ID',
        68 : 'HomoloGene ID',
        69 : 'HPA ID',
        70 : 'HPRD ID',
        71 : 'HPRD Protein Complex',
        72 : 'HPRD Protein Interactor',
        73 : 'Illumina ID',
        74 : 'IMGT/GENE-DB ID',
        75 : 'InterPro ID',
        76 : 'IPI ID',
        77 : 'KEGG Disease ID',
        78 : 'KEGG Gene ID',
        79 : 'KEGG Orthology ID',
        80 : 'KEGG Pathway ID',
        81 : 'KEGG Pathway Info',
        82 : 'KEGG Pathway Title',
        83 : 'LegioList ID',
        84 : 'Leproma ID',
        85 : 'Locus Tag',
        86 : 'MaizeGDB ID',
        87 : 'MEROPS ID',
        88 : 'MGC(ZGC/XGC) ID',
        89 : 'MGC(ZGC/XGC) Image ID',
        90 : 'MGC(ZGC/XGC) Info',
        91 : 'MGI ID',
        92 : 'MIM ID',
        93 : 'MIM Info',
        94 : 'miRBase ID',
        95 : 'NCIPID Pathway Name',
        96 : 'NCIPID Protein Complex',
        97 : 'NCIPID Protein Interactor',
        98 : 'NCIPID PTM',
        99 : 'Orphanet ID',
        100 : 'PANTHER ID',
        101 : 'Paralog - Ens Gene ID',
        102 : 'PBR ID',
        103 : 'PDB ID',
        104 : 'PeroxiBase ID',
        105 : 'Pfam ID',
        106 : 'PharmGKB Drug Info',
        107 : 'PharmGKB Gene ID',
        108 : 'PIR ID',
        109 : 'PIRSF ID',
        110 : 'PptaseDB ID',
        111 : 'PRINTS ID',
        112 : 'ProDom ID',
        113 : 'PROSITE ID',
        114 : 'PseudoCAP ID',
        115 : 'PubMed ID',
        116 : 'Reactome ID',
        117 : 'Reactome Pathway Name',
        118 : 'REBASE ID',
        119 : 'RefSeq Genomic Accession',
        120 : 'RefSeq Genomic GI',
        121 : 'RefSeq mRNA Accession',
        122 : 'RefSeq ncRNA Accession',
        123 : 'RefSeq Nucleotide GI',
        124 : 'RefSeq Protein Accession',
        125 : 'RefSeq Protein GI',
        126 : 'Rfam ID',
        127 : 'RGD ID',
        128 : 'SGD ID',
        129 : 'SMART ID',
        130 : 'STRING Protein Interactor',
        131 : 'TAIR ID',
        132 : 'Taxon ID',
        133 : 'TCDB ID',
        134 : 'TIGRFAMs ID',
        135 : 'TubercuList ID',
        136 : 'UCSC ID',
        137 : 'UniGene ID',
        138 : 'UniProt Accession',
        139 : 'UniProt Entry Name',
        140 : 'UniProt Info',
        141 : 'UniProt Protein Name',
        142 : 'UniSTS ID',
        143 : 'VectorBase Gene ID',
        144 : 'VEGA Gene ID',
        145 : 'VEGA Protein ID',
        146 : 'VEGA Transcript ID',
        147 : 'WormBase Gene ID',
        148 : 'WormPep Protein ID',
        149 : 'XenBase Gene ID',
        150 : 'ZFIN ID',
    }
    if id_list:
        return [data[i] for i in id_list]
    else:
        return ['UniProt Accession', 'GO - Biological Process', 'GO - Cellular Component', 'GO - Molecular Function']
    
def bioDBnet_request(l, headers, params, taxid=False):
    """
    Make the request to the BioDBnet website

    input :
        l (liste) : list of the value to be converted
        headers (dict) : parameter of the header of the URL request
        params (dict) : parameter of the URL request
        taxid (string) default = False : String of the taxid , default False
    output :
        json_data (dict) : JSON response as a dict of the response request
    """
    s = ",".join(l)
    params['inputValues'] = s
    if taxid:
        params['taxonId'] = taxid
    r = requests.get('https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json', headers=headers, params=params)
    json_data = json.loads(r.text)
    return json_data

    
def download_pfam(PFAMS):
    """
    Download the MSA of the UNIPROT & NCBI sequence from the PFAM familly

    input:
        PFAMS (list) : List of PFAM accession number
    output:
        PFAMS_PATH (list) : List of file with NCBI & UNIPROT sequence per accession of the PFAM in FASTA
    """
    #Cretae the PFAM Folder
    PFAMS_PATH = list()
    
    if os.path.isdir(os.path.join(os.getcwd(),"PFAMS")) == False:
        os.mkdir(os.path.join(os.getcwd(),"PFAMS"))

    #in case the list are empty
    if PFAMS:

        #request parameter
        pfam_params = {'format':'fasta',
                  'alnType':'ncbi',
                  'order':'a',
                  'case':'u',
                  'gaps':'dashes',
                  'download':'0'}

        for pfam in PFAMS:

            #Create folder for the specific PFAM
            path = os.path.join(os.getcwd(),"PFAMS",pfam)
            if os.path.isdir(path) == False:
                os.mkdir(path)

            #NCBI
            print('Download PFAM MSA alignment with NCBI sequences')
            filename_NCBI = os.path.join(path,"PFAM_{}_NCBI.fasta".format(pfam))
            if not os.path.isfile(filename_NCBI):
                resp = requests.get("http://pfam.xfam.org/family/{}/alignment/ncbi/format?".format(pfam), params=pfam_params)
                with open(filename_NCBI, "w") as f:
                    f.write(resp.content)

            #UNIPROT
            print('Download PFAM MSA alignment with UNIPROT sequences')
            filename_UNIPROT = os.path.join(path,"PFAM_{}_UNIPROT.fasta".format(pfam))
            if not os.path.isfile(filename_UNIPROT):
                resp = requests.get("http://pfam.xfam.org/family/{}/alignment/uniprot/format?".format(pfam), params=pfam_params)
                with open(filename_UNIPROT, "w") as f:
                    f.write(resp.content)
                    
            #READ THE ALIGNMENT TO EXTRACT ACCESSIONS (NCBI)
            accession_PFAM_NCBI = set([record.id.split('/')[0] for record in SeqIO.parse(filename_NCBI,'fasta')])
            accession_PFAM_UNIPROT = set([record.id.split('/')[0] for record in SeqIO.parse(filename_UNIPROT,'fasta')])
            #convert the UNIPROT to NCBI accession
            print('Convert the UNIPROT accession from the PFAM UNIPROT MSA to NCBI accessions')
            dico_uniprot_to_genbank = uniprot2ncbi_accession(accession_PFAM_UNIPROT)

            #check if UNIPROT = NCBI accesion inside the PFAM of NCBI
            print("check the UNIPROT-NCBI converted accession into the NCBI PFAM")
            uniprot_not_in_ncbi = list()
            for uniprot in dico_uniprot_to_genbank:
                list_genbank = dico_uniprot_to_genbank[uniprot]
                inside = False
                for gbk in list_genbank:
                    if gbk in accession_PFAM_NCBI:
                        inside = True
                if not inside:
                    uniprot_not_in_ncbi += [uniprot]
                    print("Adding {} to the actual NCBI Accession".format(uniprot))
            ### GET THE SAME FILE AS FASTA ###
            pfam_no_seq = os.path.join(path,'PFAM_{}_NOGAP_SEQ.fasta'.format(pfam))
            PFAMS_PATH += [pfam_no_seq]
            with open(pfam_no_seq, "w") as f:
                for records in SeqIO.parse(filename_NCBI,'fasta'):
                    iD = records.id
                    seq = str(records.seq).replace('-','')
                    f.write('>'+iD+'\n')
                    f.write(seq+'\n')
                    
                for records in SeqIO.parse(filename_UNIPROT,'fasta'):
                    iD = records.id
                    iD_split = records.id.split('/')[0]
                    #case when common accession are in the UNIPROT and NCBI
                    if iD_split not in accession_PFAM_NCBI:
                        if iD_split in uniprot_not_in_ncbi:
                            seq = str(records.seq).replace('-','')
                            f.write('>'+iD+'\n')
                            f.write(seq+'\n')
                            

            print('Create BLAST DB')
            make_blastdb([os.path.join(path,'PFAM_{}_NOGAP_SEQ.fasta'.format(pfam))])
            
    return PFAMS_PATH

def acc_to_ftp_path(acc):
    """
    from : https://github.com/ctSkennerton/scriptShed/blob/master/download_ncbi_assembly.py

    Return the FTP Path given the GCA accession number of a genome assembly
    
    input:
        acc (string) : GCA number of the NCBI Genome Assembly
    output:
        prefix (string)  : prefix of the FTP path
        accp (string)    : path of the GCA accession number
        version (string) : version of the assembly
    """
    match = re.search('(\w+)_(\d+)\.(\d)', acc)
    if match:
        prefix = match.group(1)
        accn = match.group(2)
        version = match.group(3)
        accp = re.findall('.{3}', accn)
        accp = '/'.join(accp)
        return (prefix, accp, version)
    else:
        raise ValueError("could not get FTP path from ".format(acc))


def get_files_from_types(types, base_name, ftp, path, ftp_true=False):
    """
    from https://github.com/ctSkennerton/scriptShed/blob/master/download_ncbi_assembly.py

    Download the genome
    
    input :
        types (list)      : List of NCBI type suffix file
        bas_name (string) : Basename of the genome file
        ftp (ftp object)  : ftp session object to the NCBI FTP Server
        path (string)     : output folder
    return:
        out : path of the downloaded genome
    """
    
    out = str()
    if ftp_true:

        #download
        for t in types:
            f = path+base_name+types[0]
            try:
                urllib.urlretrieve(ftp_true+'/'+base_name+t, f)
            except IOError:
                print('No proteomic file availlable for : {}'.format(base_name))
                return False
                    
            urllib.urlcleanup()
                
            #extract
            with gzip.open(f, 'rb') as f_in:
                out = f.replace('.gz','')
                with open(out, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    
            #remove file
            os.remove(f)
            return out
    
    else:
        
        for g in ftp.nlst():
            for t in types:
                
                if g == base_name + t:

                    #download
                    try:
                        urllib.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/{}/{}".format(ftp.pwd(), g), path+g)
                    except IOError:
                        print('No proteomic file availlable for : {}'.format(g))
                        return False
                    urllib.urlcleanup()
                        
                    #extract
                    with gzip.open(path+g, 'rb') as f_in:
                        out = path+g.replace('.gz','')
                        with open(out, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                            
                    #remove file
                    os.remove(path+g)
        return out

def map_types_to_file_suffix(types):
    """
    from https://github.com/ctSkennerton/scriptShed/blob/master/download_ncbi_assembly.py
    Return a list of the NCBI type suffix file
    
    input:
        type (list) : List of type
    output
        List of NCBI type suffix file
    """
    mapping = {'fna': '_genomic.fna.gz',
               'faa': '_protein.faa.gz',
               'ffn': '_cds_from_genomic.fna.gz',
               'gb': '_genomic.gbff.gz',
               'gff': '_genomic.gff.gz'
              }
    return [mapping[i] for i in types]


def download_genome(taxid, df_refseq, df_euk, df_prok):
    """
    Download the Refseq Genome  Taxid
    If no refseq Download all the proteome for each assembly
    
    input:
        taxid (list)          : list of taxid
        df_refseq (DataFrame) : pandas DataFrame of the RefSeq Gneome
        df_euk (DataFrame)    : pandas DataFrame of the eukaryotes genomes
        df_prok (DataFrame)   : pandas DataFrame of the prokaryotes genomes
    output:
        Path of the proteome for each assembly of each taxid in fasta
    """
    #get the taxid in the refseq genomes
    set1 = set(taxid)
    #get the reference genome by the given taxid
    df_selection  = df_refseq.loc[df_refseq['taxid'].isin(taxid)]
    #convert the taxid into a set (can be have more than 1 assembly in ref seq for 1 taxid)
    set2 = set(df_selection['taxid'].values.tolist())
    #check if all the given taxid are in the refseq genome taxid
    diff = set1.difference(set2)

    PATH_LIST = list()
    
    #check if the dataframe is empty
    if not df_selection.empty :
        #not empty (but can not contain all taxid)
        #get the FTP path and download the genome
        list_dict = df_selection.to_dict(orient='records')
        taxid_list = df_selection['taxid'].values.tolist()
        for d in list_dict:
            #make link
            path = d['ftp_path']
            root = path.split('/')[-1]
            f1 = path+"/"+root+'_protein.faa.gz'
            org = d['organism_name']
            strain = d['# assembly_accession']
            print('Downloading {} - {} ...'.format(org,strain))
            
            f2 = os.path.join(os.getcwd(),org.upper().replace(' ','_'),root+'_protein.faa.gz')
            PATH_LIST += [f2.replace('.gz','')]
            #create directory
            if os.path.isdir(os.path.join(os.getcwd(),org.upper().replace(' ','_'))) == False:
                os.mkdir(os.path.join(os.getcwd(),org.upper().replace(' ','_')))
            #download
            urllib.urlretrieve(f1, f2)
            #extract
            with gzip.open(f2, 'rb') as f_in:
                with open(f2.replace('.gz',''), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            #remove file
            os.remove(f2)
            
    #if df_selection dont have all taxid
    if diff:
        GCA_list = list()
        org_list = list()
        ftp_list = list()
        df_selection_euk  = df_euk.loc[df_euk['TaxID'].isin(diff)]
        #eukaryotes.txt dont have same header as the prokaryote.txt, prokaryote have 'Strain', 'Pubmed ID', 'FTP Path', 'Reference' in more
        if not df_selection_euk.empty :
            GCA_list += df_selection_euk['Assembly Accession'].values.tolist()
            org_list += df_selection_euk['#Organism/Name'].values.tolist()
            ftp_list += [False] * len(df_selection_euk.index)
        df_selection_prok = df_prok.loc[df_prok['TaxID'].isin(diff)]
        if not df_selection_prok.empty :
            GCA_list += df_selection_prok['Assembly Accession'].values.tolist()
            org_list += df_selection_prok['#Organism/Name'].values.tolist()
            ftp_list += df_selection_prok['FTP Path'].values.tolist()
        GCA_ORG_list = zip(GCA_list, org_list, ftp_list)

        if GCA_ORG_list:
            #from : https://github.com/ctSkennerton/scriptShed/blob/master/download_ncbi_assembly.py
            #Connection to the FTP server of NCBI - get only faa protein file
            types = map_types_to_file_suffix(['faa'])
            ftp = ftplib.FTP()
            ftp.connect('ftp.ncbi.nlm.nih.gov')
            ftp.login(user='anonymous',passwd='your.email@gmail.com')
            
            for accession, organism, ftp_true in GCA_ORG_list:
                #create folder
                path = os.path.join(os.getcwd(),organism.upper().replace(' ','_'),'')
                if os.path.isdir(path) == False:
                    os.mkdir(path)
                    
                prefix, accp, version = acc_to_ftp_path(accession)
                base_path = "genomes/all/{}/{}".format(prefix,accp)
                ftp.cwd(base_path)
                for f in ftp.nlst():
                    if accession in f:
                        ftp.cwd(f)
                        out = get_files_from_types(types, f, ftp, path, ftp_true=ftp_true)
                        if out:
                            PATH_LIST += [out]
                ftp.cwd("/")
        else:
            print('No proteome found for the choosen node or taxid, abort')
    return PATH_LIST

def make_blastdb(PATH_LIST):
    """
    Make the BLAST Database for each PFAM fasta file
    """
    for i in PATH_LIST:
        if FLAG:
            CREATE_NO_WINDOW = 0x08000000
            print("Make BLAST DB for {}".format(os.path.basename(i)))
            a = subprocess.call('makeblastdb -in {} -dbtype prot'.format(i), creationflags=CREATE_NO_WINDOW)
        else:
            print("Make BLAST DB for {}".format(os.path.basename(i)))
            a = subprocess.call('makeblastdb -in {} -dbtype prot'.format(i), shell=True)            

def split_fasta(fasta):
    """
    input : fasta file complete path
    output : create individual file with the fasta sequence 
    """
    FASTA_PATH = list()
    path, fasta_filename = os.path.split(fasta)
    folder_name = fasta_filename.split('.')[0]
    if os.path.isdir(os.path.join(path,folder_name+'_fasta')) == False:
        os.mkdir(os.path.join(path,folder_name+'_fasta'))
        
    for records in SeqIO.parse(fasta,'fasta'):
        iD = records.id
        seq = str(records.seq)
        fasta = os.path.join(path,folder_name+'_fasta',iD.replace('.','_')+'.fasta')
        FASTA_PATH += [fasta]
        with open(fasta,"w") as f:
            f.write('>'+records.description+'\n')
            f.write(seq+'\n')
    return FASTA_PATH

def blast(FASTA, PFAM, OUT):
    """
    Make a blast with a evalue of 10000
    
    FASTA = the path to FASTA file
    PFAM = The BLAST DB PFAM Database (protein)
    OUT = Output blast repport in XML format (outfmt = 5)
    """
    size = False
    cpt = 0
    while size == False:
        blastp_cline = NcbiblastpCommandline(query=FASTA, db=PFAM, evalue=10000,outfmt=5, out=OUT)
        cmd = str(blastp_cline)
        if FLAG:
            CREATE_NO_WINDOW = 0x08000000
            a = subprocess.call(cmd, creationflags=CREATE_NO_WINDOW)
        else:
            a = subprocess.call(cmd, shell=True)

        #check the size
        #if the file size is not correct the while loop will continue
        if os.stat(OUT).st_size != 0:
            break
        cpt +=1
        if cpt == 100:
            print('Trying 100 iteration, not working, exit for '+FATSA)
            

    
###########
# RUNNING #
###########

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''
    Welcome to PyFuncover !

    Python Function Uncover ( PyFuncover ) is a new bioinformatic tool
    able to search for protein with a specific function in a full proteome.
    The pipeline coded in python uses BLAST alignment and the sequences from
    a PFAM family as search seed.

    === REQUIREMENT ===

    BLAST
    Python dependancies : BioPython, Numpy, Matplotlib, Pandas

    USAGE :
        python PyFuncover.py -taxid [TAXID ...] -pfam [PFAM ...] --db [DB ...] --out [OUT]
        python PyFuncover.py --update

    TO UPDATE THE DATABASES :

        python PyFuncover.py --update

    === OBLIGATORY ARGUMENTS : ===

    -pfam : List of PFAM familly ID : PF####
            each separated by a blank space

        PyFuncover.py -pfam PF14651 PF#### ...
        

    -taxid: The list of TaxID for each organism you want to download a proteome
        Each separated by a space (for example Human and Yeast (S299C) taxid)

            PyFuncover.py -taxid 9606 559292

        Can be a Taxid that represent a node in the phylogenetic tree
        (Eukaryotes : 2759 ; Insecta : 50557, ...)
        He will retrieve all availlable assembly for them

    === OPTIONAL : ===

    --update :
        Download the last release of the NCBI Taxonomic Database
        Download the last RefSeq, Prokaryote and Eukaryote Genome Assembly List

        PyFuncover.py --update

    --out : Filename output
        Format are in CSV format (pandas.to_csv output)
        Default : result.csv

    --db : The list of choosen database number to retrieve data from bioDBnet :
            default : 137 45 46 47 (UNIPROT ID, GO-TERMs Databases)

            1 : Affy ID
            2 : Agilent ID
            3 : Allergome Code
            4 : ApiDB_CryptoDB ID
            5 : Biocarta Pathway Name
            6 : BioCyc ID
            7 : CCDS ID
            8 : Chromosomal Location
            9 : CleanEx ID
            10 : CodeLink ID
            11 : COSMIC ID
            12 : CPDB Protein Interactor
            13 : CTD Disease Info
            14 : CTD Disease Name
            15 : CYGD ID
            16 : dbSNP ID
            17 : dictyBase ID
            18 : DIP ID
            19 : DisProt ID
            20 : DrugBank Drug ID
            21 : DrugBank Drug Info
            22 : DrugBank Drug Name
            23 : EC Number
            24 : EchoBASE ID
            25 : EcoGene ID
            26 : Ensembl Biotype
            27 : Ensembl Gene ID
            28 : Ensembl Gene Info
            29 : Ensembl Protein ID
            30 : Ensembl Transcript ID
            31 : FlyBase Gene ID
            32 : FlyBase Protein ID
            33 : FlyBase Transcript ID
            34 : GAD Disease Info
            35 : GAD Disease Name
            36 : GenBank Nucleotide Accession
            37 : GenBank Nucleotide GI
            38 : GenBank Protein Accession
            39 : GenBank Protein GI
            40 : Gene ID
            41 : Gene Info
            42 : Gene Symbol
            43 : Gene Symbol and Synonyms
            44 : Gene Symbol ORF
            45 : Gene Synonyms
            46 : GeneFarm ID
            47 : GO - Biological Process
            48 : GO - Cellular Component
            49 : GO - Molecular Function
            50 : GO ID
            51 : GSEA Standard Name
            52 : H-Inv Locus ID
            53 : HAMAP ID
            54 : HGNC ID
            55 : HMDB Metabolite
            56 : Homolog - All Ens Gene ID
            57 : Homolog - All Ens Protein ID
            58 : Homolog - All Gene ID
            59 : Homolog - Human Ens Gene ID
            60 : Homolog - Human Ens Protein ID
            61 : Homolog - Human Gene ID
            62 : Homolog - Mouse Ens Gene ID
            63 : Homolog - Mouse Ens Protein ID
            64 : Homolog - Mouse Gene ID
            65 : Homolog - Rat Ens Gene ID
            66 : Homolog - Rat Ens Protein ID
            67 : Homolog - Rat Gene ID
            68 : HomoloGene ID
            69 : HPA ID
            70 : HPRD ID
            71 : HPRD Protein Complex
            72 : HPRD Protein Interactor
            73 : Illumina ID
            74 : IMGT/GENE-DB ID
            75 : InterPro ID
            76 : IPI ID
            77 : KEGG Disease ID
            78 : KEGG Gene ID
            79 : KEGG Orthology ID
            80 : KEGG Pathway ID
            81 : KEGG Pathway Info
            82 : KEGG Pathway Title
            83 : LegioList ID
            84 : Leproma ID
            85 : Locus Tag
            86 : MaizeGDB ID
            87 : MEROPS ID
            88 : MGC(ZGC/XGC) ID
            89 : MGC(ZGC/XGC) Image ID
            90 : MGC(ZGC/XGC) Info
            91 : MGI ID
            92 : MIM ID
            93 : MIM Info
            94 : miRBase ID
            95 : NCIPID Pathway Name
            96 : NCIPID Protein Complex
            97 : NCIPID Protein Interactor
            98 : NCIPID PTM
            99 : Orphanet ID
            100 : PANTHER ID
            101 : Paralog - Ens Gene ID
            102 : PBR ID
            103 : PDB ID
            104 : PeroxiBase ID
            105 : Pfam ID
            106 : PharmGKB Drug Info
            107 : PharmGKB Gene ID
            108 : PIR ID
            109 : PIRSF ID
            110 : PptaseDB ID
            111 : PRINTS ID
            112 : ProDom ID
            113 : PROSITE ID
            114 : PseudoCAP ID
            115 : PubMed ID
            116 : Reactome ID
            117 : Reactome Pathway Name
            118 : REBASE ID
            119 : RefSeq Genomic Accession
            120 : RefSeq Genomic GI
            121 : RefSeq mRNA Accession
            122 : RefSeq ncRNA Accession
            123 : RefSeq Nucleotide GI
            124 : RefSeq Protein Accession
            125 : RefSeq Protein GI
            126 : Rfam ID
            127 : RGD ID
            128 : SGD ID
            129 : SMART ID
            130 : STRING Protein Interactor
            131 : TAIR ID
            132 : Taxon ID
            133 : TCDB ID
            134 : TIGRFAMs ID
            135 : TubercuList ID
            136 : UCSC ID
            137 : UniGene ID
            138 : UniProt Accession
            139 : UniProt Entry Name
            140 : UniProt Info
            141 : UniProt Protein Name
            142 : UniSTS ID
            143 : VectorBase Gene ID
            144 : VEGA Gene ID
            145 : VEGA Protein ID
            146 : VEGA Transcript ID
            147 : WormBase Gene ID
            148 : WormPep Protein ID
            149 : XenBase Gene ID
            150 : ZFIN ID

    --nb : The number of parrallelized BLAST process (default : 10)
    Be carefull, high number will use lot of memory and create a stck overflow !

    '''))
    parser.add_argument('-pfam',
                        nargs='*',
                        help='The list of PFAM accession (PF#####)')
    parser.add_argument('-taxid',
                        nargs='*',
                        help='The list of TaxID for each organism')
    parser.add_argument('--db',
                        nargs='*',
                        type=int,
                        help='bioDBnet Databases Number: (default : 45 46 47 (GO-TERMs Databases))')
    parser.add_argument('--out', nargs='*', help='File output, default "result.csv"')
    parser.add_argument('--update',help='Update the NCBI Taxonomic Database and the Prokaryote, Eukaryote & RefSeq genome assembly', action='store_true', default=False)
    parser.add_argument('--nb',help='Number of parrallelized BLAST processes', type=int)
    args = parser.parse_args()
    dicoArgs = vars(args)
    #print(dicoArgs)


    if not len(sys.argv) > 1:
        parser.print_help()
        print
        print('''Exit, no argument provided, minimal use :
python PyFuncover.py -taxid [TAXID] -pfam [PFAM]''')
        exit()

    if dicoArgs['update']:
        print(dicoArgs['update'])
        print("here")
        print('Updating Taxonomic NCBI Databases and the Prokaryote, Eukaryote & RefSeq genome assembly')
        init()
        exit()

    if dicoArgs['taxid'] is None:
        parser.print_help()
        print
        print('You miss -taxid argument')
        exit()
    else:
        taxid = dicoArgs['taxid']

    if dicoArgs['pfam'] is None:
        parser.print_help()
        print
        print('You miss -pfam argument')
        exit()
    else:
        pfam = dicoArgs['pfam']
        
    if dicoArgs['db'] is None:
        biodbnet_db = [47,48,49]
    else:
        biodbnet_db = dicoArgs['db']

    if dicoArgs['out'] is None:
        output = 'results.csv'
    else:
        output = dicoArgs['db']

    if dicoArgs['nb'] is None:
        NB_BLAST_PROCESS = 10
    else:
        NB_BLAST_PROCESS = dicoArgs['nb']

    if os.path.isdir(os.path.join(os.getcwd(),'TAXONOMY')) == False:
        print("Download the reference genome list")
        init()
    else:
        print("Taxonomic Database and genome list already downloaded")
        print("Use --update to update them")
        
    print("get Taxonomic nodes...")
    dicNode = dic_Node()
    print("get Taxonomic Name...")
    dicName, dicNameInverse = dic_Name()
    print("create Lineage...")
    dicLineage = lineage(dicName,dicNode)
    print("group Taxid...")
    dicGroupedTaxid = groupedTaxId(dicLineage, dicNameInverse)
    print("Read the reference genome assembly")
    df_refseq = pd.read_csv('assembly_summary_refseq.txt',delimiter='\t', skiprows=1, dtype=str)
    print("Read the eukaryotes genome assembly")
    df_euk = pd.read_csv('eukaryotes.txt',delimiter='\t', dtype=str)
    print("Read the prokaryote genome assembly")
    df_prok = pd.read_csv('prokaryotes.txt',delimiter='\t', dtype=str)

    print("Get the list of available taxid")
    set_taxid = set()
    for i in taxid:
        print('{} : {}'.format(i,dicName[i]) )
        set_taxid.update(dicGroupedTaxid[i])
        print('Availlable species :')
        for i in set_taxid:
            print('\n\t'+i+' : '+dicName[i])
              
        
    print('Download PFAM')
    PFAMS_PATH = download_pfam(pfam)
    if not PFAMS_PATH:
        print('PyFuncover : EXIT - No PFAMs available for the user selection')
        exit()
    print('Download Genome')
    PATH_LIST = download_genome(set_taxid, df_refseq, df_euk, df_prok)

    #check if they are any genome
    if not PATH_LIST:
        print('PyFuncover : EXIT - No proteome available for the user selection')
        exit()
            
    print('Split proteome into individual fasta')

    PATH_FASTA_GENOME = list()
    for fasta in PATH_LIST:
        PATH_FASTA = split_fasta(fasta)
        PATH_FASTA_GENOME += [PATH_FASTA]


    print('BLAST')
    for PFAM in PFAMS_PATH:
        PFAM_NAME = PFAM.split(os.sep)[-2]
        for PATH_FASTA in PATH_FASTA_GENOME:

            #get the path create the xml folder
            species_folder, fasta_assembly_folder, fasta_file = PATH_FASTA[0].split(os.sep)[-3:]
            root = os.getcwd()

            xml_assembly_folder = fasta_assembly_folder.replace('fasta','xml_'+PFAM_NAME)
            xml_folder = os.path.join(root, species_folder, xml_assembly_folder)
            if os.path.isdir(xml_folder) == False:
                os.mkdir(xml_folder)

            n = NB_BLAST_PROCESS
            PATH_FASTA_CHUNKED = list(chunks(PATH_FASTA, n))
            cpt = 1
            print('BLAST on {} - {}'.format(species_folder,fasta_assembly_folder))
            for FASTA_LIST in PATH_FASTA_CHUNKED:
                thr_list = list()
                print('{}/{} ({} items) : starting'.format(cpt,len(PATH_FASTA_CHUNKED),n))
                
                for FASTA in FASTA_LIST:
                    fasta_file = os.path.split(FASTA)[-1].replace('fasta','xml')
                    OUT = os.path.join(xml_folder,fasta_file)
                    #check if the file is not already there (from a previous analysis)
                    if os.path.isfile(OUT):
                        #check if the file is not empty (occur when memory error)
                        if os.stat(OUT).st_size == 0:
                            thr = Thread(target=blast, args=[FASTA,PFAM,OUT])
                            thr.start()
                            thr_list.append(thr)
                    #if the file doesnt exist
                    if not os.path.isfile(OUT):
                        thr = Thread(target=blast, args=[FASTA,PFAM,OUT])
                        thr.start()
                        thr_list.append(thr)
                        
                for thr in thr_list:
                    thr.join()
                cpt +=1

    print('PLOTTING')
    RESULTS_LIST = list()
    for PFAM in PFAMS_PATH:
        PFAM_NAME = PFAM.split(os.sep)[-2]
        for PATH_FASTA in PATH_FASTA_GENOME:
            #get the path folder of the fasta files
            species_folder, fasta_assembly_folder, fasta_file = PATH_FASTA[0].split(os.sep)[-3:]
            root = os.getcwd()
              
            #create the graph folder
            graph_assembly_folder = fasta_assembly_folder.replace('fasta','graph_'+PFAM_NAME)
            graph_folder = os.path.join(root, species_folder, graph_assembly_folder)
            if os.path.isdir(graph_folder) == False:
                os.mkdir(graph_folder)

            #create the csv folder
            csv_assembly_folder = fasta_assembly_folder.replace('fasta','csv_'+PFAM_NAME)
            csv_folder = os.path.join(root, species_folder, csv_assembly_folder)
            if os.path.isdir(csv_folder) == False:
                os.mkdir(csv_folder)

            #make the path for xml folder
            xml_assembly_folder = fasta_assembly_folder.replace('fasta','xml_'+PFAM_NAME)
            xml_folder = os.path.join(root, species_folder, xml_assembly_folder)

            #iterate on each fasta
            for FASTA in PATH_FASTA:

                #formate the filename
                basename = os.path.basename(FASTA)
                PNG = os.path.join(graph_folder, os.path.basename(FASTA).replace('fasta','png'))
                XML = os.path.join(xml_folder, os.path.basename(FASTA).replace('fasta','xml'))
                CSV = os.path.join(csv_folder, os.path.basename(FASTA).replace('fasta','csv'))
                #get the path of the FASTA by replacing xml
                protein = SeqIO.read(FASTA, "fasta")
                protein_seq = protein.seq
                protein_len = len(protein_seq)
                
                #parse the XML BLAST REPPORT
                blast_records = SearchIO.parse(XML, 'blast-xml')
                print('Parsing '+os.path.basename(XML)+'...')
                cpt = 0
                score_list = [0] * protein_len
                for record in blast_records:
                    for hsp in record.hsps:
                        cpt += 1
                        q_s = hsp.query_start
                        sim = hsp.aln_annotation['similarity'] #middle line of the blast hsp alignment
                        for i in range(len(hsp.query.seq)):
                            if hsp.query.seq[i] == '-':
                                q_s -= 1
                                continue
                            else:
                                pos = q_s + i
                                c = sim[i]
                                if c =="+":
                                    score_list[pos] = score_list[pos] + 1
                                if c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                                    score_list[pos] = score_list[pos] + 2
                
                #plt.figure(figsize=(20,14))
                #create colormap
                
                max_score = max(score_list)
                mean_score = np.mean(score_list)
                total_score = sum(score_list)
                normalised_score = total_score / protein_len
                a = protein.description

                if not os.path.isfile(PNG):
                    hexa_colors_dict = dict()
                    for i in range(0, max_score+1):
                        cm_object = cm.RdYlGn(1.*i/max_score)
                        rgb = cm_object[:3]
                        hexa = colors.rgb2hex(rgb)
                        hexa_colors_dict[i] = hexa
                        
                    #attribute color
                    color = list()
                    for i in score_list:
                        color += [hexa_colors_dict[i]]
                    fig, ax = plt.subplots()
                    ax.yaxis.grid()
                    bars = plt.bar(range(1,len(score_list)+1,1),score_list, color=color)
                    
                    #split protein description to have nice title
                    plt.title(a, wrap = True)
                    
                    if not os.path.isfile(CSV):
                        df = pd.DataFrame(data=score_list, index=list(protein_seq), columns=['score'])
                        df.to_csv(CSV)
                        plt.savefig(PNG)
                        plt.close()

                    plt.close()
                
                #STORE RESULT
                #ACCESSION | PROTEIN NAME | SPECIES | TAXID | GCA | MAX SCORE | MEAN SCORE | TOTAL SCORE | 
                NCBI_ACCESSION = '.'.join(basename.rsplit('_',1)).replace('.fasta','')
                SPECIE = re.findall('\[(.*?)\]',protein.description)[-1]
                TAXID = dicNameInverse[SPECIE]
                GCA = "_".join(fasta_assembly_folder.split('_')[:2])
                print(NCBI_ACCESSION, a, SPECIE, TAXID, GCA, max_score, mean_score, total_score, normalised_score)
                RESULTS_LIST += [[NCBI_ACCESSION, protein.description, SPECIE, TAXID, GCA, max_score, mean_score, total_score]]

    print("Make request to BioDBnet")
    columns = ['ACCESSION','PROTEIN DESCRIPTION','SPECIES','TAXID','GCA','MAX SCORE','MEAN SCORE','TOTAL SCORE']
    df_result = pd.DataFrame(RESULTS_LIST,columns=columns)
    #formate accession without the dot
    list_of_accession = [i.split('.')[0] for i in df_result['ACCESSION']]
    #get the DB to query
    db_list = bioDBnet_database(biodbnet_db)
    data = get_biodbnet_data(list_of_accession,db_list)

    for column in data:
        dico = { df_result['ACCESSION'][row_index] : data[column][row_index] for row_index in range(len(data[column])) }
        df_result[column]= df_result['ACCESSION'].map(dico)
            
    df_result.to_csv(output, index=None)
    

if __name__ == '__main__':
    main()
