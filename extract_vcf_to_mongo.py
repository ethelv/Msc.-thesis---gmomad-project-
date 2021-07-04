import os
import vcfpy
import pymongo
import itertools

start_index = 0

vep_names = ["Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE",
             "EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position",
             "Amino_acids","Codons","Existing_variation","ALLELE_NUM","DISTANCE","STRAND","FLAGS",
             "VARIANT_CLASS","MINIMISED","SYMBOL_SOURCE","HGNC_ID","CANONICAL","TSL","APPRIS","CCDS",
             "ENSP","SWISSPROT","TREMBL","UNIPARC","GENE_PHENO","SIFT","PolyPhen","DOMAINS","HGVS_OFFSET",
             "GMAF","AFR_MAF","AMR_MAF","EAS_MAF","EUR_MAF","SAS_MAF","AA_MAF","EA_MAF","ExAAF",
             "ExAC_AFR_MAF","ExAC_AMR_MAF","ExAC_EAS_MAF","ExAC_FIN_MAF","ExAC_NFE_MAF","ExAC_OTH_MAF",
             "ExAC_SAS_MAF","CLIN_SIG","SOMATIC","PHENO","PUBMED","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS",
             "MOTIF_SCORE_CHANGE","LoF","LoF_filter","LoF_flags","LoF_info"]


def save_to_db(chunk):
    mycol.insert_many(chunk)

def normalize_vep(vep):
    vep_dictionary = {}
    vep_list = vep.split("|")
    for i in range(len(vep_names)):
        vep_dictionary[vep_names[i]] = vep_list[i]
    return vep_dictionary


CHUNK_SIZE = 5000

# Set working directory
os.chdir("/media/ethel/disk/Data")

# define file name
reader = vcfpy.Reader.from_path('gnomad.exomes.r2.1.1.sites.vcf')

# Connect to DB
myclient = pymongo.MongoClient("mongodb://192.168.198.132:27017/")
mydb = myclient["Research"]
mycol = mydb["VariantsV2"]


# counter = 0
done = 0
chunk = []
row_index = 0

# Jump to last checkpoint start point
reader = itertools.islice(reader, start_index+1, None)
for row in reader:
    if row_index < start_index:
        row_index+=1
        if row_index % 1000 == 0:
            print(str(row_index)+"/"+str(start_index))
    
    for vep in row.INFO['vep']:
        item = {}
        if type(row.ID) == list and len(row.ID)==1 :
            ID= row.ID[0]
        else:
            ID= None

        
        item["CHROM"]       = row.CHROM
        item["POS"]         = row.POS
        item["ID"]          = ID
        item["REF"]         = row.REF
        item["ALT_TYPE"]    = row.ALT[0].type
        item["ALT_VALUE"]   = row.ALT[0].value
        item["QUAL"]        = row.QUAL
        item["FILTER"]      = row.FILTER[0]
        item['vep_dictionary'] = normalize_vep(vep)
        item['vep'] = vep
        for key in row.INFO:
            if key == 'vep':
                continue
            elif type(row.INFO[key]) == list and len(row.INFO[key])==1 :
                item[key] = row.INFO[key][0]
            else:
                item[key] = row.INFO[key]
        chunk.append(item)
    row_index += 1


    if len(chunk) >= CHUNK_SIZE:
        save_to_db(chunk)
        done += len(chunk)
        print('done')
        print(done)
        print('Row index', row_index)
        chunk = []

    # counter = counter+1
    # if counter == 5001:
    # 	break

if len(chunk) > 0:
    save_to_db(chunk)
    chunk = []      
    
      
