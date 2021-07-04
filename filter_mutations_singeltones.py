#!/usr/bin/env python
# coding: utf-8

# In[186]:



from pymongo import MongoClient
import pandas as pd

client = MongoClient("mongodb://localhost:27017/")
database = client["Research"]
col_var = database["Variants"]


chrom_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"]

singleton_all = []



for chr in chrom_list: 
    
    pipeline= [
     {
        u"$match": {
            u"AC": 1.0,
            u"FILTER": u"PASS",
            u"CHROM": chr
        }
    }, 
    {
        u"$project": {
            u"AN": 1.0,
            u"vep_dictionary": 1.0,
            u"FILTER": 1.0,
            u"POS": 1.0,
            u"AN_male": 1.0,
            u"AN_female": 1.0,
            u"AC": 1.0,
            u"AC_male": 1.0,
            u"AC_female": 1.0
        }
    }, 
    {
        u"$lookup": {
            u"from": u"Coverage",
            u"localField": u"POS",
            u"foreignField": u"pos",
            u"as": u"med_coverage"
        }
    }, 
    {
        u"$unwind": {
            u"path": u"$med_coverage"
        }
    }, 
    {
        u"$match": {
            u"med_coverage.median": {
                u"$gte": 1.0
            },
            u"med_coverage.chrom": chr
        }
    }, 
    {
        u"$group": {
            u"_id": {
                u"Gene": u"$vep_dictionary.Feature",
                u"Consequence": u"$vep_dictionary.Consequence"
            },
            u"sum_AN_male": {
                u"$sum": u"$AN_male"
            },
            u"sum_AN_female": {
                u"$sum": u"$AN_female"
            },
            u"sum_AC_female": {
                u"$sum": u"$AC_female"
            },
            u"sum_AC_male": {
                u"$sum": u"$AC_male"
            },
            u"sum_AC": {
                u"$sum": u"$AC"
            },
            u"count": {
                u"$sum": 1.0
            }
        }
    }
    ]
    
    result = col_var.aggregate(pipeline)
    #list(result)
    result_df= pd.DataFrame(data = result)
    result_df['CHROM'] = str(chr)
    singleton_all.append(result_df)




all_chrom = pd.concat(singleton_all)



all_chrom= all_chrom[all_chrom._id != {}]



gene= []
Consequence= []
for i in range(len(all_chrom.iloc[:,0])):
    gene.append(all_chrom.iloc[i,0]['Gene'])
    Consequence.append(all_chrom.iloc[i,0]['Consequence'])

all_chrom['Gene']= gene
all_chrom['Consequence']= Consequence


# In[194]:


cols = all_chrom.columns.tolist()
cols = cols[-1:] + cols[:-1]
cols = cols[-1:] + cols[:-1]

all_chrom = all_chrom[cols]


all_chrom.head()



del all_chrom['_id']


import os
os.chdir("/home/ethel/Documents/all chrom")
all_chrom.to_csv("all_chr_singelton.csv")



