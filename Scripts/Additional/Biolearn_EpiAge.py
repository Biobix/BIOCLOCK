#!/usr/bin/env python3.10

import os
import sys
import feather
import pyarrow.feather as feather
import pandas as pd
from biolearn.model_gallery import ModelGallery

# Access command-line arguments passed from R
pipeline = sys.argv[1] 

os.chdir("/data/user_homes/mennovd/BIOKLOK/Results/PrePro/" + pipeline)

betas = feather.read_feather("betas_no_suffix.feather")
betas.set_index("id", inplace=True)

pheno = feather.read_feather("mAge.feather")
pheno["sex"] = pheno["sex"].astype(int)
pheno["id"] = ["X" + value for value in pheno["id"]]
pheno.set_index("id", inplace=True)

from biolearn.data_library import GeoData
data = GeoData(metadata = pheno, dnam = betas)

from biolearn.model_gallery import ModelGallery
gallery = ModelGallery()

predictors = ["Horvathv1", "Horvathv2", "GrimAgeV1", "GrimAgeV2"]

for predictor in predictors:
    print(predictor)
    if predictor in {"GrimAgeV1", "GrimAgeV2"}:
        res = ModelGallery().get(predictor).predict(data)
        res = res.drop(columns = ["Age", "Female", "AgeAccelGrim"])
    else:
        res = ModelGallery().get(predictor).predict(data)
        res = res["Predicted"]
        res.name = predictor
    if predictor == predictors[0]:
        ages = res.copy()
    else:
        ages = pd.concat([ages, res], axis = 1)
        
feather.write_feather(ages, "biolearn_mAge.feather")
