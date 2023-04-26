
import pandas as pd
import mudata as md
from copy import deepcopy
import numpy as np

from liana.testing._sample_anndata import generate_toy_mdata
from liana.method.sp._misty import misty

data_dir = "/Users/pschafer/Python_projects/misty_liana/data"

feature_map = pd.read_csv(f"{data_dir}/cytosig_map_mouse.csv", index_col=0)
interaction_tuples = list(zip(feature_map['gene'], feature_map['signature']))
predictors, targets = zip(*interaction_tuples)

mdata = md.read_h5mu(f"{data_dir}/f4hr.h5ad")
mdata = mdata[mdata.obs["rna:leiden"] != "5"]


### start: check that count matrices are not modified ###
mdata_rna = deepcopy(mdata["rna"])
mdata_cytosig = deepcopy(mdata["cytosig"])

misty(mdata, x_mod="rna", y_mod="cytosig", predictors=predictors, targets=targets, bandwidth=500, kernel="gaussian", bypass_intra=False, add_juxta=False, add_self=False)

assert np.all(mdata_rna.X.todense() ==  mdata["rna"].X.todense())
assert np.all(mdata_cytosig.X.todense() == mdata["cytosig"].X.todense())

misty(mdata, x_mod="rna", y_mod="cytosig", predictors=predictors, targets=targets, bandwidth=500, kernel="gaussian", bypass_intra=False, add_juxta=False, add_self=False, group_intra_by="leiden", group_env_by="leiden")

assert np.all(mdata_rna.X.todense() ==  mdata["rna"].X.todense())
assert np.all(mdata_cytosig.X.todense() == mdata["cytosig"].X.todense())
### end: check that count matrices are not modified ###


