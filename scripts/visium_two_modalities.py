
import pandas as pd
import mudata as md
from pathlib import Path

from liana.method.sp._misty import misty

data_dir = Path(__file__).parent / ".." / "data"
print(data_dir)

feature_map = pd.read_csv(f"{data_dir}/cytosig_map_mouse.csv", index_col=0)
interaction_tuples = list(zip(feature_map['gene'], feature_map['signature']))
predictors, targets = zip(*interaction_tuples)

mdata = md.read_h5mu(f"{data_dir}/f4hr.h5ad")
mdata = mdata[mdata.obs["rna:leiden"] != "5"]

#misty(mdata, x_mod="rna", y_mod="cytosig", predictors=predictors, targets=targets, bandwidth=500, kernel="gaussian", bypass_intra=False, add_juxta=False, add_self=True, group_intra_by=None, group_env_by=None, overwrite=True)
misty(mdata, x_mod="rna", y_mod="cytosig", predictors=predictors, targets=targets, bandwidth=500, kernel="gaussian", bypass_intra=False, add_juxta=False, add_self=True, group_intra_by="leiden", group_env_by=None, overwrite=True)
#misty(mdata, x_mod="rna", y_mod="cytosig", predictors=predictors, targets=targets, bandwidth=500, kernel="gaussian", bypass_intra=False, add_juxta=False, add_self=True, group_intra_by=None, group_env_by="leiden", overwrite=True)
#misty(mdata, x_mod="rna", y_mod="cytosig", predictors=predictors, targets=targets, bandwidth=500, kernel="gaussian", bypass_intra=False, add_juxta=False, add_self=True, group_intra_by="leiden", group_env_by="leiden", overwrite=True)
