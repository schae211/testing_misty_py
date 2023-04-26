
import squidpy as sq
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import mudata as md

from liana.method.sp._misty import misty
from liana.method.sp._misty import plot_distance_weights, plot_importance, plot_contribution, plot_performance

adata = sq.datasets.visium_hne_adata()

features = ["Dbi", "3110035E14Rik", "Resp18", "Lypd1", "1500015O10Rik", "Ngef", "Ecel1", "Kcnj13", "Csrp1", "2010300C02Rik", "Atp2b4", "Ptpn4", "Bok", "Scg2", "Rgs16", "Ogfrl1", "Itm2c", "Tuba4a", "C1ql2", "Satb2", "Ptprn", "Fzd5", "Pantr1", "Tmbim1", "Lct"]

mdata = md.MuData(data={"rna": adata})

misty(mdata=mdata, x_mod="rna", bandwidth=500, predictors=features, targets=features, overwrite=True)
misty(mdata=mdata, x_mod="rna", bandwidth=500, predictors=features, targets=features, bypass_intra=True, overwrite=True)
misty(mdata=mdata, x_mod="rna", predictors=features, targets=features, bandwidth=500, bypass_intra=False, add_juxta=False, add_para=True, add_self=False, overwrite=True)





