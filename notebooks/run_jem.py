import os
import sys
import pickle
import numpy as np
import pandas as pd

from jemm.genomic_annotation import ExonSet
from jemm.junction import JunctionCountTable
from jemm.transcript import TranscriptMeasureTable
from jemm.covariate import Contrasts, Covariate
from jemm.model import JemmLinearMixedModel

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--event_type", type=str, choices=['SE', 'RI', 'A5SS', 'A3SS'])
parser.add_argument("--sex", type=str, choices=['M', 'F'])
#parser.add_argument("--data", type=str, choices=['txr', 'jct', 'both'], default="txr")
args = parser.parse_args()

DATA_VER = 'data-V9'
EVENT_TYPE = args.event_type
SEX = args.sex
#DATA_TO_USE = args.data
DATA_TO_USE = 'both' if EVENT_TYPE=="SE" else 'txr'
MIN_RE_VAR = 0.001
OUTDIR = './%s/sex_stratified.%s/' % (DATA_VER, EVENT_TYPE)

SWAP_630_634 = True

assert EVENT_TYPE in ('SE', 'A5SS', 'A3SS', 'RI')
assert SEX in ('M', 'F')

print(EVENT_TYPE, SEX)

Jemm = JemmLinearMixedModel
# load annotations
exonset = ExonSet.from_suppa("/mnt/ceph/users/zzhang/SUPPA/index/hg38/suppa_gencodev34_%s_strict.ioe"%EVENT_TYPE, 
     cache_dir='./data/',
     event_type=EVENT_TYPE)
# contrasts and covariates
contrasts = Contrasts(name="final", levels=[
    'Control',
    'First',
    'Mid',
]) + Contrasts(name="Sex", levels=[SEX])
covs = Covariate(fp="./%s/charm_master.no_seropos.csv" % DATA_VER, sep=",",
                 index_col=0,
                 contrasts=contrasts,
                 main_effects=['final', 'plateNum', 'pid'],
                 factor_conversion={
                     'final': {
                         'First': 'final@First',
                         'Mid': 'final@Mid',
                     },
                     'plateNum': {
                         'P2R': 'plateNum@P2R',
                         'P3': 'plateNum@P3',
                         'P4': 'plateNum@P4',
                         'P5': 'plateNum@P5',
                         'P6': 'plateNum@P6',
                         'P7': 'plateNum@P7',
                         'P8': 'plateNum@P8',
                         'P9': 'plateNum@P9',
                         'P10': 'plateNum@P10',
                         'P12': 'plateNum@P12',
                         'P13': 'plateNum@P13',
                         'P14': 'plateNum@P14',
                         'P15': 'plateNum@P15',
                         'P16': 'plateNum@P16',
                         'P17': 'plateNum@P17',
                         'P18': 'plateNum@P18',
                         'P19': 'plateNum@P19',
                         'P20': 'plateNum@P20',
                         'P21': 'plateNum@P21',
                         'P22': 'plateNum@P22',
                         'P23': 'plateNum@P23',
                         'P24': 'plateNum@P24',
                         'P25': 'plateNum@P25',
                     },
                 },
                 verbose=True
             )
print(covs.formula)

# load data and build object
if EVENT_TYPE == "SE":
    JCT = JunctionCountTable.from_plaintext("./%s/compiled/jct_%s.txt" % (DATA_VER, EVENT_TYPE) )
    TXR = TranscriptMeasureTable.from_plaintext("./%s/compiled/txr_%s.txt" % (DATA_VER, EVENT_TYPE) )
else:
    JCT = pickle.load(open("./%s/compiled/jct_%s.pkl" % (DATA_VER, EVENT_TYPE), "rb"))
    TXR = pickle.load(open("./%s/compiled/txr_%s.pkl" % (DATA_VER, EVENT_TYPE), "rb"))

# there was a sample id/fastq filename swap between v7 -> v8/9
# hence, I will use 20_0634-T56 to overwrite 20_0630-T56
# zzjfrank, Aug 19, 2021
if SWAP_630_634 is True:
    JCT.data['20_0630-T56'] = JCT.data['20_0634-T56']
    TXR.data['20_0630-T56'] = TXR.data['20_0634-T56']

jem = Jemm(
    junction_measure=JCT,
    transcript_measure=TXR,
    covariates=covs,
    diff_intercept_by_measure=True,
    group_varname='pid',
    min_groupvar=MIN_RE_VAR,
    optimizer='bfgs'
)

# run tests
_ = jem.run_tests(test_type='Wald',
                  data=DATA_TO_USE,
                  force_diff_intercept=False,
                  force_rerun=True,
                  pval_adjust_method="fdr",
                  nthreads=1
                 )

jem.munge_covariates([
    'final@First',
    'final@Mid'],
    meta_name = "final2cond"
)

outfp = os.path.join(OUTDIR, "%s.%s.reg_table.tsv"%('male' if SEX=="M" else 'female', EVENT_TYPE))
df = jem.save_regression_table(outfp,
   exonset=exonset,
   annotations=None,
   order_by_covariate='final2cond',
   order_by='logP'
)
