import pandas as pd
import numpy as np
from dask.distributed import Client
import dask.dataframe as dd

client = Client(processes=False)


mps_out = '/Volumes/aps_timemachine/1105_SE_mps_dists/'
mps_dat = dd.read_csv(mps_out+'SE_MPS_*_MASTER_XMFA_OUT.csv')
mps_dat = mps_dat.drop_duplicates(subset=['b'])
mps_dat = mps_dat.reset_index(drop=True)
mps_dat = mps_dat[mps_dat['b'] != 'all']
mps_dat.compute()
mps_dat.to_csv(mps_out+'dasktest.csv', index=False)
# mps_out = '/scratch/aps376/recombo/APS156_SP_Archive/1201_SP_mps_dists/'
# mps_dat = pd.read_csv(mps_out+'SP_MPS_0_MASTER_XMFA_OUT.csv')
# for i in tqdm(np.arange(1,500)):
#     temp = pd.read_csv(mps_out+'SP_MPS_'+str(i)+'_MASTER_XMFA_OUT.csv')
#     mps_dat = mps_dat.append(temp)
# mps_dat = mps_dat.drop_duplicates(subset=['b'])
# mps_dat = mps_dat.reset_index(drop=True)
# mps_dat = mps_dat[mps_dat['b']!='all']
#
# print(len(mps_dat))
#
# mps_dat.to_csv('APS156_201202_SP_all_dists.csv', index = False)