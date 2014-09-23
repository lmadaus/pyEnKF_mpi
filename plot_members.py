#!/usr/bin/env python

# In[37]:

import cPickle
import numpy as np
import matplotlib.pyplot as plt
from frac_skill import frac_skill

from netCDF4 import Dataset

time = 80 
res = 1 
assimtype = '3d'
var = 't2'
threshold = 1.5

infile = open('time_%03d_%02dkmgrid_%s_ensemble_output.pckl' % (time, res, assimtype),'rb')
data = cPickle.load(infile)
infile.close()

mems = data[time]['%s Post Mems' % var]
print mems.shape


# In[38]:

truth = Dataset('/home/disk/pvort/lmadaus/nobackup/cm1/enkf_cm1/run/kffc_withpbl/cm1out_s.nc','r')
verify = truth.variables[var][time,:,:]
verify_regrid = verify[::5,::5]
verify_regrid = np.subtract(verify_regrid, np.mean(verify_regrid, axis=None, dtype=np.float64))
verify_regrid = np.abs(verify_regrid)


# In[ ]:

outd = {}
for mem in xrange(mems.shape[-1]):
    moddat = mems[:,:,mem]
    moddat_regrid = moddat[0:101, 0:101]
    moddat_regrid = np.subtract(moddat_regrid, np.mean(moddat_regrid, axis=None, dtype=np.float64))
    moddat_regrid = np.abs(moddat_regrid)
    FSS = frac_skill(moddat_regrid, verify_regrid, 1.0, threshold, periodic=True)
    outd[mem] = FSS
    
    


# In[33]:

values = outd[0].keys()
values.sort()
#values_actual = [v/10. for v in values]


# In[34]:
"""
# Now for plotting
plt.figure(figsize=(12,12))
for mem in outd.keys():
    plotseq = [outd[mem][k] for k in values]
    plt.plot(values,plotseq)
plt.title('Time %03d %s %dkm grid Skill Scores' % (time,var,res))
plt.xlabel('Radius (km)')
plt.ylabel('Fractions Skill Score')
plt.ylim((0,1))
plt.axhline(y=0.5,linewidth=2,c='k',linestyle='dashed')
plt.show()
"""

# In[35]:

# Now archive
try:
    outfile = open('./50mems/frac_skills.pckl','rb')
    outskills = cPickle.load(outfile)
    outfile.close()
except:
    outskills = {}
if res not in outskills.keys():
    outskills[res] = {}
if time not in outskills[res].keys():
    outskills[res][time] = {}
if var not in outskills[res][time].keys():
    outskills[res][time][var] = {}
if threshold not in outskills[res][time][var].keys():
    outskills[res][time][var][threshold] = {}
outskills[res][time][var][threshold][assimtype] = outd


# In[36]:

outfile = open('./50mems/frac_skills.pckl','wb')
cPickle.dump(outskills,outfile)
outfile.close()


# In[ ]:



