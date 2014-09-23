#/usr/bin/env python
from datetime import datetime, timedelta
import cPickle

# A list of variables used by this assimilation system


# THESE ARE ALL CURRENTLY SUPPORTED STATE VARIABLES/OBS
# SURFACE: 'alt','psfc','mslp','uwnd','vwnd','precip','temp','dewp'
# UPPER AIR: 'h_500hPa'


# List of variables to be included in the state vector
state_vars = ['t2','tsk','thpert']

# List of observations to assimilate.  Must be a subset of
# the state_vars list
assim_vars = ['t2']

# List of specific observations IDs to assimilate.  If None, all
# observation locations in the input file will be considered
assim_places = None
#assim_places = cPickle.load(open('asos_lat_lon.pickle','r')).keys()
#assim_places = ['KUIL']

# Thinning factor used to reduce observations.  If no thinning
# desired, this should be 1.
thin_factor = 1


# Forecast time of ensemble start (Ensemble initilization time)
# This is now the index of the time variable
starttime =85 

# Time of last ensemble forecast to be updated
#endtime = datetime(2014,4,29,18)

# The actual forecast hours to be included in the state vector
#fcst_hours = [-1,0,1]
fcst_hours = [0]
#fcst_hours = [0,1,2,3,4,5,6]

# What time to start the assimilation (must be included in fcst_hours)
#first_obtime = -1
first_obtime = fcst_hours[0]
#first_obtime = starttime + timedelta(hours=fcst_hours[0])

# How many subsequent hours to continue the assimlation
# If zero, only assimilate at the first_obtime
#assim_lead = len(fcst_hours)-1 if len(fcst_hours)>1 else 0
assim_lead = 0 

# Specify which observation parser to use. Set to 'MADIS' for MADIS netcdf
# observation files or "TEXT" for /home/disk/data CSV text files
# "IDEAL" to draw ideal obs from another file (see "truthdat" below)
ob_format = 'IDEAL'

# Timeis to use for verification of the forecast (must be included in fcst_hours)
#verify_time = starttime + timedelta(hours=fcst_hours[-1])

# Manally specify an inflation value here
inflation_fact = 1.0 

# Horizontal localization radius here (in km)
# Must be a tuple of two values; only the first
# value is used.  This is the localization half-width.
# Full cutoff will be twice this value
loc_rad = (20.,1.0)


# Path to location of the ensemble data
basedir = '/home/disk/trof2/lmadaus/nobackup/kffc_ensemble'

# Path to truth data
truthdat = '/home/disk/pvort/lmadaus/nobackup/cm1/enkf_cm1/run/kffc_withpbl/cm1out.nc'
truth_obsloc = 'GRID2km'


# Number of ensemble members to look for (before any lagging)
nmemtotal = 100 
# See special section below if subdirectories are different

# Number of time lags in ensemble (1 means only the most recent ensemble)
nlag = 1

# Number of ensemble members (in a list)  DO NOT EDIT
mems = range(1,nmemtotal*nlag+1)



""" Edit this section to look for different subdirectories 
This will associate each member "number" with a different subdirectory name"""
#models = ('cmcg','gasp','gfs','jma','nam','ngps','tcwb','ukmo')
#mem_dirs = dict(zip([x + 1 for x in range(len(models))], models))
#mem_dirs = dict(zip([x+1 for x in range(len(models))], ''))
mem_dirs = {}



