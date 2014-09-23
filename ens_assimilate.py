#!/usr/bin/env python
"""
This is the master run script.  It will import values from the 
NAMELIST.py file for use in the ensemble update.
Luke Madaus (Feb 2014)
"""
import cPickle
import numpy as np
from datetime import datetime, timedelta
export_flag =True 
import matplotlib
if export_flag:
    matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.axes as maxes
from mpl_toolkits.axes_grid import make_axes_locatable
from netCDF4 import Dataset
import os
from scipy.stats import nanmean
from NAMELIST import *
from mpi4py import MPI
from alt_mpi_enkf_update import mpi_filter

# Translate observations to state variables
state_equiv = {'psfc' : 'psfc',
               'alt' : 'ALT',
               'mslp' : 'MSLP',
               'u10' : 'u10',
               'v10' : 'v10',
               'precip' : 'PRECIP',
               'temp' : 't2',
               't2' : 't2',
               'thpert' : 'thpert',
               'tsk' : 'tsk',
               'cref' : 'cref',
               'dewp' : 'DPT',
               'h_500hPa' : 'H_500hPa',
               'temp_full' : 'th'}


# Plotting levels for consistency
diff_plot_levs = {'ALT' : (-2,2),
                  'PSFC' : (-2,2),
                  'MSLP' : (-2,2),
                  'PRECIP' : (-10,10),
                  'U10' : (-5,5),
                  'V10' : (-5,5),
                  'T2' : (-2,2),
                  'DPT': (-2,2),
                  'H_500hPa' : (-50,50)}

# Set dx here for now (in km)
dx = 1.0

def main():
    # Accumulate all ensemble members into the state vector
    from get_ensemble import get_ensemble
    ens, ensmeta, mems, map_shape, modlats, modlons = get_ensemble(starttime,fcst_hours,[state_equiv[v] for v in state_vars])
    # Turn ensmeta into a tuple
    ensmeta = tuple(ensmeta)
    # Quick print out of the gridpoint dimensions of the ensemble
    print "Gridpoint dimensions:", map_shape

    # Get the number of state variables and the number of members found
    Nstate = ens.shape[0]
    Nmems = ens.shape[1]
    print "Nstate:", Nstate
    print "Nmems:", Nmems
    print "Len meta:", len(ensmeta)
    # Print some example ensemble values to be sure these are not NaNs
    print "Sample ens. values:", ens[0,:]

    # Get the mean and perturbations
    xfm = np.mean(ens,axis=1,dtype=np.float64)  # Ensemble mean
    #xfmtile = np.tile(xfm,(Nmems,1)) 
    #xfmtile = np.transpose(xfmtile)
    #xfp = np.subtract(ens,xfmtile)
    xfp = np.subtract(ens,xfm[:,None]) # Broadcast the mean and subtract to get perts
    # Try inflating here
    #xfp = np.multiply(xfp,inflation_fact)


    # Start a timer to keep track of how long the next set of steps takes 
    from time import time
    start_bench = time()
    # Accumulate the observations
    # Figure out which format to use
    if ob_format == 'MADIS':
        from format_obs import format_obs
        obd = format_obs(starttime,first_obtime,first_obtime+timedelta(hours=assim_lead),assim_vars,assim_places,state_equiv)
    elif ob_format == 'TEXT':
        from format_obs import format_obs_text
        obd = format_obs_text(starttime,first_obtime,first_obtime+timedelta(hours=assim_lead),assim_vars,assim_places,state_equiv)
    elif ob_format == 'IDEAL':
        from format_ideal_obs import format_ideal_obs
        obd = format_ideal_obs(starttime,first_obtime,first_obtime+assim_lead,assim_vars,assim_places,state_equiv,modlats,modlons)
    obd_bench = time()
    print "Getting obs time:", obd_bench - start_bench
    print "Number of obs to assimilate:", len(obd.keys())
    
    # Option here to thin the number of obs used
    #print "thinning obs"
    #newob = obd.copy()
    #obkeys = obd.keys()
    #filtered = obkeys[1:3]
    ##filtered = obkeys
    #for k in obd.keys():
    #    if k not in filtered:
    #        del newob[k]
    #del obd
    #obd = newob.copy()
    #del newob
    #print "Number of obs to assimilate:", len(obd.keys())
    #new_obd = obd.copy()
    #new_obd = {}
    #keysort = obd.keys()
    #keysort.sort()
    #for k in [('59.5/59.5',100,'temp'),('79.3/19.9',100,'temp')]:
    #for k in [('79.3/19.9',100,'temp')]:
    #    new_obd[k] = obd[k]
    #    print obd[k]
    #obd = new_obd.copy()
    #raw_input()
    
    """

    # Dump out the actual observation dictionary
    #outf = open('%s_obs_assimilated.pckl' % starttime.strftime('%Y%m%d%H'),'wb')
    print "Dumping observations"
    outf = open('observations.pckl','wb')
    cPickle.dump(obd,outf)
    outf.close()

    # Dump the state
    print "Dumping state"
    outf = open('prior_state.pckl','wb')
    cPickle.dump((xfm,xfp),outf)
    outf.close()

    # Dump the state
    print "Dumping metadata"
    outf = open('ensmeta.pckl','wb')
    cPickle.dump((ensmeta,modlats,modlons,state_equiv,loc_rad[0],dx),outf)
    outf.close()




    # Read the state
    print "Reading updated state"
    outf = open('posterior_state.pckl','rb')
    xam, Xap = cPickle.load(outf)
    outf.close()


    # Import the Filter and actually do the assimilation
    from wrapper import wrap_enkf
    #inflation_fact = 1.0
    xam, Xap = wrap_enkf(xfm,xfp,obd,inflation_fact,loc_rad[0],loc_rad[1],ensmeta,state_equiv,modlats,modlons)
    # Make sure we're not full of NaNs
    #print "xam shape:", xam.shape
    print "xam nans:", np.sum(np.isnan(xam),axis=None)
    #print "Xap shape:", Xap.shape
    print "Xap nans:", np.sum(np.isnan(Xap),axis=None)
    #print "After", xam
    """
    # Now this will be an os call
    #os.system('mpirun -np 16 ./mpi_enkf_update.py')
    #os.system('mpirun -np 32 ./alt_mpi_enkf_update.py')
    xam,Xap = mpi_filter(obd,(ensmeta,modlats,modlons,state_equiv,loc_rad[0],dx),xfm,xfp,rank=0)
    end_bench = time()
    print "Assimilation time:", end_bench - obd_bench
    #comm.Disconnect()
    # Add the mean back in to get the full analysis ensemble
    #xamtile = np.tile(xam,(Nmems,1))
    #xamtile = np.transpose(xamtile)
    Xa = np.add(xam[:,None],Xap)
    print "Shape of Xap:", Xap.shape
    #print "Random Xap:", Xap[20,:]
    # Compute the spread
    spread_analysis = np.std(Xap,axis=1,dtype=np.float64)
    spread_forecast = np.std(xfp,axis=1,dtype=np.float64)
    print "Sample spread:", spread_analysis[20]

    # Archive rejected obs file if it exists
    if os.path.exists('rejected_obs.pickle'):
        os.system('mv rejected_obs.pickle %s_rejected_obs.pckl' % verify_time.strftime('%Y%m%d%H'))

    # Archive the spread for the verification time
    # Output full fields as a dictionary
    try:
        if len(fcst_hours) > 1:
            if 'cref' in state_vars:
                # We must be doing 4dassim
                infile = open('time_%03d_%02dkmgrid_4d_efa_ensemble_output.pckl' %
                       (starttime,int(truth_obsloc[4:-2])),'r')
            else:
                # We must be doing 4dassim
                infile = open('time_%03d_%02dkmgrid_4d_ensemble_output.pckl' %
                       (starttime,int(truth_obsloc[4:-2])),'r')
            #os.system('mv posterior_state.pckl time_%03d_%02dkmgrid_4d_posterior.pckl' %\
            #            (starttime, int(truth_obsloc[4:-2])),'r')
        else:
            # We must be doing 3dassim
            infile = open('time_%03d_%02dkmgrid_3d_ensemble_output.pckl' %
                       (starttime,int(truth_obsloc[4:-2])),'r')
            #os.system('mv posterior_state.pckl time_%03d_%02dkmgrid_3d_posterior.pckl' %\
            #            (starttime, int(truth_obsloc[4:-2])),'r')
        outputd = cPickle.load(infile)
        infile.close()
    except:
        outputd = {}

    def is_match(x,var,time):
        if x[-1] == var and x[-2] == time:
            return True
        else:
            return False


    for var in state_vars:
        for fhour in [starttime + f for f in fcst_hours]:
            verify_time = fhour
            verify_time_int = int(verify_time)
            #outputd[verify_time] = {}
            if verify_time not in outputd.keys():
                outputd[verify_time] = {}
            state_var_name = state_equiv[var]
            # Find locations in statemeta that match
            use_these = np.array([is_match(z, state_var_name, verify_time_int) for z in ensmeta])
            """
            prior_spread_field = [x for (x,z) in zip(list(spread_forecast),ensmeta) if z[-1] == state_var_name and z[-2] == verify_time_int]
            post_spread_field = [x for (x,z) in zip(list(spread_analysis),ensmeta) if z[-1] == state_var_name and z[-2] == verify_time_int]
            prior_mean_field = [x for (x,z) in zip(list(xfm),ensmeta) if z[-1] == state_var_name and z[-2] == verify_time_int]
            post_mean_field = [x for (x,z) in zip(list(xam),ensmeta) if z[-1] == state_var_name and z[-2] == verify_time_int]
            prior_perts_ = [x for (x,z) in zip(list(xfm),ensmeta) if z[-1] == state_var_name and z[-2] == verify_time_int]
            post_perts = [x for (x,z) in zip(list(xam),ensmeta) if z[-1] == state_var_name and z[-2] == verify_time_int]
            """
            #prior_spread_field = spread_forecast[use_these]
            print verify_time, state_var_name, spread_forecast[use_these].shape, map_shape
            try:
                print "Trying 3d"
                outputd[verify_time]['%s Prior Spread' % state_var_name] = np.reshape(spread_forecast[use_these], map_shape)
                outputd[verify_time]['%s Post Spread' % state_var_name] = np.reshape(spread_analysis[use_these], map_shape)
                outputd[verify_time]['%s Prior Mean' % state_var_name] = np.reshape(xfm[use_these], map_shape)
                outputd[verify_time]['%s Post Mean' % state_var_name] = np.reshape(xam[use_these], map_shape)
                outputd[verify_time]['%s Post Mems' % state_var_name] = \
                        np.reshape(Xa[use_these,:],(map_shape[0],map_shape[1],map_shape[2],Nmems))
                outputd[verify_time]['%s Prior Mems' % state_var_name] = \
                        np.reshape(ens[use_these,:],(map_shape[0],map_shape[1],map_shape[2],Nmems))
            except:
                print "Trying 2d"
                outputd[verify_time]['%s Prior Spread' % state_var_name] = np.reshape(spread_forecast[use_these], (map_shape[0],map_shape[1]))
                outputd[verify_time]['%s Post Spread' % state_var_name] = np.reshape(spread_analysis[use_these], (map_shape[0],map_shape[1]))
                outputd[verify_time]['%s Prior Mean' % state_var_name] = np.reshape(xfm[use_these], (map_shape[0],map_shape[1]))
                outputd[verify_time]['%s Post Mean' % state_var_name] = np.reshape(xam[use_these], (map_shape[0],map_shape[1]))
                outputd[verify_time]['%s Post Mems' % state_var_name] = np.reshape(Xa[use_these,:],(map_shape[0],map_shape[1],Nmems))
                outputd[verify_time]['%s Prior Mems' % state_var_name] = np.reshape(ens[use_these,:],(map_shape[0],map_shape[1],Nmems))

    # Dump the output dictionary
    if len(fcst_hours) > 1:
        if 'cref' in state_vars:
            # We must be doing 4dassim
            outfile = open('time_%03d_%02dkmgrid_4d_efa_ensemble_output.pckl' %
                   (starttime,int(truth_obsloc[4:-2])),'w')
        else:
            # We must be doing 4dassim
            outfile = open('time_%03d_%02dkmgrid_4d_ensemble_output.pckl' %
                   (starttime,int(truth_obsloc[4:-2])),'w')
    else:
        # We must be doing 3dassim
        outfile = open('time_%03d_%02dkmgrid_3d_ensemble_output.pckl' %
                   (starttime,int(truth_obsloc[4:-2])),'w')
    cPickle.dump(outputd,outfile)
    outfile.close()
    exit()

    # Compute the verification at the verification time
    print "Computing verification"
    from plot_output import plot_verification
    plot_verification(xfm,xam,spread_forecast,spread_analysis,ensmeta,state_vars,map_shape,obd,basedir)


if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    size = comm.size
    rank = comm.rank
    status = MPI.Status()

    if rank == 0:
        # We are the master -- go to main
        main()
    else:
        # Start up the filter and wait for instructions
        dum = mpi_filter(None, None, None, None, rank=rank)
