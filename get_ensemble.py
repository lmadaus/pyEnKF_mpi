#!/usr/bin/env python

import os
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta
import multiprocessing as mp


def main():
    init = datetime(2013,4,23,0)
    fhours = (0,6,12,18,24)
    vars = ['PSFC']
    ens, meta, mems = get_ensemble(init,fhours,vars)


# The worker class is just a process waiting
# to be assigned a task
class worker(mp.Process):
    def __init__(self, task_queue, result_queue):
        mp.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
    def run(self):
        proc_name = self.name
        while True:
            # Go to the task queue and get a new task
            next_task = self.task_queue.get()
            if next_task is None:
                # If there are no other tasks, then we are done
                break
            result = next_task()
            self.result_queue.put(result)
        return

# Here's what each ensemble-getting task will actually be
class get_mem(object):
    # Initialize whatever ensemble member we are using
    def __init__(self, basedir,starttime,ensmem,ensplace,vardict,state_vars,fdates,elev):
        self.basedir = basedir
        self.ensmem = ensmem
        self.ensplace = ensplace
        self.starttime = starttime
        self.vardict = vardict
        self.state_vars = state_vars
        self.fdates = fdates
        self.elev = elev
    def __call__(self):
        # Get the data here
        statevect = get_data(self.basedir,self.starttime,self.ensmem,self.ensplace,self.vardict,self.state_vars,self.fdates,False,self.elev)
        return (self.ensplace,statevect)



def get_data(basedir,init,m,mplace,vardict,state_vars,fdates,make_ensmeta,elev):
    mlist = []
    ensmeta = []
    # Need an outer loop here for each date in fdates
    for fdate in fdates:
        # Check for this again later
        try:
            indata = Dataset('/'.join([basedir, 'mem%03d_cm1out.nc' % m]),'r')
        except:
            print "Could not open member %s for init date %s and forecast time %s" %\
                    (m,init)
        # Get a list of the indexes for these forecast dates
        """
        model_times = indata_sfc.variables['time'][:]
        model_dates = [datetime(1970,1,1,0) + timedelta(seconds=d) for d in model_times]
        model_indexes = [model_dates.index(f) for f in fdates]
        """
        # Edit here since the forecast index is the init time
        fdex = fdate 
        ftime = fdate
        #for ftime,fdex in zip(fdates,model_indexes):
        #with ftime and fdate:
        curshape = None
        if True:
            if make_ensmeta:
                print "    ", ftime
            for v in state_vars:
                if v in vardict.keys():
                    try:
                        vardata = np.squeeze(indata.variables[vardict[v]][fdex,:,:,:])
                    except:
                        vardata = np.squeeze(indata.variables[vardict[v]][fdex,:,:])
                elif v == 'psfc':
                    #vardata = np.squeeze(indata.variables['PRES_meansealevel'][fdex,:,:])
                    vardata = np.squeeze(indata.variables['prs'][fdex,0,:,:])
                    #vardata = vardata / 100.
                    # Convert to altimeter
                    #vardata = np.power((np.power((vardata/100. - 0.3),0.190284)+np.multiply(elev,8.4228807E-5)),1/0.190284)
                    # Convert back to Pascals
                    #vardata = vardata * 100.
                else:
                    print "Error: could not open variable:", v
                    exit(1)
                #except:
                #    print "Could not load variable %s from member %s" % (v,m)
                #    exit(1)

                
                length = np.product(vardata.shape)
                varlist = np.reshape(vardata,length)
                # Put in the right units, if needed
                if v in ['psfc','ALT','MSLP']:
                    varlist = varlist / 100.
                #elif v.startswith('H_'):
                #    # Divide geopotential height by 9.81 to get raw meters
                #    varlist = varlist / 9.81
                #elif v in ('T2','DPT'):
                #    varlist = (varlist - 273.) * 9./5. + 32
                mlist = mlist + list(varlist)
                if make_ensmeta:
                    if len(vardata.shape) == 2:
                        rows = range(vardata.shape[0])
                        cols = range(vardata.shape[1])
                        cmat,rmat = np.meshgrid(cols,rows)
                        rowlist = np.reshape(rmat,length)
                        collist = np.reshape(cmat,length)
                        for r,c in zip(list(rowlist),list(collist)):
                            ensmeta.append((r,c,0,fdate,v))
                    else:
                        # Three dimensional
                        rows = vardata.shape[1]
                        cols = vardata.shape[2]
                        levs = vardata.shape[0]
                        grids = np.mgrid[0:levs, 0:rows, 0:cols]
                        lmat = grids[0]
                        rmat = grids[1]
                        cmat = grids[2]
                        levlist = np.reshape(lmat,length)
                        rowlist = np.reshape(rmat,length)
                        collist = np.reshape(cmat,length)
                        for l,r,c in zip(list(levlist),list(rowlist),list(collist)):
                            ensmeta.append((r,c,l,fdate,v))

                    # Grab mod lats and mod lnos
                    modlats = indata.variables['yh'][:] # These are now in kilometers
                    modlons = indata.variables['xh'][:] # These are now in kilometers
                    if curshape is None:
                        curshape = vardata.shape
                    elif len(vardata.shape) > curshape:
                        curshape = vardata.shape
        indata.close()
    if make_ensmeta:
        return np.array(mlist), ensmeta, curshape, modlats, modlons
    else:
        return np.array(mlist)




def get_ensemble(init,fhours,state_vars):
    from NAMELIST import basedir, mems,nmemtotal , mem_dirs
    vardict = {
               't2'   : 't2',
               'thpert'   : 'thpert',
               'tsk'   : 'tsk',
               'u10'  : 'u10',
               'v10'  : 'v10',
               'cref' : 'cref',
               'th'   : 'th'}

    # Here is the base directory
    #basedir = './ecmwf_grib'
    # And here are the ensemble members
    #mems = range(1,51) 

    # Build a list of the datetimes we'll need to grab
    #fdates = [init + timedelta(hours=n) for n in fhours]
    fdates = [init+f for f in fhours]

    # Change this into a list of wrfout filenames
    #fnames = ['%s' % n.strftime('%Y-%m-%d_%H:%M:%S') for n in fdates]

    # Get the terrain from a separate file
    #terrin = Dataset('./terrain_template.nc','r')
    #elev = np.divide(np.squeeze(terrin.variables['HGT'][0,:,:]),1.00)
    #terrin.close()
    elev = []

    # Loop through each ensemble member and load the file
    #print "Member %s" % m
    # Recall that mem_dirs is a dictionary of directories to each member "number"
    print "Building first member"
    #statevect, ensmeta, shape_vardata, modlats, modlons = get_data(basedir,init,mem_dirs[1],1,vardict,state_vars,fdates,True,elev)
    statevect, ensmeta, shape_vardata, modlats, modlons = get_data(basedir,init,1,1,vardict,state_vars,fdates,True,elev)


    # Now make the ensemble template
    ens = np.zeros((len(mems),len(ensmeta)))
    ens[0,:] = statevect 
    print np.shape(ens)
    print np.shape(ensmeta), ensmeta[-10:]

    # Set up the pool for everything else
    print "Beginning MP call for rest of ensemble"
    # CALL HERE
    num_process = mp.cpu_count() 
    tasks = mp.Queue()
    results = mp.Queue()
    processes = [worker(tasks,results) for k in xrange(num_process)]

    # Start up the processes
    for w in processes:
        w.start()
    # Populate the Queue
    for k in mems[1:]:
        curlag = k/(nmemtotal+1)
        curensnum = k%(nmemtotal+1)
        print k, curlag, curensnum
        if curlag > 0:
            curensnum += 1
        tasks.put(get_mem(basedir,init,k,k,vardict,state_vars,fdates,elev))
        #if k > 20:
        #    # Use previous ensemble
        #    tasks.put(get_mem(basedir,init-timedelta(hours=24),k-50,k,vardict,state_vars,fdates,elev))
        #else:
        #    # Use current ensemble
        #    tasks.put(get_mem(basedir,init,k,k,vardict,state_vars,fdates,elev))
    # Add in None tasks to kill all processes
    for k in xrange(num_process):
        tasks.put(None)
    # Start collecting the enseble members
    num_left = len(mems[1:])
    while num_left > 0:
        result = results.get()
        ens[mems.index(result[0]),:] = result[1]
        num_left = num_left-1
        print "Remaining:", num_left



    #ens = np.array(ens)
    ens = np.transpose(ens)
    return ens, ensmeta, mems, shape_vardata, modlats, modlons


            






if __name__ == '__main__':
    main()
