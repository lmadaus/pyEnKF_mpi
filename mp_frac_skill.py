#!/usr/bin/env python
import multiprocessing as mp
import numpy as np

class worker(mp.Process):
    def __init__(self, task_queue, result_queue):
        mp.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
    def run(self):
        proc_name = self.name
        while True:
            # Get the next task
            next_task = self.task_queue.get()
            if next_task is None:
                break
            result = next_task()
            self.result_queue.put(result)
        return

class compute_frac(object):
    def __init__(self, pt, n, truestate, modstate, thresh, radii, ypts, xpts,\
                 periodic):
        self.pt = pt
        self.n = n
        self.truestate = truestate
        self.modstate = modstate
        self.thresh = thresh
        self.radii = radii
        self.ypts = ypts
        self.xpts = xpts
        self.periodic = periodic
    def __call__(self):
        n,true_frac,model_frac = compute_fracs(self.pt, self.n,\
                            self.truestate, self.modstate, self.thresh,\
                            self.radii, self.ypts, self.xpts, self.periodic)
        return (n, true_frac, model_frac)


def compute_fracs(pt, n, truestate, modstate, thresh, radii, ypts, xpts, periodic):
    # Compute the distance to all other points based on periodicity
    import numpy as np
    dist = find_periodic_min_distance(pt,ypts,xpts,periodic)
    # Now loop through all radii
    truth_fracs = []
    model_fracs = []
    for c,r in enumerate(radii):
        r = float(r)
        truth_local = np.sum(np.ma.masked_where(dist >= r, truestate) >= thresh[0])
        model_local = np.sum(np.ma.masked_where(dist >= r, modstate) >= thresh[1])
        npts = np.sum(np.where(dist <= r, 1,0))
        truth_frac = truth_local / float(npts)
        model_frac = model_local / float(npts)
        #if n == 1:
        #    print r, truth_local, model_local, npts
        truth_fracs.append(truth_frac)
        model_fracs.append(model_frac)
        #truths[c,n] = truth_frac
        #models[c,n] = model_frac
    return (n,truth_fracs,model_fracs)


def mp_frac_skill(modelfield, truthfield, dx, threshold, radius=None, periodic=False):
    """ Function to compute the fraction skill score of Roberts and Lean (2007)
    / Roberts (2008) for a model field and truth field.  Assumes that the model
    and truth field have already been re-gridded to be the same shape/grid.  If
    "None" is specified as the localization radius, will compute FSS for all
    integral multiples of dx from 5*dx (~minimum resolvable scale) through 50*dx
    """

    import numpy as np
    import multiprocessing as mp

    # Check to see that our shapes are identical
    if modelfield.shape != truthfield.shape:
        print "ERROR: frac_skill requires model field and truth field to be on the same grid!"
        return

    # Get dimensions
    ydim,xdim = modelfield.shape
    ypts, xpts = np.meshgrid(np.arange(ydim), np.arange(xdim))
    # Set up the verification points, depending on periodicity
    if periodic:
        # Can include all points
        verify_points = [(y,x) for y in range(ydim)[::2] for x in \
                         range(xdim)[::2]]
    else:
        # Stay away from the boundaries
        if radius is None:
            max_radius = 50
        else:
            max_radius = radius
        ydim_use = range(ydim)[max_radius:-max_radius]
        xdim_use = range(xdim)[max_radius:-max_radius]
        verify_points = [(y,x) for y in ydim_use for x in xdim_use]

    # Set radii to test
    if radius is None:
        radii = [3,4,5,6,7,8,9,10,12,15,20,25]
    else:
        radii = [radius]

    # Test our thresholds
    try:
        use_threshold = (float(threshold[0]), float(threshold[1]))
    except:
        # Must be a percentile
        use_threshold = (np.percentile(truthfield,float(threshold[0][:-3]),axis=None),
                         np.percentile(modelfield,float(threshold[0][:-3]),axis=None))



    truths = np.zeros((len(radii),len(verify_points))) 
    models = np.zeros((len(radii),len(verify_points))) 
    print "Computing FSS at all points"

    # Spin up mp here
    num_process = mp.cpu_count()
    tasks = mp.Queue()
    results = mp.Queue()
    processes = [worker(tasks,results) for k in xrange(num_process)]
    print "Using %d processes..." % num_process
    # Create a manager to hold shared arrays
    manager = mp.Manager()

    # Start the processes
    for w in processes:
        w.start()
    # Populate the queue
    for n,pt in enumerate(verify_points):
        tasks.put(compute_frac(pt, n, truthfield, modelfield, use_threshold,\
                                radii, ypts, xpts, periodic))
    # Put the kill pills
    for k in xrange(num_process):
        tasks.put(None)
    num_left = len(verify_points)
    while num_left > 0:
        result = results.get()
        n, truth_frac, model_frac = result
        truths[:,n] = truth_frac
        models[:,n] = model_frac
        num_left = num_left - 1
        if num_left % 1000 == 0:
            print "Remaining:", num_left


    # Now go back and compute the total FSS
    outd = {}
    num = len(verify_points)
    #FBS = 1/float(num) * np.sum([(o-m)**2 for o,m in zip(truths[r],models[r])])
    #FBS_worst = 1/float(num) * (np.sum([o**2 for o in truths[r]]) + \
    #                            np.sum([m**2 for m in models[r]]))
    FBS = 1.0 / num * np.sum(np.power(np.subtract(truths, models), 2), axis=1)
    FBS_worst = 1.0 / num * np.add(np.sum(np.power(truths, 2), axis=1),
                                       np.sum(np.power(models, 2), axis=1))


    # Output in terms of raw distance (r*dx)
    for c,r in enumerate(radii):
        outd[r*dx] = 1.0 - (FBS[c] / FBS_worst[c])
    return outd






def find_periodic_min_distance(pt,ygrid,xgrid,isperiodic):
    """ Return an array of the minimum distance between a point and all other
    points in gridpoint dimensions """
    import numpy as np
    dist = np.sqrt(np.add(np.power(ygrid-pt[0],2),np.power(xgrid-pt[1],2)))
    if isperiodic:
        ny,nx = ygrid.shape
        # Nine possible locations to test
        for testpt in [pt, (pt[0]+ny,pt[1]-nx), (pt[0]+ny,pt[1]),
                       (pt[0]+ny,pt[1]+ny), (pt[0],pt[1]-nx), (pt[0],pt[1]+nx),
                       (pt[0]-ny,pt[1]-nx), (pt[0]-ny,pt[1]),
                       (pt[0]-ny,pt[1]+nx)]:
            
            newdist = np.sqrt(np.add(np.power(ygrid-testpt[0],2),np.power(xgrid-testpt[1],2)))
            dist = np.minimum(dist,newdist)
    return dist



if __name__ == '__main__':
    # Actually try this here
    import numpy as np
    #from frac_skill_cython import frac_skill_cython
    threshold = ('95pct','95pct') # Anomaly to try and resolve
    dx = 1.0 # Here, dx in km
    modtim = 95 # Model time step to verify
    var = 't2'
    method = '3d'
    #gridres = 10
    grids = [10,4,2,1]
    outd = {}
    for gridres in grids:
        print "GRID:", gridres
        filename = './50mems/time_%03d_%02dkmgrid_%s_ensemble_output.pckl' % (modtim,
                                                                     gridres,
                                                                     method)
        # Try loading the files
        import cPickle
        from netCDF4 import Dataset
        modfil = open(filename,'r')
        moddict = cPickle.load(modfil)
        modfil.close()
        moddat = moddict[modtim]['%s Post Mean' % var]
        # Now truth
        truth = Dataset('/home/disk/pvort/lmadaus/nobackup/cm1/enkf_cm1/run/kffc_withpbl/cm1out_s.nc','r')
        verify = truth.variables[var][modtim,:,:]
        # Need to re-grid this to the same grid
        verify_regrid = verify[::5,::5]
        moddat_regrid = moddat[0:101, 0:101]
        print "Mod shape:", moddat_regrid.shape
        print "Verify shape:", verify_regrid.shape
        # Remove the means
        moddat_regrid = np.subtract(moddat_regrid,np.mean(moddat_regrid,axis=None,dtype=np.float64))
        verify_regrid = np.subtract(verify_regrid,np.mean(verify_regrid,axis=None,dtype=np.float64))
        print "VAR:", np.var(moddat_regrid, axis=None, dtype=np.float64)
        moddat_regrid = np.abs(moddat_regrid)
        verify_regrid = np.abs(verify_regrid)
        moddat_regrid = moddat_regrid.astype(np.float64)
        verify_regrid = verify_regrid.astype(np.float64)
        # Now the FSS
        #FSS = frac_skill(moddat_regrid, verify_regrid, dx, threshold, periodic=True)
        FSS = mp_frac_skill(moddat_regrid, verify_regrid, dx, threshold, periodic=True)
        #FSS = frac_skill_cython(moddat_regrid, verify_regrid, dx, threshold, periodic=True)
        outd[gridres] = FSS

    # Now plot
    print "Plotting"
    import matplotlib
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8,8))
    for res in grids:
        radii = outd[res].keys()
        radii.sort()
        scores = [outd[res][r] for r in radii]
        plt.plot(radii, scores, linewidth=3,label='%02dkm' % res)
    plt.axhline(y=0.5, linewidth=2, linestyle='dotted',c='k')
    plt.ylim((0,1))
    plt.legend(loc=0)
    plt.ylabel('Fractions Skill Score')
    plt.xlabel('Radius (km)')
    plt.title('Time %03d %s>%s obs grid FSS' % (modtim, var, str(threshold[0])))
    plt.grid()
    plt.show()
    #plt.savefig('fss_%s_%s.png' % (var, method),bbox_inches='tight')
