#!/usr/bin/env python

from mpi4py import MPI
import numpy as np
from scipy import linalg
import cPickle
from util_wrap import stat_sig_local, cyth_efficient_localization

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), *named)
    return type('Enum', (), enums)

def mpi_filter(obs, ensmetas, xbm, Xbp, rank):
    comm = MPI.COMM_WORLD # Get communicator
    #rank = comm.rank # Which processor is this

    # Rank 0 process is collecting and distributing
    if rank == 0:
        xam, Xap = master(obs, ensmetas, xbm, Xbp, rank)
        return xam, Xap
    else:
        worker(rank)
        return None

def master(obs, ensmetas, xbm, Xbp,rank=0):
    """ Function performed by the master node """

    """
    # Broadcast ensmeta to all workers
    print "Broadcasting..."
    statemeta = ensmetas[0]
    modlats = ensmetas[1]
    modlons = ensmetas[2]
    state_equiv = ensmetas[3]
    localization=ensmetas[4]
    metadict = dict(zip(statemeta, range(len(statemeta))))
    stateys = np.array([x[0] for x in statemeta])
    statexs = np.array([x[1] for x in statemeta])
    latlist = np.array(modlats, dtype=np.float64)
    lonlist = np.array(modlons, dtype=np.float64)
    comm.bcast(metadict, root=0)
    comm.bcast(stateys, root=0)
    comm.bcast(statexs, root=0)
    comm.bcast(latlist, root=0)
    comm.bcast(lonlist, root=0)
    comm.bcast(state_equiv, root=0)
    comm.bcast(localization, root=0)
    print "Done"
    """
    comm = MPI.COMM_WORLD # Get communicator
    rank = comm.rank # Which processor is this
    size = comm.size # How many processors we have
    status = MPI.Status() # Status object
    tags = enum('READY','DONE','DONEOB','EXIT','START','AWAITSTATE','UPDATESTATE')
    # Now start figuring out the tasks
    obnames = obs.keys()
    tasks = len(obnames)
    task_index = 0
    num_workers = size - 1 # Rank 0 process not included as a worker
    closed_workers = 0
    state_hold = False
    awaiting_state = []
    obmeta = []
    obdat = []
    # Allocate the analysis state
    xam = np.zeros(xbm.shape)
    Xap = np.zeros(Xbp.shape)

    # Figure out the maximum y grid point and x grid point
    maxy = len(ensmetas[1])
    maxx = len(ensmetas[2])
    state_equiv = ensmetas[3]
    statemeta = ensmetas[0]
    metadict = dict(zip(statemeta, range(len(statemeta))))

    # Split the state by the number of workers
    lenstate = Xbp.shape[0]
    Nstate = lenstate
    # Don't add one here
    split_len = lenstate / (num_workers)
    state_sent = 1

    # Make a state array for the ensemble obs estimates
    obXbp = np.zeros((tasks,Xbp.shape[1]))
    obXbm = np.zeros(tasks)

    # Start by computing the observation priors
    print "Processing observations for assimilation..."
    for obnum, ob in enumerate(obnames):
        if obnum % 200 ==0:
            print "    Ob number:", obnum
        H = np.zeros((1,Nstate))
        weight_dict = obs[ob]['weights']
        weight_sum = np.sum([1./weight_dict[k] for k in weight_dict.keys()])
        locy = 0
        locx = 0
        for k in weight_dict.keys():
            statevar = (k[0], k[1], obs[ob]['zgrid'], ob[1],
                        state_equiv[ob[2]])
            state_loc = metadict[statevar]
            #print "   ",statevar, state_loc
            H[0,state_loc] = (1./weight_dict[k])/weight_sum
            locx += k[1] * (1./weight_dict[k])/weight_sum
            locy += k[0] * (1./weight_dict[k])/weight_sum
        obyname,obxname = ob[0].split('/')
        #if obnum % 200 == 0:
        #    print "       locx:", locx, " ob0:", float(obxname), " obx:", obs[ob]['xgrid']
        #    print "       locy:", locy, " ob0:", float(obyname), " oby:", obs[ob]['ygrid']
        # Now compute prior
        obXbp[obnum,:] = np.dot(H,Xbp)
        obXbm[obnum] = np.dot(H,xbm)
        #print 'Check:', obs[ob]['value'], obXbm[obnum]
        # Format is y,x,z,time,value,error
        #obmeta.append((obs[ob]['ygrid'],obs[ob]['xgrid'],0,ob[1],obs[ob]['value'],obs[ob]['oberr']))
        obmeta.append((float(obyname),float(obxname),0,ob[1],obs[ob]['value'],obs[ob]['oberr']))




    print "Beginning MPI filter with %d workers" % num_workers
    while closed_workers < num_workers:
        # Get the state from a worker
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG,
                         status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        #print "Received tag!", source, tag
        if tag == tags.DONE:
            # Have updated their chunk of the state, reconstruct
            startdex, len_chunk, chunk_xam, chunk_Xap = data
            xam[startdex:startdex+len_chunk] = chunk_xam
            Xap[startdex:startdex+len_chunk,:] = chunk_Xap
            # Out of obs and states so send the worker a signal to shut down
            comm.send(None, dest=source, tag=tags.EXIT)
            print "Worker", source, "has submitted updated state"
        elif tag == tags.READY:
            if state_sent <= num_workers:
                # Send chunks of state to workers
                chunk = {}
                state_start = (state_sent-1) * split_len
                chunk['state start'] = state_start
                if state_sent == num_workers:
                    # If we didn't divide equally, the last worker will pick up
                    # the extra slack
                    xbm_chunk = xbm[state_start:]
                    Xbp_chunk = Xbp[state_start:,:]
                    meta_chunk = ensmetas[0][state_start:]
                else:
                    xbm_chunk = xbm[state_start:state_start + split_len]
                    Xbp_chunk = Xbp[state_start:state_start + split_len,:]
                    meta_chunk = ensmetas[0][state_start:state_start + split_len]
                chunk['len state'] = len(meta_chunk)
                #print "split_len:", split_len
                #print "Sending worker", source, "chunk of length",len(meta_chunk)
                # Add in the obs
                chunk['meta'] = meta_chunk + tuple(obmeta)
                #print "XBM_CHUNK:", xbm_chunk.shape
                #print "obXbm:", obXbm.shape
                chunk['xbm'] = np.concatenate((xbm_chunk, obXbm), axis=0)
                chunk['Xbp'] = np.concatenate((Xbp_chunk, obXbp), axis=0)
                # Other parameters
                chunk['localization'] = ensmetas[4]
                chunk['dx'] = ensmetas[5]
                chunk['maxy'] = maxy
                chunk['maxx'] = maxx
                comm.send(chunk, dest=source, tag=tags.UPDATESTATE)
                state_sent += 1
            else:
                comm.send(None, dest=source, tag=tags.EXIT)
        elif tag == tags.EXIT:
            print("Worker %d exited." % source)
            closed_workers += 1
        """
        elif tag == tags.DONE:
            # Worker is done.  Get the updated state
            xbm, Xbp = data
            # Let us know the state is no longer held
            state_hold = False
            #print("Returned state from worker %d" % source)
            # Send the state to the next worker
            if len(awaiting_state) != 0:
                next_worker = awaiting_state[0]
                awaiting_state.pop(0)
                comm.send((xbm,Xbp), dest=next_worker, tag=tags.UPDATESTATE)
                state_hold = True
        """
    print "End of master loop"
    # Dump the posterior state here
    #outfile = open('posterior_state.pckl','wb')
    #cPickle.dump((xam,Xap), outfile)
    #outfile.close()
    print "Done."
    return xam, Xap

def worker(rank):
    """ This function defines what the worker does """
    comm = MPI.COMM_WORLD # Get communicator
    rank = comm.rank # Which processor is this
    size = comm.size # How many processors we have
    name = MPI.Get_processor_name()
    status = MPI.Status() # Status object
    tags = enum('READY','DONE','DONEOB','EXIT','START','AWAITSTATE','UPDATESTATE')
    #print "Worker %d on %s is running" % (rank, name)

    while True:
        # Start by sending a ready tag
        comm.send(None, dest=0, tag=tags.READY)
        #print "Worker %d is ready" % rank
        #print "Worker %d has localization %f" % (rank, localization)
        datachunk = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()

        # Options here
        if tag == tags.UPDATESTATE:
            print "Worker %d on %s has state" % (rank, name)
            # Unpack our data chunk
            start_index = datachunk['state start']
            # Remember---this length is the length of the state vector WITHOUT
            # the added observation estimates
            Nstate = datachunk['len state']
            state_meta = datachunk['meta']
            xbm = datachunk['xbm']
            Xbp = datachunk['Xbp']
            # Other parameters
            loc = datachunk['localization']
            dx = datachunk['dx']
            maxx = datachunk['maxx']
            maxy = datachunk['maxy']

            # Figure out where observations start
            obs_start = Nstate
            Nfull = len(state_meta)
            # Get lists of the y and x values for all points in the state
            # (including obs) for localizing
            yvals = np.array([v[0] for v in state_meta])
            xvals = np.array([v[1] for v in state_meta])
            #if rank == 2:
            #    print "NSTATE:", Nstate, "Nfull:", Nfull, "obs_start:", obs_start
            for obnum in xrange(obs_start,Nfull):
                # Set up the empty H matrix
                H = np.zeros((1,Nfull))
                if obnum % 200 == 0 and rank == 2:
                    print "   On ob:", obnum-obs_start
                # Set up H for the ob --- simple
                H[0,obnum] = 1
                obval = state_meta[obnum][4]
                oberr = state_meta[obnum][5]
                oby = state_meta[obnum][0]
                obx = state_meta[obnum][1]
                # Figure out localization
                obloc = get_chunk_ideal_loc(yvals,xvals,oby,obx,loc,dx,maxx=maxx,
                                    maxy=maxy, periodic=True)
                # Now update the state
                xbm, Xbp = enkf_update(xbm, Xbp, np.array([[obval]]), H, oberr,\
                                       obloc, 1.0)
            # When done, return the state
            comm.send((start_index,Nstate,xbm[:Nstate],Xbp[:Nstate]), dest=0, tag=tags.DONE)

            # Check if we are within 3 stds from mean
            #Ye = np.dot(H, np.add(xbm[:, None], Xbp))
            #mye = np.mean(Ye,1)
            #sye = np.std(Ye,1)
            #if (obval >= (mye - (3*sye))) and (obval <= (mye + (3*sye))):
            #if True:
            #    xam, Xap = enkf_update(xbm, Xbp, np.array([[obval]]), H, oberr,\
            #                           obloc, 1.0)
            #else:
            #    # Return the original state if we aren't assimilating this ob
            #    xam = xbm
            #    Xap = Xbp
        elif tag == tags.EXIT:
            break
    comm.send(None, dest=0, tag=tags.EXIT)

    return None









def enkf_update(xbm, Xbp, Y, H, R, loc, inflate):
    # Originator: G. J. Hakim

    # Input variables are
    # Xbp => Array of ensemble estimates of state perturbations from mean (num_state_vars x num_ens_mems)
    # xbm => ensemble mean (vector of length num_state_vars)
    # Y => Observations (num_observations x num_state_vars)
    # H => Translation matrix from state vector to obs
    # R => Observation error covariance matrix
    # loc => Localization
    #print "Y shape:", np.shape(Y)
    #print "H shape:", np.shape(H)
    #print "R shape:", np.shape(R)
    #print "loc shape:", np.shape(loc)
    #raw_input()

    
    Nens = np.shape(Xbp)[1]   # Number of ensemble members
    Ndim = np.shape(Xbp)[0]   # Number of variables in state vector
    Nobs = np.shape(Y)[0]     # Total number of observations
    
    # Make a big vector of all of the ensemble members
    # estimates of state translated into observation
    # space (by H)
    #print "Mean", np.tile(xbm,(Nens,1))
    #print np.transpose(np.tile(xbm,(Nens,1))) + Xbp
    #print H
    #Ye = np.dot(H,np.transpose(np.tile(xbm,(Nens,1))) + Xbp)
    Ye = np.dot(H,np.add(xbm[:,None],Xbp))
    #For now, obs are just the first value
    #Ye = np.tile(xbm[0],Nens) + Xbp[0,:]
    #print "After:", Ye
    #raw_input()

    # The ensemble mean of the model estimate (in obs space)
    mye = np.mean(Ye,1)

    # Now loop over all observations
    for ob in xrange(Nobs):
        # Remove the mean from the model estimate of ob
        ye = Ye[ob,:] - mye[ob]
        #print "ye", ye
        # Find the variance among the ensemble members
        varye = np.var(ye)
        # And find the observation error variance from the R matrix
        obs_err = R[ob,ob]

        # Find the innovation --the difference between the ob value
        # and the ensemble mean estimate of the ob
        innov = Y[ob] - mye[ob]

        # Now find the innovation variance -- the sum of the variance of the ob
        # and the varaiance of the ensemble estimate
        # This goes into the denominator of the Kalman gain
        kdenom = (varye + obs_err)

        # The numerator of the Kalman gain is the covariance between
        # the ensemble members and the obs-transformed ensemble members
        kcov = np.dot(Xbp,np.transpose(ye)) / (Nens-1)

        # Option to inflate the covariances by a certain factor
        if inflate != None:
            kcov = inflate * kcov

        # Option to localize the gain
        if loc != None:
            kcov = np.multiply(kcov,np.transpose(loc[ob,:]))
            #kcov = np.dot(kcov,np.transpose(loc[ob,:]))
   
        # Compute the Kalman gain
        kmat = np.divide(kcov,kdenom)
        #print "kmat", kmat.shape
        #print "innov", innov.shape

        # Now do the updates
        # First update the mean
        #xam = xbm + np.dot(np.dot(H,kmat),innov)
        #print "kmat", np.shape(kmat)
        #print "innov", np.shape(innov)
        #xam = xbm + np.dot(kmat,innov)
        xam = xbm + np.multiply(kmat,innov)

        # And each ensemble member perturbation
        # this is the step in the Kalman filter equations
        beta = 1./(1. + np.sqrt(obs_err/(varye+obs_err)))
        kmat = np.multiply(beta,kmat)

        ye = np.array(ye)[np.newaxis]
        kmat = np.array(kmat)[np.newaxis]

        Xap = Xbp - np.dot(kmat.T, ye)

    # Return the analysis mean and perturbations
    return xam, Xap

def get_chunk_ideal_loc(statey,statex,oby,obx,cov_half,dx,\
                        maxx=None, maxy=None, periodic=True,use_boxcar=False):
    """ Function to compute localization assuming a rectangual grid (dx) """
    if maxx is None:
        maxx = np.max(statex)
        maxy = np.max(statey)

    # Compute the distance from the point to all state points
    alldist = dx * np.sqrt(np.add(np.power(np.subtract(statey,oby),2),\
                               np.power(np.subtract(statex,obx),2)))
    if periodic:
        for pt in [np.array([oby+maxy,obx]),
                  np.array([oby-maxy,obx]), np.array([oby,obx+maxx]),
                  np.array([oby,obx-maxx]), np.array([oby+maxy, obx+maxx]),
                  np.array([oby+maxy, obx-maxx]),
                  np.array([oby-maxy,obx+maxx]),np.array([oby-maxy, obx-maxx])]:
            # Compare to distance and use the minimum
            dist1 = np.sqrt(np.add(np.power(np.subtract(statey,pt[0]),2),np.power(np.subtract(statex,pt[1]),2)))
            alldist = np.minimum(alldist,dist1)
    # Now compute the localization
    dist = np.divide(alldist,cov_half)
    dist = np.where((dist == 0.0), 0.00001, dist)
    #covfact = np.zeros(statey.shape)

    if use_boxcar:
        # Simply make the weight one within twice the half-width
        covfact = np.where((dist <= 2), 1.0, 0.0)
    else:
        # Gaspari cohn
        covfact = np.where(dist <= 1, (((-0.25*dist+0.5)*dist+0.625)*dist-(5.0/3.0)) * np.power(dist,2) + 1,\
            ((((dist/12.0-0.5)*dist+0.625)*dist + (5.0/3.0)) * dist - 5.0) * dist + 4.0 - (2.0/(3.0*dist)))
        # Now apply the cutoff
        covfact = np.where((dist <= 2),covfact,0.0)

    # Return the covweight in the right format
    covweight = np.ones((1,statey.shape[0]))
    covweight[0,:] = covfact
    return covweight











        
if __name__ == '__main__':
    # Set up the MPI comm
    comm = MPI.COMM_WORLD # Get communicator
    size = comm.size # How many processors we have
    rank = comm.rank # Which processor is this
    status = MPI.Status() # Status object
    # Load the observations
    if rank == 0:
        try:
            print "Loading observations"
            obfile = open('observations.pckl','rb')
            obs = cPickle.load(obfile)
            obfile.close()
        except:
            print "ERROR: Unable to load observations file observations.pckl"
            exit(1)
        # Try to load the mean and perturbations
        try:
            print "Loading prior state"
            ensfile = open('prior_state.pckl','rb')
            xbm, Xbp = cPickle.load(ensfile)
            ensfile.close()
        except:
            print "ERROR: Unable to load ensemble prior state prior_state.pckl"
            exit(1)

        # Try to load the ensemble metadata 
        try:
            print "Loading ensemble meta"
            metafile = open('ensmeta.pckl','rb')
            ensmetas = cPickle.load(metafile)
            metafile.close()
        except:
            print "ERROR: Unable to load ensemble metadata ensmeta.pckl"
            exit(1)
        xam, Xap = mpi_filter(obs, ensmetas, xbm, Xbp)
    else:
        dum1, dum2 = mpi_filter(None, None, None, None)
