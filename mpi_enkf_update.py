#!/usr/bin/env python

from mpi4py import MPI
import numpy as np
from scipy import linalg
import cPickle
from util_wrap import stat_sig_local, cyth_efficient_localization

# Set up the MPI comm
comm = MPI.COMM_WORLD # Get communicator
size = comm.size # How many processors we have
rank = comm.rank # Which processor is this
status = MPI.Status() # Status object

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), *named)
    return type('Enum', (), enums)
tags = enum('READY','DONE','EXIT','START','AWAITSTATE','UPDATESTATE')



def main():
    # Rank 0 process is collecting and distributing
    if rank == 0:
        master(rank)
    else:
        worker(rank)


def master(rank=0):
    """ Function performed by the master node """
    # Load the observations
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
    # Now start figuring out the tasks
    obnames = obs.keys()
    tasks = len(obnames)
    task_index = 0
    num_workers = size - 1 # Rank 0 process not included as a worker
    closed_workers = 0
    state_hold = False
    awaiting_state = []
    print "Beginning MPI filter with %d workers" % num_workers
    while closed_workers < num_workers:
        # Get the state from a worker
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG,
                         status=status)
        source = status.Get_source()
        tag = status.Get_tag()

        if tag == tags.DONE:
            # The worker has processed the localization and H for an ob and is
            # ready for a new one
            # Update the state
            if state_hold == False:
                #comm.send((xbm,Xbp), dest=source, tag=tags.UPDATESTATE)
                state_hold = True
                # Take loc and H and update the state
                obloc, H, obval, oberr = data
                # Check if we are within 3 stds from mean
                Ye = np.dot(H, np.add(xbm[:, None], Xbp))
                mye = np.mean(Ye,1)
                sye = np.std(Ye,1)
                #if (obval >= (mye - (3*sye))) and (obval <= (mye + (3*sye))):
                if True:
                    xbm, Xbp = enkf_update(xbm, Xbp, np.array([[obval]]), H, oberr,\
                                           obloc, 1.0)
                state_hold = False
            #else:
            #    awaiting_state.append(source)
        elif tag == tags.READY:
            # Worker is ready having finished a task, so send it a new ob
            if task_index < tasks:
                #print "Sending ob to %d" % source
                if task_index % 50 == 0:
                    print "On ob number:", task_index
                obname = obnames[task_index]
                #print "sending ob:", obname
                obd = obs[obname]
                obd['obtime'] = obname[1]
                obd['obtype'] = obname[2]
                comm.send(obd, dest=source, tag=tags.START)
                task_index += 1
            else:
                # Out of obs, so send the worker a signal to shut down
                comm.send(None, dest=source, tag=tags.EXIT)
        elif tag == tags.EXIT:
            #print("Worker %d exited." % source)
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
    outfile = open('posterior_state.pckl','wb')
    cPickle.dump((xbm,Xbp), outfile)
    outfile.close()
    print "Done."

def worker(rank):
    """ This function defines what the worker does """
    name = MPI.Get_processor_name()
    print "Worker %d on %s is running" % (rank, name)
    # Get the broadcasted ensemble metadata
    try:
        #print "Loading ensemble meta"
        metafile = open('ensmeta.pckl','rb')
        ensmetas = cPickle.load(metafile)
        metafile.close()
    except:
        print "ERROR: Unable to load ensemble metadata ensmeta.pckl"
        exit(1)


    # Broadcast ensmeta to all workers
    statemeta = ensmetas[0]
    modlats = ensmetas[1]
    modlons = ensmetas[2]
    state_equiv = ensmetas[3]
    localization=ensmetas[4]
    metadict = dict(zip(statemeta, range(len(statemeta))))
    stateys = np.array([x[0] for x in statemeta])
    statexs = np.array([x[1] for x in statemeta])
    lonlist,latlist = np.meshgrid(modlons, modlats)
    #metadict = {}
    #state_equiv = {}
    #stateys = []
    #statexs = []
    #latlist = []
    #lonlist = []
    #localization=0.0
    Nstate = len(metadict.keys())

    while True:
        # Start by sending a ready tag
        comm.send(None, dest=0, tag=tags.READY)
        #print "Worker %d is ready" % rank
        #print "Worker %d has localization %f" % (rank, localization)
        obd = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()

        # Options here
        if tag == tags.START:
            #print "Worker %d has ob" % rank
            # Compute H for this ob
            H = np.zeros((1,Nstate))
            weight_dict = obd['weights']
            weight_sum = np.sum([1./weight_dict[k] for k in weight_dict.keys()])
            for k in weight_dict.keys():
                statevar = (k[0], k[1], obd['zgrid'], obd['obtime'],
                            state_equiv[obd['obtype']])
                state_loc = metadict[statevar]
                #print "   ",statevar, state_loc
                H[0,state_loc] = (1./weight_dict[k])/weight_sum
            outloc = cyth_efficient_localization(Nstate, stateys, statexs,
                                                 latlist, lonlist, obd['lat'],
                                                 obd['lon'], localization, 0)
            obloc = np.ones((1, Nstate))
            obloc[0,:] = np.array(outloc)

            oberr = obd['oberr']
            obval = obd['value']



            # Tell the master node we are ready for the state
            #print "Worker %d awaits state" % rank
            comm.send((obloc,H,obval, oberr), dest=0, tag=tags.DONE)
            """
            # Get the state when it is sent
            xbm, Xbp = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
            #print "Worker %d HAS state" % rank
            # Update the state
            # Check if we are within 3 stds from mean
            Ye = np.dot(H, np.add(xbm[:, None], Xbp))
            mye = np.mean(Ye,1)
            sye = np.std(Ye,1)
            #if (obval >= (mye - (3*sye))) and (obval <= (mye + (3*sye))):
            if True:
                xam, Xap = enkf_update(xbm, Xbp, np.array([[obval]]), H, oberr,\
                                       obloc, 1.0)
            #else:
            #    # Return the original state if we aren't assimilating this ob
            #    xam = xbm
            #    Xap = Xbp
            # Send the new state
            comm.send((xam,Xap), dest=0, tag=tags.DONE)
            """
        elif tag == tags.EXIT:
            break
    comm.send(None, dest=0, tag=tags.EXIT)










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
    Ye = np.dot(H,np.transpose(np.tile(xbm,(Nens,1))) + Xbp)
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

        
if __name__ == '__main__':
    main()
