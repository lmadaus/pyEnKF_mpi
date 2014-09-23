#!/usr/bin/env python

import numpy as np
from scipy import linalg

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

        
        
