#!/usr/bin/env python


def format_ideal_obs(ensemble_start, firsttime, lasttime, obtypes, siteids,
                     state_equiv,modlats,modlons):
    """ Use the state metadata to format these observations """
    from netCDF4 import Dataset
    import numpy as np
    import os
    import cPickle
    from bisect import bisect

    # Get the file to derive truth from
    from NAMELIST import truthdat, truth_obsloc

    # Build a list of times
    obtimes = range(ensemble_start+firsttime,ensemble_start+lasttime+1)
    #obtimes = [ensemble_start + firsttime]

    obd = {}
    print "Processing observations"
    print "Deriving from TRUTH file:", truthdat
    for time in obtimes:
        try:
            obfil = Dataset(truthdat, 'r')
        except:
            print "Unable to open truth file for observations:", truthdat
            exit(1)

        # Set up the dimensions
        try:
            xcoords = obfil.variables['xh']
            ycoords = obfil.variables['yh']
        except:
            xcoords = obfil.variables['lon']
            ycoords = obfil.variables['lat']

        dy = ycoords[2]-ycoords[1]
        dx = xcoords[2]-xcoords[1]
        nx = xcoords.shape[0]
        ny = ycoords.shape[0]

        # Figure out how we are going to grab these observations
        if truth_obsloc.startswith('GRID'):
            truth_spacing = float(truth_obsloc[4:-2])
            # This is the grid spacing in kilometers.  Find out how many
            # gridpoints in this would be
            grdpt_space = int(truth_spacing/dx)+1
            xpts = range(0,nx,grdpt_space)[:-1]
            ypts = range(0,ny,grdpt_space)[:-1]
            # Make a list of the combinations
            truth_points = [(y,x) for y in ypts for x in xpts]
        else:
            # In the future, allow for specific specification of locations
            # For now, just exit as unrecognizable
            print "Unrecognized truth obs locations: trugh_obsloc=", truth_obsloc
            exit(1)


        #  Now loop through each ob type
        for obtype in obtypes:
            obfield = obfil.variables[state_equiv[obtype]]
            # Loop through all ob
            for pt in truth_points:
                oby,obx = pt
                if obtype in ['psfc']:
                #    obval = np.mean(obfield[time,1:3,oby,obx],dtype=np.float64)
                    obval /= 100.0
                obval = obfield[time,oby,obx]
                oberr = 1.0
                # IMPORTANT -- perturb the ob value in accordance with the ob
                # error variance
                obval = np.random.normal(loc=obval, scale=np.sqrt(oberr))

                # Get the "kilometer" coordinates
                obykm = ycoords[oby]
                obxkm = xcoords[obx]
                # Compute weighting
                single_point = False
                # Do a four-point stencil for weighting based on the ensemble
                # "lats" and "lons" provided

                right = bisect(modlons,obxkm)
                top = bisect(modlats,obykm)
                left = right - 1
                bottom = top - 1

                if left < 0:
                    left = left + len(modlons) 
                if bottom < 0:
                    bottom = bottom + len(modlats)

                weight_dict = {}
                single_point = False
                for point in ((bottom,left), (bottom,right),\
                              (top, right), (top, left)):
                    
                    ensykm = modlats[point[0]]
                    ensxkm = modlons[point[0]]
                    dist = np.sqrt((ensykm - obykm)**2 + (ensxkm-obxkm)**2)
                    if dist > 0.001 and not single_point:
                        weight_dict[point] = dist
                    else:
                        # If we're right on the point (should often be the case)
                        # then just include that point
                        weight_dict = {}
                        weight_dict[point] = 1.0
                        single_point = True
                obz = 0 # For now
                # Now write the ob
                obname = ('%3.1f/%3.1f' % (obykm, obxkm), time, obtype)
                obd[obname] = {}
                obd[obname]['value'] = obval
                obd[obname]['oberr'] = oberr * np.eye(2,2)
                obd[obname]['lon'] = obxkm
                obd[obname]['lat'] = obykm
                obd[obname]['xgrid'] = obx
                obd[obname]['ygrid'] = oby
                obd[obname]['zgrid'] = obz
                obd[obname]['weights'] = weight_dict
                #print obname
                #print obd[obname]
        obfil.close()
    return obd


                        


                

