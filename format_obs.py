#!/usr/bin/env python




def compute_cov_factor(z_in,c):
    # Based on DART module cov_cutoff_mod.f90
    # Computes the Gaspari Cohn weight to this point
    # given the distance from the observation to this
    # point (z_in) and the covariance cutoff distance
    # reqested (c)

    z = abs(z_in)
    # If we're greater than twice the cutoff distance
    # away, than just return zero
    if z >= 2*c:
        cov_weight = 0.0
    elif z <= c:
        # Compute the fraction of the cutoff we are away
        r = z/c
        # Keep only low order terms of GC
        cov_weight = (((-0.25*r + 0.5)*r+0.625)*r-5./3.)*(r**2)+1
    else:
        r = z/c
        cov_weight = ((((r/12. - 0.5) * r + 0.625) * r + 5./3.) * r - 5.0) * r\
            + 4.0 - 2.0/(3.0 * r)
    return cov_weight


def format_obs_text(ensemble_start,firsttime,lasttime,obtypes,siteids,state_equiv):
    # Use the state metadata to format the observation operators
    # This particular function only parses surface obs
    # from the local /home/disk/data directory
    #from surface_parse_bufkit import obs_parser
    from scipy.io import netcdf
    import numpy as np
    from datetime import datetime, timedelta
    import os
    import cPickle
    #from get_coords import get_coords
    from wrftools import latlon_to_ij
    print "Loading geo data"
    # This parameter sets how far from the boundary to ignore obs
    ob_ignore_halo = 5


    # Get basic geographical data from the template file
    terrfile = netcdf.netcdf_file('terrain_template.nc','r')
    lat_array = terrfile.variables['XLAT'][0,:,:]
    lon_array = terrfile.variables['XLONG'][0,:,:]
    #lonarr,latarr = np.meshgrid(lonlist,latlist)
    minlat = np.min(lat_array,axis=None)
    maxlat = np.max(lat_array,axis=None)
    minlon = np.min(lon_array,axis=None)
    maxlon = np.max(lon_array,axis=None)
    print "minlon:", minlon
    print "maxlon:", maxlon
    maxy,maxx = np.shape(lat_array)


    # Load the elevation
    #terrfile = Dataset('terrain.nc','r')
    model_elev = terrfile.variables['HGT'][0,:,:]
    #terrfile.close()


    # Build a list of times
    #firsttime = firsttime + timedelta(hours=6)  # REMOVE THIS LATER
    curtime = firsttime
    timelist = []
    while curtime <= lasttime:
        timelist.append(curtime)
        # Change this for variable time frequency later! LEM
        curtime = curtime + timedelta(hours=12)
    obd = {}        

    # Figure out if we have surface, upper air or both obs
    obs_from_sfc = ['alt','psfc','mslp','uwnd','vwnd','precip','temp','dewp']
    obs_from_upa = ['h_500hPa']

    # Sort this out
    sfc_obs = [f for f in obtypes if f in obs_from_sfc]
    upa_obs = [f for f in obtypes if f in obs_from_upa]

    print "Processing obs"
    print "Using LOCAL TEXT obs"
    print sfc_obs
    print upa_obs
    # Process these separately
    if len(sfc_obs) > 0:
        print "Working on surface obs."
        for time in timelist:
            for file in ['saous/%s.saous' % time.strftime('%Y%m%d%H'), 'saocn/%s.saocn' % time.strftime('%Y%m%d%H'),\
                        'cman/%s.cman' % time.strftime('%Y%m%d%H'), 'buoyfxd/%s.buoyfxd' % time.strftime('%Y%m%d%H')]:
                filename = '/home/disk/data/surface/decoded/%s' % (file)
                print filename
                try:
                    obfil = open(filename,'r')
                except:
                    print "File not found!"
                    continue
                oblines = obfil.readlines()
                oblist = xrange(len(oblines))
                obfil.close()

                if file.endswith('buoyfxd'):
                    maritime = True
                else:
                    maritime = False
                #except:
                #    print 'ERROR: time missing'
                #    print siteid, time
                #    continue
                for obline,obnum in zip(oblines,oblist):
                    if obnum % 500 == 0:
                        print "Processing ob:", obnum

                    obvals = obline.split(',')
                    try:
                        siteid = obvals[0][1:-1]
                    except:
                        continue
                    if siteid == '':
                        continue
                    # Line here to exclude certain obs
                    if siteids != None:
                        # If we're only assimilating certain sites, then ignore other sites
                        if siteid.upper() not in siteids:
                            continue

                    # Check to see if we're in the domain
                    try:
                        sitelat = float(obvals[5])
                        sitelon = float(obvals[6])
                        siteelev = float(obvals[7])
                    except:
                        continue
                    # WRF wants negative longitudes for west
                    #if sitelon < 0:
                    #    sitelon = 360. + sitelon
                    if (sitelat > maxlat) or (sitelat < minlat) or (sitelon > maxlon) or (sitelon < minlon):
                        continue

                    # Check the observation time
                    obhour = int(obvals[3])
                    obmin = int(obvals[4])
                    # Date comes from the filename
                    # check if we have to go forward or back a day
                    if time.hour == 0 and obhour == 23:
                        obdate = time - timedelta(days=1)
                    elif time.hour == 23 and obhour == 0:
                        obdate = time + timedelta(days=1)
                    else:
                        obdate = time

                    obstime = datetime(obdate.year, obdate.month, obdate.day, obhour, obmin)
                    # Check that we're within 15 minutes of the requested time
                    if (abs((obstime-time).seconds) + abs((obstime-time).days*86400)) > 15*60:
                        continue

                    for obtype in sfc_obs:
                        # Format this ob's name
                        obname = (siteid,int(time.strftime('%Y%m%d%H')),obtype)
                        # Now get the ob value
                        # Try-except will exclude observations that
                        # don't have this ob type
                        try:
                            if obtype == 'temp':
                                #obval = curob.tmpc * 9./5 + 32
                                # Surface temperature in Kelvin
                                obval = ((float(obvals[13])-32.) * 5./9.) + 273.
                                oberr = 1. 
                            elif obtype == 'wspd':
                                #obval = curob.sknt * 0.51444
                                # Wind speed in meters per second
                                obval = float(obvals[16]) * 0.51444
                                oberr = 1. 
                            elif obtype == 'precip':
                                # Convert to mm from hundredths of an inch
                                obval = float(obvals[21]) * 25.4/100.
                                oberr = 2. 
                            elif obtype == 'dewp':
                                # Dewpoint in Kelvin
                                #obval = curob.dwpc * 9./5 + 32
                                obval = ((float(obvals[14])-32.) * 5./9.) + 273.
                                oberr = 2. 
                            elif obtype == 'alt':
                                # Altimeter in hPa
                                """
                                if maritime:
                                    psfc = float(obfil.variables['stationPress'][obnum])/100.
                                    obval = (((psfc-0.3)**0.190284) + (siteelev * 8.4228807E-5))**(1/0.190284)
                                    if np.isnan(obval):
                                        # Use sea-level pressure as last resort and set elevation to zero
                                        obval = float(obfil.variables['seaLevelPress'][obnum])/100.
                                        siteelev = 0.
                                else:
                                """
                                obvaltemp = obvals[18]
                                try:
                                    int(obvaltemp)
                                except:
                                    continue
                                if int(obvaltemp) > 500:
                                    obval = float('2'+obvaltemp) * 33.86/100.
                                else:
                                    obval = float('3'+obvaltemp) * 33.86/100.
                                oberr = 1.
                                # Sanity check after calculation
                                if obval < 800 or obval > 1100:
                                    continue
                            #elif obtype == 'psfc':
                            #    # Sea-level pressure in hPa
                            #    obval = float(obfil.variables['stationPress'][obnum])/100.
                            #    oberr = 1. 
                            #elif obtype == 'mslp':
                            #    # Sea-level pressure in hPa
                            #    obval = float(obfil.variables['seaLevelPress'][obnum])/100.
                            #    oberr = 1. 
                            elif obtype == 'uwnd':
                                wdir = float(obvals[15])
                                wspd = float(obvals[16]) * 0.51444
                                obval = -1 * wspd * np.sin(4.0 * np.arctan(1.0)/180. * wdir)
                                oberr = 1.5
                            elif obtype == 'vwnd':
                                wdir = float(obvals[15])
                                wspd = float(obvals[16]) * 0.51444
                                obval = -1 * wspd * np.cos(4.0 * np.arctan(1.0)/180. * wdir)
                                oberr = 1.5 

                            # Get the closest i,j point to ob 
                            obx,oby = latlon_to_ij(terrfile,sitelat,sitelon)
                            if obx == None or oby == None:
                                # Ob is not in domain.  Continue
                                continue
                            # Check to be sure that the surrounding points are also at least 
                            # ob_ignore_halo gridpoints into the domain so we can interpolate
                            if obx <= ob_ignore_halo or abs(maxx-obx) <= ob_ignore_halo or\
                                oby <= ob_ignore_halo or abs(maxy-oby) <= ob_ignore_halo:
                                continue



                            weight_dict = {}
                            terrain_fail = False
                            single_point = False
                            #for point in ((obx_plus,oby_plus),(obx_plus,oby_minus),(obx_minus,oby_plus),(obx_minus,oby_minus)):
                            # Do a 9-point stencil from this point and weight accordingly
                            pointlist = [(obx,oby),(obx-1,oby), (obx+1,oby),(obx,oby-1),(obx-1,oby-1),(obx-1,oby+1),(obx,oby+1),\
                                (obx-1,oby+1),(obx+1,oby-1)]
                            for point in pointlist:
                                # Compute distance to this point
                                # remembering array dimensions go j,i
                                pointlon = lon_array[point[1],point[0]]
                                pointlat = lat_array[point[1],point[0]]
                                dist = np.sqrt((sitelat-pointlat)**2 + (sitelon-pointlon)**2)
                                # If we just happen to fall right on a grid point, don't interpolate
                                if dist > 0.001 and not single_point:
                                    weight_dict[point] = dist
                                else:
                                    weight_dict = {}
                                    weight_dict[point] = 1.0
                                    single_point = True

                                # Terrain check here (currently 400m)
                                modelev = model_elev[point[1],point[0]]
                                if abs(siteelev - modelev) > 400. and not terrain_fail:
                                    # Add to list of rejected obs
                                    #try:
                                    #    rejectfile = open('rejected_obs.pickle','r')
                                    #    rejects = cPickle.load(rejectfile)
                                    #    rejectfile.close()
                                    #except:
                                    #    rejects = []
                                    #rejects.append(obname[0] + ' ' +  str(obval) + ' ' + str(modelev) + ' ' + str(siteelev) + ' TERRAIN')
                                    terrain_fail = True

                            if terrain_fail:
                                #rejectfile = open('rejected_obs.pickle','w')
                                #rejects = cPickle.dump(rejects,rejectfile)
                                #rejectfile.close()
                                continue


                            obz = 0

                            # Now write the ob
                            if np.isnan(obval):
                                #print "Badob!"
                                continue
                            obd[obname] = {}
                            obd[obname]['value'] = obval
                            obd[obname]['oberr'] = oberr * np.eye(2,2)
                            obd[obname]['lat'] = sitelat
                            obd[obname]['lon'] = sitelon
                            obd[obname]['xgrid'] = obx
                            obd[obname]['ygrid'] = oby
                            obd[obname]['zgrid'] = obz
                            obd[obname]['weights'] = weight_dict

                        except:
                            continue
                obfil.close()
                print "   Obs so far:", len(obd.keys())


    terrfile.close()
    return obd





def format_obs(ensemble_start,firsttime,lasttime,obtypes,siteids,state_equiv):
    # Use the state metadata to format the observation operators
    #from surface_parse_bufkit import obs_parser
    from netCDF4 import Dataset
    import numpy as np
    from datetime import datetime, timedelta
    import os
    import cPickle
    #from get_coords import get_coords
    from wrftools import latlon_to_ij
    from scipy.io import netcdf
    print "Loading geo data"
    # This parameter sets how far from the boundary to ignore obs
    ob_ignore_halo = 5
    # This sets the base time for computing the obs time
    epoch = datetime(1970,1,1,0)


    # Get basic geographical data from the template file
    terrfile = netcdf.netcdf_file('terrain_template.nc','r')
    lat_array = terrfile.variables['XLAT'][0,:,:]
    lon_array = terrfile.variables['XLONG'][0,:,:]
    #lonarr,latarr = np.meshgrid(lonlist,latlist)
    minlat = np.min(lat_array,axis=None)
    maxlat = np.max(lat_array,axis=None)
    minlon = np.min(lon_array,axis=None)
    maxlon = np.max(lon_array,axis=None)
    print "minlon:", minlon
    print "maxlon:", maxlon
    maxy,maxx = np.shape(lat_array)


    # Load the elevation
    #terrfile = Dataset('terrain.nc','r')
    model_elev = terrfile.variables['HGT'][0,:,:]
    #terrfile.close()


    # Build a list of times
    #firsttime = firsttime + timedelta(hours=6)  # REMOVE THIS LATER
    curtime = firsttime
    timelist = []
    while curtime <= lasttime:
        timelist.append(curtime)
        # Change this for variable time frequency later! LEM
        curtime = curtime + timedelta(hours=12)
    obd = {}        

    # Figure out if we have surface, upper air or both obs
    obs_from_sfc = ['alt','psfc','mslp','uwnd','vwnd','precip','temp','dewp']
    obs_from_upa = ['h_500hPa']

    # Sort this out
    sfc_obs = [f for f in obtypes if f in obs_from_sfc]
    upa_obs = [f for f in obtypes if f in obs_from_upa]

    print "Processing obs"
    print sfc_obs
    print upa_obs
    # Process these separately
    if len(sfc_obs) > 0:
        print "Working on surface obs."
        for time in timelist:
            for file in ['%s_maritime' % time.strftime('%Y%m%d_%H00'), '%s_metar' % time.strftime('%Y%m%d_%H00')]:
                filename = './obs/%s' % (file)
                print filename
                try:
                    obfil = Dataset(filename,'r')
                except:
                    print "File not found!"
                    continue
                oblist = xrange(len(obfil.dimensions['recNum']))
                if file.endswith('maritime'):
                    maritime = True
                else:
                    maritime = False
                #except:
                #    print 'ERROR: time missing'
                #    print siteid, time
                #    continue
                for obnum in oblist:
                    if obnum % 500 == 0:
                        print "Processing ob:", obnum
                    try:
                        siteid = ''.join(filter(None,obfil.variables['stationName'][obnum]))
                    except:
                        continue
                    if siteid == '':
                        continue
                    # Line here to exclude certain obs
                    if siteids != None:
                        # If we're only assimilating certain sites, then ignore other sites
                        if siteid.upper() not in siteids:
                            continue

                    # Check to see if we're in the domain
                    sitelat = float(obfil.variables['latitude'][obnum])
                    sitelon = float(obfil.variables['longitude'][obnum])
                    siteelev = float(obfil.variables['elevation'][obnum])
                    # WRF wants negative longitudes for west
                    #if sitelon < 0:
                    #    sitelon = 360. + sitelon
                    if (sitelat > maxlat) or (sitelat < minlat) or (sitelon > maxlon) or (sitelon < minlon):
                        continue
                    
                    # Check the observation time
                    obstime = epoch + timedelta(obfil.variables['timeObs'][obnum])
                    # Check that we're within 15 minutes of the requested time
                    if (abs((obstime-time).seconds) + abs((obstime-time).days*86400)) > 15*60:
                        continue

                    for obtype in sfc_obs:
                        # Format this ob's name
                        obname = (siteid,int(time.strftime('%Y%m%d%H')),obtype)
                        # Now get the ob value
                        if True:
                        #try:
                            if obtype == 'temp':
                                #obval = curob.tmpc * 9./5 + 32
                                # Surface temperature in Kelvin
                                obval = float(obfil.variables['temperature'][obnum])
                                oberr = 1. 
                            elif obtype == 'wspd':
                                #obval = curob.sknt * 0.51444
                                # Wind speed in meters per second
                                obval = float(obfil.variables['windSpeed'][obnum])
                                oberr = 1. 
                            elif obtype == 'precip':
                                obval = np.sqrt(float(obfil.variables['precip6Hour'][obnum]))
                                oberr = 0.1 
                            elif obtype == 'dewp':
                                # Dewpoint in Kelvin
                                #obval = curob.dwpc * 9./5 + 32
                                obval = float(obfil.variables['dewpoint'][obnum])
                                oberr = 2. 
                            elif obtype == 'alt':
                                # Altimeter in hPa
                                if maritime:
                                    psfc = float(obfil.variables['stationPress'][obnum])/100.
                                    obval = (((psfc-0.3)**0.190284) + (siteelev * 8.4228807E-5))**(1/0.190284)
                                    if np.isnan(obval):
                                        # Use sea-level pressure as last resort and set elevation to zero
                                        obval = float(obfil.variables['seaLevelPress'][obnum])/100.
                                        siteelev = 0.
                                else:
                                    obval = float(obfil.variables['altimeter'][obnum])/100.
                                oberr = 1.
                                # Sanity check after calculation
                                if obval < 800 or obval > 1100:
                                    continue
                            elif obtype == 'psfc':
                                # Sea-level pressure in hPa
                                obval = float(obfil.variables['stationPress'][obnum])/100.
                                oberr = 1. 
                            elif obtype == 'mslp':
                                # Sea-level pressure in hPa
                                obval = float(obfil.variables['seaLevelPress'][obnum])/100.
                                oberr = 1. 
                            elif obtype == 'uwnd':
                                wdir = float(obfil.variables['windDir'][obnum])
                                wspd = float(obfil.variables['windSpeed'][obnum])
                                obval = -1 * wspd * np.sin(4.0 * np.arctan(1.0)/180. * wdir)
                                oberr = 1.5
                            elif obtype == 'vwnd':
                                wdir = float(obfil.variables['windDir'][obnum])
                                wspd = float(obfil.variables['windSpeed'][obnum])
                                obval = -1 * wspd * np.cos(4.0 * np.arctan(1.0)/180. * wdir)
                                oberr = 1.5 

                            
                            """
                            THIS BLOCK TO FIGURE OUT WHAT THE I-J COORDINATES OF THE
                            OB LOCATION IN THE DOMAIN IS REPLACED BY EXTERNAL latlon_to_ij
                            FUNCTION SPECIFICALLY FOR WRF
                            # Sort out the H vector
                            #H = np.zeros((1,len(statemeta)))
                            #obx,oby = get_coords(latarr,lonarr,curob.lat,curob.lon)
                            #latdfound = False
                            #for l in latlist:
                            #    if sitelat < l and not latdfound:
                            #        oby_plus = latlist.index(l)
                            #        oby_minus = oby_plus-1
                            #        latdfound = True
                            #londfound = False
                            #for l in lonlist:
                            #    if sitelon < l and not londfound:
                            #        if sitelon < 180. and l >= 180.:
                            #            continue
                            #        obx_plus = lonlist.index(l)
                            #        obx_minus = obx_plus-1
                            #        londfound = True
                            #if sitelon < 180.:
                            #    print lonlist
                            #    print sitelon, "obx_plus:",obx_plus, "obx_minus:",obx_minus
                            #    raw_input()
                            # Now sort out the weighting for linear interpolation
                            # We have the four nearest coordinates
                            # Check for wrapping:
                            #if obx_minus < 0:
                            #    obx_minus = len(lonlist) + obx_minus
                            # Ignore ob if off domain
                            if obx_minus < 0 or obx_plus >= len(lonlist) or oby_minus < 0 or oby_plus >= len(latlist):
                                continue
                            """
                            obx,oby = latlon_to_ij(terrfile,sitelat,sitelon)
                            if obx == None or oby == None:
                                # Ob is not in domain.  Continue
                                continue
                            # Check to be sure that the surrounding points are also at least 
                            # ob_ignore_halo gridpoints into the domain so we can interpolate
                            if obx <= ob_ignore_halo or abs(maxx-obx) <= ob_ignore_halo or\
                                oby <= ob_ignore_halo or abs(maxy-oby) <= ob_ignore_halo:
                                continue



                            weight_dict = {}
                            terrain_fail = False
                            single_point = False
                            #for point in ((obx_plus,oby_plus),(obx_plus,oby_minus),(obx_minus,oby_plus),(obx_minus,oby_minus)):
                            # Do a 9-point stencil from this point and weight accordingly
                            pointlist = [(obx,oby),(obx-1,oby), (obx+1,oby),(obx,oby-1),(obx-1,oby-1),(obx-1,oby+1),(obx,oby+1),\
                                (obx-1,oby+1),(obx+1,oby-1)]
                            for point in pointlist:
                                # Compute distance to this point
                                # remembering array dimensions go j,i
                                pointlon = lon_array[point[1],point[0]]
                                pointlat = lat_array[point[1],point[0]]
                                dist = np.sqrt((sitelat-pointlat)**2 + (sitelon-pointlon)**2)
                                # If we happen to fall exactly on a grid point, don't interpolate
                                if dist > 0.001 and not single_point:
                                    weight_dict[point] = dist
                                else:
                                    weight_dict = {}
                                    weight_dict[point] = 1.0
                                    single_point = True
                                # Terrain check here (currently 400m)
                                modelev = model_elev[point[1],point[0]]
                                if abs(siteelev - modelev) > 400. and not terrain_fail:
                                    # Add to list of rejected obs
                                    #try:
                                    #    rejectfile = open('rejected_obs.pickle','r')
                                    #    rejects = cPickle.load(rejectfile)
                                    #    rejectfile.close()
                                    #except:
                                    #    rejects = []
                                    #rejects.append(obname[0] + ' ' +  str(obval) + ' ' + str(modelev) + ' ' + str(siteelev) + ' TERRAIN')
                                    terrain_fail = True

                            if terrain_fail:
                                #rejectfile = open('rejected_obs.pickle','w')
                                #rejects = cPickle.dump(rejects,rejectfile)
                                #rejectfile.close()
                                continue


                            obz = 0

                            # Now write the ob
                            if np.isnan(obval):
                                #print "Badob!"
                                continue
                            obd[obname] = {}
                            obd[obname]['value'] = obval
                            obd[obname]['oberr'] = oberr * np.eye(2,2)
                            obd[obname]['lat'] = sitelat
                            obd[obname]['lon'] = sitelon
                            obd[obname]['xgrid'] = obx
                            obd[obname]['ygrid'] = oby
                            obd[obname]['zgrid'] = obz
                            obd[obname]['weights'] = weight_dict

                        #except:
                        #    continue
                obfil.close()
                print "   Obs so far:", len(obd.keys())

    if len(upa_obs) > 0:
        print "Working on upper air obs."
        for time in timelist:
            filename = './obs/%s_upa' % (time.strftime('%Y%m%d_%H00'))
            print filename
            obfil = Dataset(filename,'r')
            #try:
            if True:
                oblist = xrange(len(obfil.dimensions['recNum']))

            #except:
            #    print 'ERROR: time missing'
            #    print siteid, time
            #    continue
            for obnum in oblist:
                if obnum % 200 == 0:
                    print "Processing ob:", obnum
                try:
                    siteid = ''.join(filter(None,obfil.variables['staName'][obnum]))
                except:
                    print "Site ID not usable"
                    continue
                if siteid == '':
                    print "Site ID is blank"
                    continue
                if siteids != None:
                    # If we're only assimilating certain sites, then ignore other sites
                    if siteid.upper() not in siteids:
                        continue

                # Check to see if we're in the domain
                sitelat = float(obfil.variables['staLat'][obnum])
                sitelon = float(obfil.variables['staLon'][obnum])
                if sitelon < 0:
                    sitelon = 360. + sitelon
                if (sitelat > maxlat) or (sitelat < minlat) or (sitelon > maxlon) or (sitelon < minlon):
                    print "SITE NOT IN BOX"
                    print maxlat, minlat, sitelat
                    print maxlon, minlon, sitelon
                    print "_____________________"
                    raw_input()
                    continue
                for obtype in upa_obs:
                    # Format this ob's name
                    obname = (siteid,int(time.strftime('%Y%m%d%H')),obtype)
                    # Now get the ob value
                    # Separate the obtype into the height and the press level
                    obvar,oblev = obtype.split('_')
                    oblev = float(oblev[:-3])
                    # Find where this level is
                    plevs = obfil.variables['prMan'][obnum,:]
                    plevs = list(plevs)
                    try:
                        # See if this level is in the profile
                        oblevdex = plevs.index(oblev)
                    except:
                        # Otherwise go on to the next one
                        print "Observation lacks this vertical level"
                        continue
                    #if True:

                    try:
                        if obvar in ('t','T'):
                            # Temperature in Kelvin
                            obval = float(obfil.variables['tpMan'][obnum,oblevdex])
                            oberr = 2. 
                        elif obtype in ('v','V'):
                            wdir = float(obfil.variables['wdMan'][obnum,oblevdex])
                            wspd = float(obfil.variables['wsMan'][obnum,oblevdex])
                            obval = -1 * wspd * np.cos(4.0 * np.arctan(1.0)/180. * wdir)
                            oberr = 5. 
                        elif obtype in ('U','u'):
                            wdir = float(obfil.variables['wdMan'][obnum,oblevdex])
                            wspd = float(obfil.variables['wsMan'][obnum,oblevdex])
                            obval = -1 * wspd * np.sin(4.0 * np.arctan(1.0)/180. * wdir)
                            oberr = 5. 
                        elif obvar in ('H','h'):
                            # Geopotential height
                            obval = float(obfil.variables['htMan'][obnum,oblevdex])
                            oberr = 10. 


                        # Sort out the H vector
                        #H = np.zeros((1,len(statemeta)))

                        #obx,oby = get_coords(latarr,lonarr,curob.lat,curob.lon)
                        latdfound = False
                        for l in latlist:
                            if sitelat < l and not latdfound:
                                oby_plus = latlist.index(l)
                                oby_minus = oby_plus-1
                                latdfound = True
                        londfound = False
                        for l in lonlist:
                            if sitelon < 180. and l >= 180.:
                                continue
                            if sitelon < l and not londfound:
                                obx_plus = lonlist.index(l)
                                obx_minus = obx_plus-1
                                londfound = True
                        # Now sort out the weighting for linear interpolation
                        # We have the four nearest coordinates
                        # Check for wrapping:
                        #if obx_minus < 0:
                        #    obx_minus = len(lonlist) + obx_minus
                        # Ignore ob if off domain
                        if obx_minus < 0 or obx_plus >= len(lonlist) or oby_minus < 0 or oby_plus >= len(latlist):

                            #print "Observation too close to edge of domain"
                            #print "sitelon:", sitelon, "startlon:", lonlist[0], "endlon:",lonlist[-1]
                            #print lonlist
                            #print "obxplus:", obx_plus, "obxminus:", obx_minus
                            continue

                        weight_dict = {}
                        for point in ((obx_plus,oby_plus),(obx_plus,oby_minus),(obx_minus,oby_plus),(obx_minus,oby_minus)):
                            # Compute distance to this point
                            pointlon = lonlist[point[0]]
                            pointlat = latlist[point[1]]
                            dist = np.sqrt((sitelat-pointlat)**2 + (sitelon-pointlon)**2)
                            weight_dict[point] = dist

                        obz = 0
                        # Now write the ob
                        if np.isnan(obval):
                            #print "Bad observation value!"
                            continue
                        obd[obname] = {}
                        obd[obname]['value'] = obval
                        obd[obname]['oberr'] = oberr * np.eye(2,2)
                        obd[obname]['lat'] = sitelat
                        obd[obname]['lon'] = sitelon
                        obd[obname]['xgrid'] = obx_plus
                        obd[obname]['ygrid'] = oby_plus
                        obd[obname]['zgrid'] = obz
                        obd[obname]['weights'] = weight_dict
                        if obnum % 200 == 0:
                            print obname, obval

                    except:
                        continue
            obfil.close()




    terrfile.close()
    return obd





