#!/usr/bin/env python

from datetime import datetime, timedelta
import os

# Make a list of dates
starttimes = range(70,105,5)
os.system('touch RUN_IN_PROGRESS')

# Now loop through each and run
for start in starttimes:
    if not os.path.exists('RUN_IN_PROGRESS'):
        print "ABORTING RUN DUE TO REMOVAL OF RUN_IN_PROGRESS"
        exit(1)
    # Rewrite the namelist
    print "Rewriting time", start
    oldnml = open('NAMELIST.py','r')
    newnml = open('new_NAMELIST.py','w')
    for line in oldnml:
        if line.startswith('starttime'):
            print >>newnml, 'starttime = %d' % (start)
        else:
            print >>newnml, line[:-1]
    oldnml.close()
    newnml.close()
    os.system('mv new_NAMELIST.py NAMELIST.py')
    # Finally, run
    os.system('./ens_assimilate.py')


