#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import numpy as np
from pylab import *
import Sybase
import re
import os
import time

taskdata = os.path.join(os.environ.get('SKA', '/proj/sot/ska'),
                       'data', 'fid_drift_mon')

def get_fid_stats(cur, det):
    cur.execute('select id_num, id_string, tstart, ang_y_med, ang_z_med, sim_z_offset'
                ' FROM fid_stats'
                ' WHERE proc_status IS NULL AND id_string LIKE @det'
                , {'@det': det + '%'})
    vals = cur.fetchall()
    return np.rec.fromrecords(vals, names=[x[0] for x in cur.description])

def plotfids(detstats, det):
    figure(1, figsize=(6,4.5))
    clf()
    subplot(2, 1, 1)
    ylabel('Y offset (arcsec)')
    maxyear = int(time.strftime('%Y')) + 1
    axis([1999, maxyear, -25, 10])
    title(det + ' Fid Drift')
    grid()
    subplot(2, 1, 2)
    xlabel('Year')  # MET (years)')
    ylabel('Z offset (arcsec)')
    axis([1999, maxyear, -25, 10])
    grid()
    for fid in fids[det]:
        print 'Fid %s %d' % (det, fid)
        fidstats = detstats[detstats['id_num'] == fid]
        year = fidstats['tstart'] / 86400. / 365.25 + 1998.0
        normmask = np.logical_and(year > 2002.0, year < 2003.)
        fidstatsnorm = fidstats[normmask]
        if len(fidstatsnorm) < 3:
            print 'Not enough readouts for the median'
            continue
        y0 = np.median(fidstatsnorm['ang_y_med'])
        z0 = np.median(fidstatsnorm['ang_z_med'] + fidstatsnorm.sim_z_offset - sim_z_nom)
        subplot(2, 1, 1)
        plot(year-year0, fidstats['ang_y_med'] - y0, 
             ',', markerfacecolor=fidcolor[fid],
             scaley=False, scalex=False)
        subplot(2, 1, 2)
        plot(year-year0, fidstats['ang_z_med'] - z0 + fidstats['sim_z_offset'] - sim_z_nom,
             ',', markerfacecolor=fidcolor[fid],
             scaley=False, scalex=False)
    savefig(os.path.join(taskdata, 'drift_%s.png' % re.sub(r'-', '_', det.lower())))

##colormap = [!col.red,!col.magenta,!col.green,!col.black,!col.blue,!col.light_red]
year0 = 49075200 / 86400. / 365.25 + 1998.0  # year0 is year of launch Jul-23-1999
year0 = 0.
dets = ('ACIS-S', 'ACIS-I', 'HRC-S', 'HRC-I')
fids = {'ACIS-S': [1,2,3,4,5,6],
        'ACIS-I': [1,2,3,4,5,6],
        'HRC-I': [7,8,9,10],
        'HRC-S': [11,12,13,14],
        }
fidcolor = ' bgrcmkbgrcbgrc'

db = Sybase.connect('SYBASE', 'aca_ops', 'password_here', 'aca')
c = db.cursor()

for det in dets:
    # Some filtering here?
    detstats = get_fid_stats(c, det)
    sim_z_nom = np.median(detstats['sim_z_offset']) # median for this detector
    plotfids(detstats, det)

db.close()

