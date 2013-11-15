import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Chandra.Time import DateTime

from calc_abs_cel_pointing import calc_all_dets

def mission_stats_by_year():
    drs = []
    y = []
    for year in range(1999, int(DateTime().frac_year) + 1):
        print "Year {}".format(year)
        y.append(year)
        drs.append(calc_all_dets('{}:001:00:00:00.000'.format(year),
                                 '{}:001:00:00:00.000'.format(year + 1),
                                 plot=False))
    plt.plot(y, drs, '.', markersize=20, color='blue')
    plt.ylabel('Radius of fid light offset from median (arcsec)')
    plt.title('All detector 99th percentile fid light offset radius')
    plt.xlabel('Year')
    plt.grid()
