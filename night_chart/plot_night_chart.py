#!/usr/bin/env python3
import datetime
import numpy as np
from astropy.coordinates import SkyCoord
import ephem as ep
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.colors as mcolors


def plot_nightchart(longitude, latitude, elevation, correct_ut, xlabel,
    figname):

    # set observing place
    myplace = ep.Observer()
    myplace.long = longitude
    myplace.lat  = latitude
    myplace.elevation = elevation

    # set objects
    sun  = ep.Sun()
    moon = ep.Moon()

    # set range of plot dates
    date0 = datetime.datetime(2019,1,1,0,0,0)
    datetime_lst = [date0 + datetime.timedelta(days=i) for i in range(365)]
    datenums = mdates.date2num(datetime_lst)

    # initialize figure
    fig = plt.figure(dpi=150,figsize=(10,14))
    ax = fig.add_axes([0.12, 0.08, 0.82, 0.86])

    # set x range
    x1 = date0 - datetime.timedelta(hours=9)
    x2 = date0 + datetime.timedelta(hours=9)

    # plot sun rise and sun set time
    myplace.horizon=0
    dt1_lst, dt2_lst = [], []
    for idt, dt in enumerate(datetime_lst):
        myplace.date = dt
        t1 = myplace.previous_setting(sun, use_center=True)
        t2 = myplace.next_rising(sun, use_center=True)
        ut_offset1 = correct_ut(t1.datetime())
        ut_offset2 = correct_ut(t2.datetime())
        dt1 = t1.datetime() + ut_offset1 - dt
        dt2 = t2.datetime() + ut_offset2 - dt
        dt1_lst.append(date0 + dt1)
        dt2_lst.append(date0 + dt2)
    ax.plot(dt1_lst, datenums, '-', color='C3', alpha=0.5)
    ax.plot(dt2_lst, datenums, '-', color='C3', alpha=0.5)

    # plot civilized twilight time
    myplace.horizon='-6'
    dt1_lst, dt2_lst = [], []
    for idt, dt in enumerate(datetime_lst):
        myplace.date = dt
        t1 = myplace.previous_setting(sun, use_center=True)
        t2 = myplace.next_rising(sun, use_center=True)
        ut_offset1 = correct_ut(t1.datetime())
        ut_offset2 = correct_ut(t2.datetime())
        dt1 = t1.datetime() + ut_offset1 - dt
        dt2 = t2.datetime() + ut_offset2 - dt
        dt1_lst.append(date0 + dt1)
        dt2_lst.append(date0 + dt2)
    ax.fill_betweenx(datenums, dt1_lst, dt2_lst, color='C0', alpha=0.2, lw=0)

    # plot astronomical twilight time
    myplace.horizon='-18'
    dt1_lst, dt2_lst = [], []
    for idt, dt in enumerate(datetime_lst):
        myplace.date = dt
        t1 = myplace.previous_setting(sun, use_center=True)
        t2 = myplace.next_rising(sun, use_center=True)
        ut_offset1 = correct_ut(t1.datetime())
        ut_offset2 = correct_ut(t2.datetime())
        dt1 = t1.datetime() + ut_offset1 - dt
        dt2 = t2.datetime() + ut_offset2 - dt
        dt1_lst.append(date0 + dt1)
        dt2_lst.append(date0 + dt2)
    ax.fill_betweenx(datenums, dt1_lst, dt2_lst, color='C0', alpha=0.2, lw=0)

    # get times of full moons
    fullmoon_lst = []
    i = 0
    while(True):
        dt = date0 + datetime.timedelta(days=i)
        myplace.date = dt
        moon.compute(myplace)
        if dt > datetime_lst[-1]:
            break
        t = ep.next_full_moon(dt)
        found = False
        for _t in fullmoon_lst:
            deltat = _t - t.datetime()
            if abs(deltat.total_seconds())<24*3600:
                found = True
        if not found and t.datetime()<datetime_lst[-1]:
            fullmoon_lst.append(t.datetime())
        i += 15

    # get transit time of each full moon
    transit_lst = []
    for t in fullmoon_lst:
        myplace.date = t + datetime.timedelta(hours=5)
        transit_dt = myplace.previous_transit(moon).datetime()
        ut_offset = correct_ut(transit_dt)
        delta_dt = transit_dt + ut_offset - datetime.datetime(transit_dt.year,
                                                  transit_dt.month,
                                                  transit_dt.day, 0,0,0)
        if delta_dt.total_seconds()>12*3600:
            delta_dt = delta_dt - datetime.timedelta(hours=24)

        transit_lst.append(date0 + delta_dt)

    # plot the full moons
    ax.scatter(transit_lst, fullmoon_lst, color='C1', alpha=0.6, s=30, lw=0)

    # plot RA lines
    for ra in np.arange(0, 24, 3):
        # generate a fake star with (RA, Dec) = (ra, 0.0)
        eq = SkyCoord(ra=ra*15, dec=0., frame='icrs', unit='deg')
        string = 'NONAME,f|M|g2,%s,%s,5.0,2000'%(
                eq.ra.to_string(unit='h', sep=':'),
                eq.dec.to_string(unit='deg',sep=':'))
        star = ep.readdb(string)
        transit_lst = []
        real_transit_lst = []
        dt_lst = []

        transit_data = []
        prev_transit = None
        for dt in datetime_lst:
            myplace.date = dt + datetime.timedelta(hours=12)
            star.compute(myplace)

            real_transit_dt = myplace.previous_transit(star).datetime()
            delta_dt = real_transit_dt - dt
            transit = date0 + delta_dt

            if prev_transit is None or \
                abs((transit - prev_transit).total_seconds()) < 12*3600:
                transit_lst.append(transit)
                real_transit_lst.append(real_transit_dt)
                dt_lst.append(dt)
            else:
                item = (transit_lst, dt_lst, real_transit_lst)
                transit_data.append(item)
                # reset the list
                transit_lst = []
                real_transit_lst = []
                dt_lst = []
                prev_transit = None
            prev_transit = transit
        if len(dt_lst)>0:
            item = (transit_lst, dt_lst, real_transit_lst)
            transit_data.append(item)

        # plot
        for transit_lst, dt_lst, real_transit_lst in transit_data:
            ut_offset_lst = [correct_ut(dt) for dt in real_transit_lst]
            newtransit_lst = [t1+to for t1, to in zip(transit_lst, ut_offset_lst)]
            ax.plot(newtransit_lst, dt_lst, '--', color='k', alpha=0.2, lw=1)

            # add texts of RA in hour on the top
            if dt_lst[0] == datetime_lst[0] and x1 < newtransit_lst[0] < x2:
                ax.text(newtransit_lst[0],
                        dt_lst[0] + datetime.timedelta(days=6),
                        '%02dh'%ra, color='k', alpha=0.6)

            # add texts of RA in hour on the bottom
            if dt_lst[-1] == datetime_lst[-1] and x1 < newtransit_lst[-1] < x2:
                ax.text(newtransit_lst[-1],
                        dt_lst[-1] - datetime.timedelta(days=2),
                        '%02dh'%ra, color='k', alpha=0.6)

            # add texts of RA in hour on the left
            dist_lst = np.array([abs((x1-t).total_seconds()) for t in newtransit_lst])
            if dist_lst.min() < 2*3600:
                i = dist_lst.argmin()
                if i>0 and i<len(dist_lst):
                    ax.text(x1+datetime.timedelta(minutes=8), dt_lst[i],
                        '%02dh'%ra, color='k', alpha=0.6)

            # add texts of Dec in hour on the right
            dist_lst = np.array([abs((x2-t).total_seconds()) for t in newtransit_lst])
            if dist_lst.min() < 2*3600:
                i = dist_lst.argmin()
                if i>0 and i<len(dist_lst):
                    ax.text(x2-datetime.timedelta(minutes=30), dt_lst[i],
                        '%02dh'%ra, color='k', alpha=0.6)


    # set x and y axis to date
    ax.xaxis_date()
    ax.yaxis_date()

    # set ticks
    ax.xaxis.set_major_locator(mdates.HourLocator())
    ax.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=np.arange(0,10,60)))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    ax.yaxis.set_major_locator(mdates.MonthLocator())
    ax.yaxis.set_minor_locator(mdates.DayLocator(bymonthday=[1, 11, 21]))
    ax.yaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))

    # set x and y ranges
    ax.set_xlim(x1, x2)
    ax.set_ylim(datenums[-1], datenums[0])

    ax.grid(True, ls='--', alpha=0.5, which='both')
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel('Date', fontsize=14)

    fig.savefig(figname)

def plot_wst():

    def correct_ut(dt):
        # correct UT to Central European Time in 2019
        # returns +2 or +1 (CET - UT)
        # so CET = CET + correct_ut(dt)
        if dt > datetime.datetime(2019, 3, 31, 2, 0, 0) and \
           dt < datetime.datetime(2019, 10,27, 3, 0, 0):
            return datetime.timedelta(hours=2)
        else:
            return datetime.timedelta(hours=1)

    longitude = '12:00:43.4'
    latitude  = '47:42:13.1'
    elevation = 1950
    xlabel = 'Central European Time (CET)'
    figname = 'night_chart_wendelstein.pdf'
    plot_nightchart(longitude, latitude, elevation, correct_ut, xlabel, figname)

def plot_xinglong():

    def correct_ut(dt):
        # correct UT to Beijing Time
        # returns +8 (BJ - UT)
        return datetime.timedelta(hours=8)

    longitude = '117:34:28.35'
    latitude  = '40:23:45.36'
    elevation = 900
    xlabel = 'Beijing Time (UTC + 8)'
    figname = 'night_chart_xinglong.png'
    plot_nightchart(longitude, latitude, elevation, correct_ut, xlabel, figname)

if __name__=='__main__':
    plot_wst()
    plot_xinglong()
