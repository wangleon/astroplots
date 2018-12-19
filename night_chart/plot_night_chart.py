#!/usr/bin/env python3
import datetime
import numpy as np
from astropy.coordinates import SkyCoord
import ephem as ep
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.colors as mcolors

def calc_sun(myplace, horizon, local_datetime_lst, correct_ut, date0):

    sun  = ep.Sun()

    # plot sun rise and sun set time
    myplace.horizon=horizon
    y_lst = []
    dt1_lst, dt2_lst = [], []
    for idt, dt in enumerate(local_datetime_lst):
        ut_offset = correct_ut(dt)
        ut = dt - ut_offset
        myplace.date = ut
        t1 = myplace.previous_setting(sun, use_center=True)
        t2 = myplace.next_rising(sun, use_center=True)
        midnight = datetime.datetime.combine(dt.date(), datetime.time.min)
        dt1 = t1.datetime() + ut_offset - midnight
        dt2 = t2.datetime() + ut_offset - midnight

        if dt1 > datetime.timedelta(hours=12):
            dt1 -= datetime.timedelta(hours=24)
        elif dt1 < datetime.timedelta(hours=-12):
            dt1 += datetime.timedelta(hours=24)
        if dt2 > datetime.timedelta(hours=12):
            dt2 -= datetime.timedelta(hours=24)
        elif dt2 < datetime.timedelta(hours=-12):
            dt2 += datetime.timedelta(hours=24)

        # fix daylight saving time problem
        if len(dt1_lst) > 0:
            ddt = (date0 + dt1 - dt1_lst[-1]).total_seconds()
            if ddt > 50*60:
                y_lst.append(dt)
                dt1_lst.append(date0 + dt1 - datetime.timedelta(hours=1))
                dt2_lst.append(date0 + dt2 - datetime.timedelta(hours=1))
            if ddt < -50*60:
                y_lst.append(dt)
                dt1_lst.append(date0 + dt1 + datetime.timedelta(hours=1))
                dt2_lst.append(date0 + dt2 + datetime.timedelta(hours=1))

        y_lst.append(dt)
        dt1_lst.append(date0 + dt1)
        dt2_lst.append(date0 + dt2)
    return y_lst, dt1_lst, dt2_lst
    

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
    #local_datetime_lst = [date0 + datetime.timedelta(hours=2*i) for i in range(365*12)]
    local_datetime_lst = [date0 + datetime.timedelta(days=i) for i in range(365)]
    datenums = mdates.date2num(local_datetime_lst)

    # initialize figure
    fig = plt.figure(dpi=150,figsize=(10,14))
    ax = fig.add_axes([0.12, 0.08, 0.82, 0.86])

    # set x range
    x1 = date0 - datetime.timedelta(hours=9)
    x2 = date0 + datetime.timedelta(hours=9)

    plot_sun  = True
    plot_twg  = True
    plot_moon = True
    plot_ra   = True

    if plot_sun:
        # plot sun rise and sun set time
        y_lst, dt1_lst, dt2_lst = calc_sun(myplace, 0, local_datetime_lst, correct_ut, date0)
        artist_sun, = ax.plot(dt1_lst, mdates.date2num(y_lst), '-', color='C3', alpha=0.7)
        ax.plot(dt2_lst, mdates.date2num(y_lst), '-', color='C3', alpha=0.7)

    if plot_twg:
        # plot civilized twilight time
        y_lst, dt1_lst, dt2_lst = calc_sun(myplace, '-6', local_datetime_lst, correct_ut, date0)
        artist_ctwg = ax.fill_betweenx(mdates.date2num(y_lst), dt1_lst, dt2_lst, color='C0',
                alpha=0.2, lw=0)

        # plot astronomical twilight time
        y_lst, dt1_lst, dt2_lst = calc_sun(myplace, '-18', local_datetime_lst, correct_ut, date0)
        ax.fill_betweenx(mdates.date2num(y_lst), dt1_lst, dt2_lst, color='C0', alpha=0.2, lw=0)

        # plot a fake twilight artist
        artist_atwg = ax.fill_betweenx(mdates.date2num(y_lst)-1e3,dt1_lst, dt2_lst,
                color='C0', alpha=0.36, lw=0) 

    if plot_moon:
        # get times of full moons
        fullmoon_lst = []
        i = 0
        while(True):
            dt = date0 + datetime.timedelta(days=i)
            if dt > local_datetime_lst[-1]:
                break
            myplace.date = dt
            moon.compute(myplace)
            t = ep.next_full_moon(dt)
            found = False
            for _t in fullmoon_lst:
                deltat = _t - t.datetime()
                if abs(deltat.total_seconds())<24*3600:
                    found = True
            if not found and t.datetime()<local_datetime_lst[-1]:
                fullmoon_lst.append(t.datetime())
            i += 15
 
        # get transit time of each full moon
        transit_lst = []
        for t in fullmoon_lst:
            myplace.date = t + datetime.timedelta(hours=5)
            transit_dt = myplace.previous_transit(moon).datetime()
            ut_offset = correct_ut(transit_dt)
            midnight = datetime.datetime.combine(transit_dt.date(), datetime.time.min)
            delta_dt = transit_dt + ut_offset - midnight
            if delta_dt.total_seconds()>12*3600:
                delta_dt -= datetime.timedelta(hours=24)
 
            transit_lst.append(date0 + delta_dt)
 
        # plot the full moons
        artist_moon = ax.scatter(transit_lst, fullmoon_lst, color='C1', alpha=0.6,
                s=30, lw=0)

    if plot_ra:
        # plot RA lines
        for ra in np.arange(0, 24, 2):
            # generate a fake star with (RA, Dec) = (ra, 0.0)
            eq = SkyCoord(ra=ra*15, dec=0., frame='icrs', unit='deg')
            string = 'NONAME,f|M|g2,%s,%s,5.0,2000'%(
                    eq.ra.to_string(unit='h', sep=':'),
                    eq.dec.to_string(unit='deg',sep=':'))
            star = ep.readdb(string)
            transit_lst = []
            dt_lst = []
 
            transit_data = []
            prev_transit = None
            for dt in local_datetime_lst:
                midnight = datetime.datetime.combine(dt.date(), datetime.time.min)
                ut_offset = correct_ut(dt)
                ut = dt - ut_offset
                myplace.date = ut + datetime.timedelta(hours=12)
                # add 12 hours to aviod the transit time jump problem
                star.compute(myplace)
 
                real_transit_dt = myplace.previous_transit(star).datetime()
                delta_dt = real_transit_dt + ut_offset - midnight
                if delta_dt > datetime.timedelta(hours=12):
                    delta_dt -= datetime.timedelta(hours=24)
                elif delta_dt < datetime.timedelta(hours=-12):
                    delta_dt += datetime.timedelta(hours=24)
                transit = date0 + delta_dt
 
                if prev_transit is None or \
                    abs((transit - prev_transit).total_seconds()) < 12*3600:

                    # fix daylight saving time problem
                    if prev_transit is not None and \
                        (transit - prev_transit).total_seconds() > 50*60:
                        transit_lst.append(transit - datetime.timedelta(hours=1))
                        dt_lst.append(dt)
                    if prev_transit is not None and \
                        (transit - prev_transit).total_seconds() < -50*60:
                        transit_lst.append(transit + datetime.timedelta(hours=1))
                        dt_lst.append(dt)

                    transit_lst.append(transit)
                    dt_lst.append(dt)
                else:
                    item = (transit_lst, dt_lst)
                    transit_data.append(item)
                    # reset the list
                    transit_lst = []
                    dt_lst = []
                    prev_transit = None
                prev_transit = transit
            if len(dt_lst)>0:
                item = (transit_lst, dt_lst)
                transit_data.append(item)
 
            # plot RA transit lines (also local sidereal time)
            for transit_lst, dt_lst in transit_data:
 
                artist_ra, = ax.plot(transit_lst, dt_lst, '--', color='k', alpha=0.2, lw=1)
 
                # add texts of RA in hour on the top
                if dt_lst[0] == local_datetime_lst[0] and x1 < transit_lst[0] < x2:
                    ax.text(transit_lst[0] - datetime.timedelta(minutes=15),
                            dt_lst[0] + datetime.timedelta(days=6),
                            '%02dh'%ra, color='k', alpha=0.6)
 
                # add texts of RA in hour on the bottom
                if dt_lst[-1] == local_datetime_lst[-1] and x1 < transit_lst[-1] < x2:
                    ax.text(transit_lst[-1] - datetime.timedelta(minutes=15),
                            dt_lst[-1] - datetime.timedelta(days=2),
                            '%02dh'%ra, color='k', alpha=0.6)
 
                # add texts of RA in hour on the left
                dist_lst = np.array([abs((x1-t).total_seconds()) for t in transit_lst])
                if dist_lst.min() < 2*3600:
                    i = dist_lst.argmin()
                    if i>0 and i<len(dist_lst):
                        ax.text(x1+datetime.timedelta(minutes=8), dt_lst[i],
                            '%02dh'%ra, color='k', alpha=0.6)
 
                # add texts of Dec in hour on the right
                dist_lst = np.array([abs((x2-t).total_seconds()) for t in transit_lst])
                if dist_lst.min() < 2*3600:
                    i = dist_lst.argmin()
                    if i>0 and i<len(dist_lst):
                        ax.text(x2-datetime.timedelta(minutes=40), dt_lst[i],
                            '%02dh'%ra, color='k', alpha=0.6)


    # set x and y axis to date
    ax.xaxis_date()
    ax.yaxis_date()

    # set ticks
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=np.arange(0,24,2)))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=np.arange(0,24,1)))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    ax.yaxis.set_major_locator(mdates.MonthLocator())
    ax.yaxis.set_minor_locator(mdates.DayLocator(bymonthday=[1, 11, 21]))
    ax.yaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))

    # set x and y ranges
    ax.set_xlim(x1, x2)
    ax.set_ylim(datenums[-1], datenums[0])

    ax.grid(True, ls=':', alpha=0.6, lw=1,  which='major')
    ax.grid(True, ls=':', alpha=0.5, lw=0.5, which='minor')
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel('Date', fontsize=14)

    ax.legend((artist_sun, artist_moon, artist_ctwg, artist_atwg, artist_ra),
            ('Sunsets & Sunrises', 'Transits of Full Moon', 'Civilized Twilight', 'Astronomical Twilight', 'Local Sidereal Time'),
            loc = 'center left',
            )

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
    figname = 'night_chart_wendelstein.png'
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
