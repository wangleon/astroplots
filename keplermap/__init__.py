#!/usr/bin/env python
import os,math
import numpy as np
import scipy.interpolate as intp
import xml.dom.minidom as md

import matplotlib as mpl
from matplotlib.axes import Axes
from matplotlib.patches import Polygon, Circle
from matplotlib.collections import PatchCollection

def get_kepler_center(season=1):
    ra0_sum = 0.0
    de0_sum = 0.0
    pth = os.path.dirname(__file__)
    filename = '%s/data/KeplerFOVseason%d.kml'%(pth,season)
    dom = md.parse(filename)
    root = dom.documentElement
    pms = root.getElementsByTagName('Placemark')
    bound_lst = []
    for pm in pms:
        name = pm.getElementsByTagName('name')[0]
        modname = name.childNodes[0].data
        if 'Mod.out 13.' in modname:
            coord = pm.getElementsByTagName('coordinates')[0]
            text = coord.childNodes[0].data
            g = text.split()
            ra_sum = 0.
            de_sum = 0.
            for i in range(4):
                k = g[i].split(',')
                ra  = float(k[0])+180.0
                dec = float(k[1])
                ra_sum += ra
                de_sum += dec
            ra_center = ra_sum/4.0
            de_center = de_sum/4.0
            ra0_sum += ra_center
            de0_sum += de_center
    return (ra0_sum/4.0,de0_sum/4.0)

def get_kepler_field(season=1):
    pth = os.path.dirname(__file__)
    filename = '%s/data/KeplerFOVseason%d.kml'%(pth,season)
    dom = md.parse(filename)
    root = dom.documentElement
    pms = root.getElementsByTagName('Placemark')
    ccd_lst = []
    for pm in pms:
        name = pm.getElementsByTagName('name')[0]
        coord = pm.getElementsByTagName('coordinates')[0]
        text = coord.childNodes[0].data
        g = text.split()
        coord_lst = []
        for c in g[0:-1]:
            k = c.split(',')
            ra  = float(k[0])+180.0
            dec = float(k[1])
            coord_lst.append([ra,dec])
        coord_lst = np.array(coord_lst)
        ccd_lst.append(coord_lst)
    return ccd_lst

class KeplerMap(Axes):
    def __init__(self, fig, rect=None, **kwargs):

        if rect == None:
            rect = [0.1,0.1,0.8,0.8]
        self.season = kwargs.pop('season',1)
        self.ra_direction = kwargs.pop('ra_direction',1)
        # center coordinate in (RA, Dec)
        self.center = get_kepler_center(season=self.season)
        # self.center must be set before Axes.__init__, because Axes.__init__
        # calls axes.cla(), which in turns calls set_ylim(), which is
        # overwritten in this class.

        self.ref_dec1, self.ref_dec2 = 20., 50.

        Axes.__init__(self, fig, rect, **kwargs)

        # prepare for the draw the ticks on the top & right axis
        self.twiny = self.twiny()
        self.twinx = self.twinx()

        #self.set_aspect(1.0)

        self.set_ra_ticks(5)
        self.set_dec_ticks(5)

        self.set_lim(self.center[0]-11, self.center[0]+11,
                     self.center[1]-9,  self.center[1]+9)
        #self.grid()

    def __lambert(self,ra,dec):
        '''
        Lambert projection
        see http://en.wikipedia.org/wiki/Lambert_conformal_conic_projection
        '''
        p     = math.pi
        ln    = np.log
        power = np.power

        sin = np.sin
        cos = np.cos
        tan = np.tan
        cot = lambda(x): 1./tan(x)
        sec = lambda(x): 1./cos(x)

        ra   = np.deg2rad(ra)
        dec  = np.deg2rad(dec)
        ra0  = np.deg2rad(self.center[0])
        dec0 = np.deg2rad(self.center[1])

        dec1 = np.deg2rad(self.ref_dec1)
        dec2 = np.deg2rad(self.ref_dec2)

        n = ln(cos(dec1)*sec(dec2))/ln(tan(p/4+dec2/2)*cot(p/4+dec1/2))
        F = cos(dec1)*power(tan(p/4+dec1/2),n)/n
        rho0 = F*power(cos(p/4+dec0/2)/sin(p/4+dec0/2),n)
        rho  = F*power(cos(p/4+dec/2)/sin(p/4+dec/2),n)
        x = rho*sin(n*(ra-ra0))
        y = rho0 - rho*cos(n*(ra-ra0))
        return x,y


    def __reverse_lambert(self,x,y):

        p     = math.pi
        ln    = np.log
        power = np.power

        sin = np.sin
        cos = np.cos
        tan = np.tan
        cot = lambda(x): 1./tan(x)
        sec = lambda(x): 1./cos(x)

        ra0  = np.deg2rad(self.center[0])
        dec0 = np.deg2rad(self.center[1])

        dec1 = np.deg2rad(self.ref_dec1)
        dec2 = np.deg2rad(self.ref_dec2)

        n = ln(cos(dec1)*sec(dec2))/ln(tan(p/4+dec2/2)*cot(p/4+dec1/2))
        F = cos(dec1)*power(tan(p/4+dec1/2),n)/n
        rho0 = F*power(cos(p/4+dec0/2)/sin(p/4+dec0/2),n)
        
        rho = np.sign(n)*np.sqrt(x**2+(rho0-y)**2)
        theta = np.arctan(x/(rho0-y))
        dec = 2.*np.arctan(np.power((F/rho),(1./n))) - p/2
        ra = ra0 + theta/n
        return np.rad2deg(ra), np.rad2deg(dec)

    def plot(self, *args, **kwargs):
        '''
        a wrap of matplotlib.axes.plot() method
        '''
        x1,x2 = self.get_xlim()
        y1,y2 = self.get_ylim()
        x, y = self.__lambert(args[0], args[1])
        super(KeplerMap, self).plot(x, y, *args[2:], **kwargs)
        self.set_xlim(x1,x2)
        self.set_ylim(y1,y2)

    def scatter(self, *args, **kwargs):
        x1,x2 = self.get_xlim()
        y1,y2 = self.get_ylim()
        x, y = self.__lambert(args[0], args[1])
        super(KeplerMap, self).scatter(x, y, *args[2:], **kwargs)
        self.set_xlim(x1,x2)
        self.set_ylim(y1,y2)

    def text(self, *args, **kwargs):
        '''
        a wrap of matplotlib.axes.text() method
        '''
        x1,x2 = self.get_xlim()
        y1,y2 = self.get_ylim()
        x, y = self.__lambert(args[0], args[1])
        super(KeplerMap, self).text(x, y, *args[2:], **kwargs)
        self.set_xlim(x1,x2)
        self.set_ylim(y1,y2)


    def set_lim(self, ra1, ra2, dec1, dec2):
        self.ra1,  self.ra2  = ra1,  ra2
        self.dec1, self.dec2 = dec1, dec2
        _,y1 = self.__lambert(self.center[0],dec1)
        _,y2 = self.__lambert(self.center[0],dec2)
        x1 = self.__find_x(ra1, y1, dec1)
        x2 = self.__find_x(ra2, y1, dec1)
        self.set_xlim(x1,x2)
        self.set_ylim(y1,y2)

        self.__refresh_ticks()

    def set_ra_ticks(self, major, minor=None):
        self.ra_major = abs(major)
        self.ra_minor = minor

    def set_dec_ticks(self, major, minor=None):
        self.dec_major = abs(major)
        self.dec_minor = minor

    def __refresh_ticks(self):
        x1,x2 = self.get_xlim()
        y1,y2 = self.get_ylim()

        # refresh RA major ticks
        self.xmajor_ticks = []
        if x1 < x2:
            self.ra_majors = np.arange(0,360+1e-7,self.ra_major)
        else:
            self.ra_majors = np.arange(360,0-1e-7,-self.ra_major)
        xticklabels = []
        for ra in self.ra_majors:
            x = self.__find_x(ra, y1, self.dec1)
            if x == None:
                continue
            elif (x > x1) == (x < x2):
                self.xmajor_ticks.append(x)
                xticklabels.append(str(ra)+u'\xb0')
        self.set_xticks(self.xmajor_ticks)
        self.get_xaxis().set_tick_params(direction='out')
        self.set_xticklabels(xticklabels)

        # refresh RA minor ticks
        self.xminor_ticks = []
        self.ra_minors = []
        if self.ra_minor != None:
            if x1 < x2:
                lst = np.arange(0, 360+1e-7, self.ra_minor)
            else:
                lst = np.arange(360, 0-1e-7, -self.ra_minor)
            for ra in lst:
                if ra not in self.ra_majors:
                    self.ra_minors.append(ra)
                    x = self.__find_x(ra, y1, self.dec1)
                    if x == None:
                        continue
                    elif (x > x1) == (x < x2):
                        self.xminor_ticks.append(x)

        # set the ticks on the top x-axis
        self.twiny.set_xlim(x1,x2)
        self.twiny.set_ylim(y1,y2)
        self.xmajor_ticks2 = []
        xticklabels2 = []
        for ra in self.ra_majors:
            x = self.__find_x(ra, y2, self.dec2)
            if x == None:
                continue
            elif (x > x1) == (x < x2):
                self.xmajor_ticks2.append(x)
                xticklabels2.append(str(ra)+u'\xb0')
        self.twiny.set_xticks(self.xmajor_ticks2)
        self.twiny.get_xaxis().set_tick_params(direction='out')
        self.twiny.set_xticklabels(xticklabels2)


        # refresh DEC ticks
        self.ymajor_ticks = []
        if y1 < y2:
            self.dec_majors = np.arange(-90,90+1e-7,self.dec_major)
        else:
            self.dec_majors = np.arange(90,-90-1e-7,-self.dec_major)
        yticklabels = []
        for dec in self.dec_majors:
            if abs(abs(dec)-90)<1e-4:
                continue
            y = self.__find_y(dec, x1, self.ra1)
            if y == None:
                continue
            elif (y > y1) == (y < y2):
                self.ymajor_ticks.append(y)
                yticklabels.append(str(dec)+u'\xb0')
        self.set_yticks(self.ymajor_ticks)
        self.get_yaxis().set_tick_params(direction='out')
        self.set_yticklabels(yticklabels)

        # refresh DEC minor ticks
        self.yminor_ticks = []
        self.dec_minors = []
        if self.dec_minor !=None:
            if y1 < y2:
                lst = np.arange(-90,90+1e-7,self.dec_minor)
            else:
                lst = np.arange(90,-90-1e-7,-self.dec_minor)
            for dec in lst:
                if abs(abs(dec)-90)<1e-4:
                    continue
                elif dec not in self.dec_majors:
                    self.dec_minors.append(dec)
                    y = self.__find_y(dec, x1, self.ra1)
                    if y == None:
                        continue
                    elif (y > y1) == (y < y2):
                        self.yminor_ticks.append(y)

        # refresh right DEC major ticks
        self.twinx.set_xlim(x1,x2)
        self.twinx.set_ylim(y1,y2)
        self.ymajor_ticks2 = []
        yticklabels2 = []
        for dec in self.dec_majors:
            if abs(abs(dec)-90)<1e-4:
                continue
            y = self.__find_y(dec, x2, self.ra2)
            if y == None:
                continue
            elif (y > y1) == (y < y2):
                self.ymajor_ticks2.append(y)
                yticklabels2.append(str(dec)+u'\xb0')
        self.twinx.set_yticks(self.ymajor_ticks2)
        self.twinx.get_yaxis().set_tick_params(direction='out')
        self.twinx.set_yticklabels(yticklabels2)

    def plot_kepler_field(self, *args, **kwargs):
        '''plot kepler fields'''
        ccd_lst = get_kepler_field()
        for ccd in ccd_lst:
            self.plot_polygon(ccd, *args, **kwargs)

    def plot_polygon(self, eqcoords, **kwargs):
        xy = []
        for coord in eqcoords:
            x,y = self.__lambert(coord[0], coord[1])
            xy.append([x,y])
        polygon = Polygon(np.array(xy), True, **kwargs)
        self.add_patch(polygon)

    def grid(self, b=None, which='major',axis='both',**kwargs):
        color     = kwargs.pop('color',mpl.rcParams['grid.color'])
        alpha     = kwargs.pop('alpha',mpl.rcParams['grid.alpha'])
        linestyle = kwargs.pop('linestyle',
                    kwargs.pop('ls',mpl.rcParams['grid.linestyle']))
        linewidth = kwargs.pop('linewidth',
                    kwargs.pop('lw',mpl.rcParams['grid.linewidth']))

        mpl.rcParams['axes.grid'] = False

        x1,x2 = self.get_xlim()
        y1,y2 = self.get_ylim()

        if b in [None,True,'on'] and axis in ['both','y']:
            # plot horizantal lines
            for dec in self.dec_majors:
                if abs(abs(dec)-90)<1e-6:
                    continue
                ra_lst = np.arange(0,360,2.)
                dec_lst = np.zeros_like(ra_lst)+dec
                self.plot(ra_lst,dec_lst,color=color,
                          ls=linestyle,lw=linewidth,alpha=alpha)

        # plot vertical lines
        if b in [None,True,'on'] and axis in ['both','x']:
            for ra in self.ra_majors:
                dec_lst = np.arange(-88,89,2.)
                ra_lst = np.zeros_like(dec_lst)+ra
                self.plot(ra_lst,dec_lst,color=color,
                          ls=linestyle,lw=linewidth,alpha=alpha)

        self.set_xlim(x1,x2)
        self.set_ylim(y1,y2)

    def __find_x(self, ra, y, dec0):
        dec1 = dec0
        niter = 0
        max_iter = 10
        while(True):
            niter += 1
            x1, y1   = self.__lambert(ra,dec1)
            ra1,dec1 = self.__reverse_lambert(x1,y)
            if abs(ra-ra1)<1e-4:
                break
            if niter >= max_iter:
                return None
        return x1

    def __find_y(self, dec,x,ra0):
        ra1 = ra0
        niter = 0
        max_iter = 10
        while(True):
            niter += 1
            x1, y1   = self.__lambert(ra1,dec)
            ra1,dec1 = self.__reverse_lambert(x,y1)
            if abs(dec-dec1)<1e-4:
                break
            if niter >= max_iter:
                return None
        return y1

    def plot_field_stars(self, V_limit=None, *args, **kwargs):
        '''
        plot background stars in Kepler field
        '''
        path = os.path.dirname(__file__)
        filename = '%s/data/field_stars.dat'%path

        ra_lst   = []
        dec_lst  = []
        vmag_lst = []

        infile = open(filename)
        for row in infile:
            row = row.strip()
            if len(row)==0 or row[0] in '#%':
                continue
            g = row.split('|')
            HD        = int(g[0])
            HR        = int(g[1])
            Flamsteed = g[2].strip().decode('utf-8')
            Bayer     = g[3].strip().decode('utf-8')
            ra        = float(g[4])
            dec       = float(g[5])
            vmag      = float(g[6])
            bv        = float(g[7])
            ra_lst.append(ra)
            dec_lst.append(dec)
            vmag_lst.append(vmag)
        infile.close()
        vmag_lst = np.array(vmag_lst)

        s_lst = kwargs.pop('s',(8-vmag_lst)**2)
        c_lst = kwargs.pop('c',['k' for v in ra_lst])

        self.scatter(ra_lst, dec_lst, *args, s=s_lst, c=c_lst, **kwargs)

    def plot_clusters(self):

        path = os.path.dirname(__file__)
        filename = '%s/data/clusters.dat'%path

        ngc_lst  = []
        ra_lst   = []
        dec_lst  = []
        radius_lst = []

        infile = open(filename)
        for row in infile:
            row = row.strip()
            if len(row)==0 or row[0] in '#%':
                continue
            g = row.split('|')
            ngc    = int(g[0])
            ra     = float(g[1])
            dec    = float(g[2])
            radius = float(g[3])/60.
            self.plot_circle(ra,dec, radius, color='k', ls=':',lw=0.5)
        infile.close()
        
    def plot_circle(self,ra0,dec0,radius,*args,**kwargs):
        '''
        plot a circle with given radius in Kepler field
        '''
        #pseudo (ra,dec) of lamost field
        ra_lst = np.arange(0,360.001,5)
        dec_lst = np.zeros_like(ra_lst)+90.0-radius

        # real (ra,dec) of lamost field
        ra2_lst  = []
        dec2_lst = []
        for i in range(len(ra_lst)):
            ra  = ra_lst[i]
            dec = dec_lst[i]
            ra2,dec2 = self.__rotate_sphere(ra,dec,dec0,ra0)
            ra2_lst.append(ra2)
            dec2_lst.append(dec2)

        #x_lst, y_lst = self.__lambert(ra2_lst,dec2_lst)
        self.plot(ra2_lst, dec2_lst, *args, **kwargs)

    def __rotate_sphere(self,ra,dec,DELTA_NGP,L_NCP):
        '''rotate shpere'''
        ALPHA_NGP = 0.0
        AI     = (90. - DELTA_NGP)/180.*math.pi
        ALPHA0 = (ALPHA_NGP + 90.)/180.*math.pi
        LON0   = (L_NCP - 90.)/180.*math.pi

        MIN    = 1e-10

        ra  = ra/180.*math.pi
        dec = dec/180.*math.pi

        ALat = math.sin(dec)*math.cos(AI) - \
               math.cos(dec)*math.sin(AI)*math.sin(ra-ALPHA0)
        BLat = math.sqrt( 1 - ALat*ALat )
        L1   = (math.sin(dec)*math.sin(AI) + \
                math.cos(dec)*math.cos(AI)*math.sin(ra-ALPHA0)) / BLat
        L2   = (math.cos(dec)*math.cos(ra-ALPHA0)) / BLat
        if L1>1.0: L1=1.0
        if abs(L1) < MIN:
            if L2 > 0.0: L0 = 0.0
            else:        L0 = math.pi
        else:
            if L1 >= 0.0:
                if L2 >= 0.0:
                    L0 = math.asin(L1)
                else:
                    L0 = math.pi - math.asin(L1)
            else:
                if L2 >= 0.0:
                    L0 = 2.0*math.pi + math.asin(L1)
                else:
                    L0 = math.pi - math.asin(L1)

        L = L0 + LON0
        #if L >= 2.0*math.pi:
        #    L = L - 2.0*math.pi
        Lat = math.asin(ALat)
        l = L/math.pi*180.
        b = Lat/math.pi*180.

        return (l,b)

