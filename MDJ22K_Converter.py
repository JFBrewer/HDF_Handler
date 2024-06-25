#!/usr/bin/env python

import math
import numpy as np
import datetime

def mj2k2gregorianday(mj2kday):
    mjd2k_conv = 51544
    mjd_conv = 2400000.5
    igreg = 2299161

    julday = mj2kday + mjd_conv + mjd2k_conv
    jday = math.floor(julday + 0.5)
    fracday = julday - jday + 0.5

    if jday >= igreg:
        jalpha = math.floor((jday - 1867216.25) / 36524.25)
        ja = jday + 1 + jalpha - math.floor(0.25 * jalpha)
    else:
        ja = jday

    jb = ja + 1524
    jc = math.floor((jb - 122.1) / 365.25)
    jd = math.floor(365.25 * jc)
    je = math.floor((jb - jd) / 30.6001)
    day = math.floor(jb - jd - math.floor(30.6001 * je))
    month = je - 1 if je <= 13 else je - 13
    year = jc - 4715 if month < 3 else jc - 4716
    
    time = (datetime.datetime(2000, 1, 1) + datetime.timedelta(days=fracday)).strftime('%H:%M:%S')
    date_string = f"{year}/{str(month).zfill(2)}/{str(day).zfill(2)} {time}"
    
    return date_string

def mj2k2gregorianday_vectorized(mj2kday):
    # This doesn't quite work for reasons unknown.
    mjd2k_conv = 51544
    mjd_conv = 2400000.5
    igreg = 2299161
    # Conversion to Julian day
    julday = mj2kday + mjd_conv + mjd2k_conv
    # Get fractional part of the Julian day
    jday = np.floor(julday + 0.5)
    fracday = julday - jday + 0.5
    jalpha = np.floor((jday - 1867216.25) / 36524.25)
    ja = jday + 1 + jalpha - np.floor(0.25 * jalpha)
    jb = np.where(jday >= igreg, ja, jday)
    jc = np.floor((jb - 122.1) / 365.25)
    jd = np.floor(365.25 * jc)
    je = np.floor((jb - jd) / 30.6001)
    day = np.floor(jb - jd - np.floor(30.6001 * je))
    # months are in the range [1, 12]
    ii = np.where(je <= 13, je - 1, je - 13)
    ii = np.where(ii < 3, ii + 4715, ii + 4716)
    year = np.where(ii < 3, jc - 4715, jc - 4716)
    leap = np.mod(year, 4)
    # time of the day
    time = [(datetime.datetime(2000, 1, 1) + datetime.timedelta(days=fday)).strftime('%H:%M:%S') for fday in fracday]
    #time = [datetime.datetime.strptime(str(fday), '%H:%M:%S') for fday in fracday]
    str_month = np.array(['{:02d}'.format(month) for month in ii])
    str_day = np.array(['{:02d}'.format(day) for day in day])
    # date in string format 'yyyy/mm/dd HH:MM:SS'
    date_string = np.core.defchararray.add(np.array(['{:04d}'.format(year) for year in year]),
                                            '/', str_month, '/', str_day, ' ', time)
    return date_string
