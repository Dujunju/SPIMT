"""
--PhotSim.py--
构建观测场景，生成仿真图象
最后修改：2021-01-06
"""
from abc import ABC

import numpy as np
import numpy.random as rd
import pandas as pd
from pandas import DataFrame
from skyfield.api import Topos, EarthSatellite, load
from astropy.nddata import CCDData
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.io.fits as fits
from astropy.time import Time, TimeDelta
from astropy.table import Table
import math
from astropy.stats import SigmaClip
ts = load.timescale()
from photutils import Background2D, SExtractorBackground
from scipy import interpolate
from astropy.stats import gaussian_fwhm_to_sigma
from astroquery.xmatch import XMatch


class Station(Topos):
    def __init__(self, *args, **kwargs):
        super(Station, self).__init__(*args, **kwargs)


class Target(EarthSatellite):
    def __init__(self, *args, magnitude=15, **kwargs):
        super(Target, self).__init__(*args, **kwargs)
        self.magnitude = magnitude


class Stars(DataFrame, ABC):
    def __init__(self, *args, **kwargs):
        super(Stars, self).__init__(*args, **kwargs)


class Telescope:
    def __init__(self, ra_sigma=0.02, dec_sigma=0.02, fov=0.2, zero_point=2.4, k=0.17):
        self.ra_sigma = ra_sigma
        self.dec_sigma = dec_sigma
        self.fov = fov*u.deg
        self.zero_point = zero_point
        self.k = k


class CCD:
    def __init__(self, shape=(2048, 2048), gain=1.63, plate_const=[], cd=[], pixel_size=0.3533, pixel_scale=1):
        self.shape = shape
        self.gain = gain
        self.plate_const = np.array(plate_const)
        self.cd = cd
        self.pixel_size=pixel_size
        self.pixel_scale = pixel_scale


class SynFits:

    def __init__(self, station=None, target=None, telescope=None, CCD=None, stars=None,
                 time_init=None, exposure=None, tracking_mode='target', A=None, D=None, dA=0, dD=0):
        self.station = station
        self.target = target
        self.telescope = telescope
        self.ccd = CCD
        self.stars = stars
        self.time_init = Time(time_init, scale='utc') if time_init else None
        self.exposure = TimeDelta(exposure, format='sec') if exposure else None
        self.tracking_mode = tracking_mode
        self.station2target = target - station if (target and station) else None
        self.time_mid = time_init + 0.5*exposure if (time_init and exposure) else None
        self.A = A
        self.D = D
        self.radec_init = None
        self.radec_mid = None
        self.radec_end = None
        self.dA = dA
        self.dD = dD

    def add_station(self, *args, **kwargs):
        print('add station ...')
        self.station = Station(*args, **kwargs)

    def add_target(self, *args, **kwargs):
        print('add target ...')
        self.target = Target(*args, **kwargs)

    def add_telescope(self, *args,  **kwargs):
        print('add telescope ...')
        self.telescope = Telescope(*args, **kwargs)

    def add_ccd(self, *args, **kwargs):
        print('add ccd ...')
        self.ccd = CCD(*args, **kwargs)

    def set_setting(self, seeing=2, time_init=None, exposure=None, tracking_mode=None, A=None, D=None, dA=0, dD=0):
        print('set setting ...')
        self.seeing = seeing
        self.time_init = Time(time_init, scale='utc')
        self.exposure = TimeDelta(exposure, format='sec')
        self.tracking_mode = tracking_mode
        self.time_mid = self.time_init + 0.5*self.exposure
        self.A = A
        self.D = D
        self.dA = dA
        self.dD = dD

    def add_stars(self, mag_range=None, **kwargs):
        print('add_stars ...')

        self.station2target = self.target - self.station

        # 计算曝光中间时刻(time_mid)，望远镜的指向（A_mid, D_mid）
        if self.tracking_mode == 'target':
            A_mid, D_mid, distance = self.station2target.at(ts.from_astropy(self.time_mid)).radec()
            A_mid = A_mid._degrees+self.dA
            D_mid = D_mid._degrees+self.dD
        elif self.tracking_mode == 'sidereal':
            A_mid, D_mid = self.A, self.D
        elif self.tracking_mode == 'parking':
            A_mid, D_mid = self.A + 15 * 0.5 * self.exposure.value/3600, self.D
        else:
            print('The tracking mode can not be recognized!')

        # 利用Simbad，获取视场中的背景恒星
        radec_mid = SkyCoord(A_mid, D_mid, unit=(u.deg, u.deg))
        self.radec_init = self.station2target.at(ts.from_astropy(self.time_init)).radec()
        self.radec_mid = radec_mid

        # 利用Vizier，获取视场中的背景恒星
        custom_visizer = Vizier()
        custom_visizer.ROW_LIMIT = -1
        stars = custom_visizer.query_region(radec_mid, radius=self.telescope.fov/2, catalog='UCAC4', **kwargs)
        stars = XMatch.query(cat1=stars[0], cat2='vizier:I/350/gaiaedr3', max_distance=1*u.arcsec, colRA1='RAJ2000', colDec1='DEJ2000')
        stars = stars.to_pandas()
        stars = stars[(stars['Vmag'] > mag_range[0]) & (stars['Vmag'] < mag_range[1])]
        stars['RAJ2000'] = stars['ra'] + stars['pmra']*(self.time_mid.jyear-2016)/1000/60/60/np.cos(np.deg2rad(stars['dec']))
        stars['DEJ2000'] = stars['dec'] + stars['pmdec']*(self.time_mid.jyear-2016)/ 1000 / 60 / 60
        stars = stars[['UCAC4', 'RAJ2000', 'DEJ2000', 'Vmag']]
        stars = stars.rename(columns={'UCAC4': 'ID', 'RAJ2000': 'RA', 'DEJ2000': 'DEC', 'Vmag': 'FLUX_V'})
        stars = stars.dropna()

        # 把空间目标添加进恒星列表
        stars = stars.append([{'FLUX_V': self.target.magnitude, 'ID': 'TARGET'}])
        self.stars = Stars(stars)

    def add_photon(self):
        print('add photons ...')
        self.station2target = self.target - self.station

        # 计算大气质量
        alt_mid, az_mid, distance_mid = self.station2target.at(ts.from_astropy(self.time_mid)).altaz()
        airmass = 1 / np.cos(0.5 * np.pi - alt_mid.radians)

        # 计算曝光时间内，到达焦平面的光子数
        self.stars['FLUX'] = 10 ** (-0.4 * (self.stars['FLUX_V']
                                                  - 25 + self.telescope.zero_point
                                                  + self.telescope.k * airmass))
        self.stars['FLUX_exp'] = np.random.poisson(self.stars['FLUX']*self.exposure.value)

        # 记录光子的信息
        photons = pd.DataFrame()
        photons['ID'] = self.stars['ID']
        photons['FLUX_exp'] = self.stars['FLUX_exp']
        photons['T_IJ'] = photons['FLUX_exp'].map(lambda x: np.random.uniform(0, self.exposure.value, x))
        self.photons = pd.DataFrame({'ID': photons['ID'].repeat(photons['T_IJ'].str.len()),
                                     'T_IJ': np.concatenate(photons['T_IJ'].values)})

    def tracking(self, inter=True):
        print('tracking ... ')

        # 将光子按到达时间排序
        photons = self.photons.sort_values(by='T_IJ')

        # 生成二维布朗运动
        photons['DT'] = photons['T_IJ'].diff(1).fillna(0)
        photons['DT_CUM'] = photons['DT'].cumsum()
        photons['DRA'] = photons['DT'].map(lambda x: self.telescope.ra_sigma * np.random.normal(0, x**0.5, 1).item())
        photons['DDEC'] = photons['DT'].map(lambda x: self.telescope.dec_sigma * np.random.normal(0, x**0.5, 1).item())
        photons['DRA'].iloc[0] = self.dA
        photons['DDEC'].iloc[0] = self.dD
        photons['DRA_CUM'] = photons['DRA'].cumsum()
        photons['DDEC_CUM'] = photons['DDEC'].cumsum()

        if self.tracking_mode == 'target':
            if inter:
                t = ts.from_astropy(self.time_init+TimeDelta(np.arange(0, self.exposure.value+0.1, 0.1), format='sec'))
                RA_, Dec_, dis = self.station2target.at(t).radec()

                RA_tck = interpolate.splrep(np.arange(0, self.exposure.value+0.1, 0.1), RA_._degrees, s=0)
                RA_Tar = interpolate.splev(photons['T_IJ'], RA_tck, der=0)

                Dec_tck = interpolate.splrep(np.arange(0, self.exposure.value + 0.1, 0.1), Dec_._degrees, s=0)
                Dec_Tar = interpolate.splev(photons['T_IJ'], Dec_tck, der=0)
            else:
                t = ts.from_astropy(self.time_init+TimeDelta(photons['T_IJ'].values, format='sec'))
                RA_Tar = []
                Dec_Tar = []
                num = len(photons) // 10000
                for i in range(num + 1):
                    print(i, '//', num)
                    RA_, Dec_, dis = self.station2target.at(t[i * 10000:(i + 1) * 10000]).radec()
                    RA_Tar.extend(RA_._degrees)
                    Dec_Tar.extend(Dec_._degrees)

            photons['RA_TAR'] = RA_Tar
            photons['DEC_TAR'] = Dec_Tar

            photons['A'] = photons['RA_TAR'] + photons['DRA_CUM']
            photons['D'] = photons['DEC_TAR'] + photons['DDEC_CUM']

        elif self.tracking_mode == 'sidereal':

            t = ts.from_astropy(self.time_init+TimeDelta(photons['T_IJ'].loc[photons['ID']=='TARGET'].values, format='sec'))
            RA_, Dec_, dis = self.station2target.at(t).radec()
            photons['RA_TAR'] = 0
            photons['DEC_TAR'] = 0
            photons['RA_TAR'].loc[photons['ID'] == 'TARGET'] = RA_._degrees
            photons['DEC_TAR'].loc[photons['ID'] == 'TARGET'] = Dec_._degrees
            photons['A'] = self.A + photons['DRA_CUM']
            photons['D'] = self.D + photons['DDEC_CUM']

        elif self.tracking_mode == 'parking':
            t = ts.from_astropy(self.time_init+TimeDelta(photons['T_IJ'].loc[photons['ID']=='TARGET'].values, format='sec'))
            RA_, Dec_, dis = self.station2target.at(t).radec()
            photons['RA_TAR'] = 0
            photons['DEC_TAR'] = 0
            photons['RA_TAR'].loc[photons['ID'] == 'TARGET'] = RA_._degrees
            photons['DEC_TAR'].loc[photons['ID'] == 'TARGET'] = Dec_._degrees
            photons['A'] = self.A + 15*photons['T_IJ']/3600 + photons['DRA_CUM']
            photons['D'] = self.D + photons['DDEC_CUM']
        else:
            print('The tracking mode can not be recognized!')

        photon_map = pd.merge(self.stars[['RA', 'DEC', 'ID']], photons, on='ID')

        photon_map['RA'] = np.where(photon_map['ID'] == 'TARGET', photon_map['RA_TAR'], photon_map['RA'])
        photon_map['DEC'] = np.where(photon_map['ID'] == 'TARGET', photon_map['DEC_TAR'], photon_map['DEC'])

        T1 = np.cos(np.deg2rad(photon_map['DEC'])) * np.sin(np.deg2rad(photon_map['RA'] - photon_map['A']))
        T2 = np.sin(np.deg2rad(photon_map['DEC'])) * np.cos(np.deg2rad(photon_map['D'])) - \
             np.cos(np.deg2rad(photon_map['DEC'])) * np.sin(np.deg2rad(photon_map['D'])) * \
             np.cos(np.deg2rad(photon_map['RA'] - photon_map['A']))
        T3 = np.sin(np.deg2rad(photon_map['DEC'])) * np.sin(np.deg2rad(photon_map['D'])) + \
             np.cos(np.deg2rad(photon_map['DEC'])) * np.cos(np.deg2rad(photon_map['D'])) * \
             np.cos(np.deg2rad(photon_map['RA'] - photon_map['A']))
        xi = T1 / T3
        eta = T2 / T3

        photon_map['XI'] = xi
        photon_map['ETA'] = eta

        # 理想坐标到测量坐标
        xy = np.dot(np.linalg.inv(self.ccd.cd), np.rad2deg(photon_map[['XI', 'ETA']].values.T))
        photon_map['X'] = xy[0, :]
        photon_map['Y'] = xy[1, :]

        photon_map['X_SE'] = photon_map['X'] + rd.normal(0, self.seeing, len(photon_map['X']))
        photon_map['Y_SE'] = photon_map['Y'] + rd.normal(0, self.seeing, len(photon_map['Y']))

        self.photon_map = photon_map

    def rendering(self, bias=0, background=0, dx=0, dy=0):
        print('rendering ...')

        img_photon = np.zeros(self.ccd.shape)
        x_size = self.ccd.shape[0]
        y_size = self.ccd.shape[1]

        for x, y in zip(self.photon_map['X_SE'], self.photon_map['Y_SE']):
            x = x + dx
            y = y + dy
            if (-0.5 < x < x_size) and (-0.5 < y < y_size):
                img_photon[int(y), int(x)] = img_photon[int(y), int(x)] + 1

        self.img_photon = img_photon
        self.image_adu = (img_photon/self.ccd.gain + background)

    def writeto(self, img_photon='photon.fits', img_adu='adu.fits', header=None, **kwargs):
        print('write to ...')

        hdu = fits.PrimaryHDU(self.img_photon.astype(np.float32), header=header)
        hdu.writeto(img_photon, **kwargs)

        hdu = fits.PrimaryHDU(self.image_adu.astype(np.float32), header=header)
        hdu.writeto(img_adu, **kwargs)


if __name__ == '__main__':
   pass
