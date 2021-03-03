"""
spimt.py
A python code for simulating photometric images of moving target with photon-mapping.
Author: Junju Du
E-mail: dujunju@mail.sdu.edu.cn
Last update: 2021-03-03
"""

from abc import ABC
import numpy as np
import numpy.random as rd
import pandas as pd
from pandas import DataFrame
from skyfield.api import Topos, EarthSatellite, load
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.io.fits as fits
from astropy.time import Time, TimeDelta
from scipy import interpolate
from astroquery.xmatch import XMatch
ts = load.timescale()


class Station(Topos):
    """
    A ground station class.
    """
    def __init__(self, *args, **kwargs):
        """
        :param args: the arguments of skyfield.api.Topos.
        :param kwargs: the keywords of skyfield.api.Topos.
        """
        super(Station, self).__init__(*args, **kwargs)


class Target(EarthSatellite):
    """
    A moving target class.
    """
    def __init__(self, *args, magnitude=15, **kwargs):
        """
        :param args: the arguments of skyfield.api.EarthSatellite.
        :param magnitude: the standard magnitude of the moving target.
        :param kwargs: the keywords of skyfield.api.EarthSatellite.
        """
        super(Target, self).__init__(*args, **kwargs)
        self.magnitude = magnitude


class Stars(DataFrame, ABC):
    """
    Field stars class.
    """
    def __init__(self, *args, **kwargs):
        """
        :param args: the arguments of the pandas.DataFrame.
        :param kwargs: the keywords of the pandas.DataFrame.
        """
        super(Stars, self).__init__(*args, **kwargs)


class Telescope:
    """
    A telescope class.
    """
    def __init__(self, ra_sigma=0.02, dec_sigma=0.02, fov=0.2, zero_point=2.4, k=0.17):
        """
        :param ra_sigma: the scale parameter of Brownian motion in right ascension, in pixel/second.
        :param dec_sigma: the scale parameter of Brownian motion in declination, in pixel/second.
        :param fov: the field of view of telescope, in degree.
        :param zero_point: the zero-point of telescope.
        :param k: the first order extinction coefficient.
        """
        self.ra_sigma = ra_sigma
        self.dec_sigma = dec_sigma
        self.fov = fov * u.deg
        self.zero_point = zero_point
        self.k = k


class CCD:
    """
    A CCD class.
    """
    def __init__(self, shape=(2048, 2048), gain=1.63, plate_const=None, cd=None, pixel_size=0.3533, pixel_scale=1):
        """
        :param shape: the dimensions of the CCD.
        :param gain: the gain of the CCD.
        :param plate_const: the plate constants of the CCD.
        :param cd: the WCS (CD1_1, CD1_2, CD2_1, CD2_2) of the CCD.
        :param pixel_size: the pixel size of the CCD.
        :param pixel_scale: the pixel scale of the CCD.
        """
        self.shape = shape
        self.gain = gain
        self.plate_const = np.array(plate_const)
        self.cd = cd
        self.pixel_size = pixel_size
        self.pixel_scale = pixel_scale


class SynFits:
    """
    An observation scene.
    """

    def __init__(self, station=None, target=None, telescope=None, ccd=None, stars=None,
                 time_init=None, exposure=None, tracking_mode='target', A0=None, D0=None, eA=0, eD=0):
        """
        :param station: the ground station.
        :param target: the moving target.
        :param telescope: the telescope.
        :param CCD: the CCD.
        :param stars: the field stars.
        :param time_init: the initial time of exposure.
        :param exposure: the exposure time.
        :param tracking_mode: the tracking mode of telescope.
        :param A0: the right ascension of telescope pointing.
        :param D0: the declination of telescope pointing.
        :param eA: the initial error in right ascension of telescope pointing.
        :param eD: the initial error in declination of telescope pointing.
        """
        self.station = station
        self.target = target
        self.telescope = telescope
        self.ccd = ccd
        self.stars = stars
        self.time_init = Time(time_init, scale='utc') if time_init else None
        self.exposure = TimeDelta(exposure, format='sec') if exposure else None
        self.tracking_mode = tracking_mode
        self.station2target = target - station if (target and station) else None
        self.time_mid = time_init + 0.5 * exposure if (time_init and exposure) else None
        self.A0 = A0
        self.D0 = D0
        self.radec_init = None
        self.radec_mid = None
        self.radec_end = None
        self.eA = eA
        self.eD = eD

    def add_station(self, *args, **kwargs):
        """
        Add a ground station to scene.
        :param args: the arguments of Station.
        :param kwargs: the keywords of Station.
        :return:
        """
        print('add station ...')
        self.station = Station(*args, **kwargs)

    def add_target(self, *args, **kwargs):
        """
        Add a moving target to scene.
        :param args: the arguments of Target.
        :param kwargs: the keyword of Target.
        :return:
        """
        print('add target ...')
        self.target = Target(*args, **kwargs)

    def add_telescope(self, *args, **kwargs):
        """
        Add a telescope to scene.
        :param args: the arguments of Telescope.
        :param kwargs: the keywords of Telescope.
        :return:
        """
        print('add telescope ...')
        self.telescope = Telescope(*args, **kwargs)

    def add_ccd(self, *args, **kwargs):
        """
        Add a CCD to scene.
        :param args: the arguments of CCD.
        :param kwargs: the keywords of CCD.
        :return:
        """
        print('add ccd ...')
        self.ccd = CCD(*args, **kwargs)

    def set_setting(self, seeing=2, time_init=None, exposure=None, tracking_mode=None, A0=None, D0=None, eA=0, eD=0):
        """
        Set the observation setting.
        :param seeing: the sigma of PSF. The PSF is a Gaussian function.
        :param time_init: the initial time of exposure.
        :param exposure: the exposure time.
        :param tracking_mode: the tracking mode of telescope.
        :param A0: the right ascension of telescope pointing.
        :param D0: the declination of telescope pointing.
        :param eA: the initial error in right ascension of telescope pointing.
        :param eD: the initial error in declination of telescope pointing.
        :return:
        """
        print('set setting ...')
        self.seeing = seeing
        self.time_init = Time(time_init, scale='utc')
        self.exposure = TimeDelta(exposure, format='sec')
        self.tracking_mode = tracking_mode
        self.time_mid = self.time_init + 0.5 * self.exposure
        self.A0 = A0
        self.D0 = D0
        self.eA = eA
        self.eD = eD

    def add_stars(self, mag_range=None, **kwargs):
        """
        Add field stars to scene.
        :param mag_range: the mag range of field stars.
        :param kwargs:
        :return:
        """
        print('add field stars ...')
        self.station2target = self.target - self.station

        # calculate the mid-exposure (time_mid)，and the corresponding telescope's pointing (A_mid, D_mid).
        if self.tracking_mode == 'target':
            A_mid, D_mid, distance = self.station2target.at(ts.from_astropy(self.time_mid)).radec()
            A_mid = A_mid._degrees + self.eA
            D_mid = D_mid._degrees + self.eD
        elif self.tracking_mode == 'sidereal':
            A_mid, D_mid = self.A0, self.D0
        elif self.tracking_mode == 'parking':
            A_mid, D_mid = self.A0 + 15 * 0.5 * self.exposure.value / 3600, self.D0
        else:
            print('The tracking mode can not be recognized!')

        # Obtain the field stars from Vizier.
        # Catalog is 'UCAC4'.
        # fileter is 'V'.
        radec_mid = SkyCoord(A_mid, D_mid, unit=(u.deg, u.deg))
        self.radec_init = self.station2target.at(ts.from_astropy(self.time_init)).radec()
        self.radec_mid = radec_mid

        custom_visizer = Vizier()
        custom_visizer.ROW_LIMIT = -1
        stars = custom_visizer.query_region(radec_mid, radius=self.telescope.fov / 2, catalog='UCAC4', **kwargs)
        stars = XMatch.query(cat1=stars[0], cat2='vizier:I/350/gaiaedr3', max_distance=1 * u.arcsec, colRA1='RAJ2000',
                             colDec1='DEJ2000')
        stars = stars.to_pandas()
        stars = stars[(stars['Vmag'] > mag_range[0]) & (stars['Vmag'] < mag_range[1])]
        stars['RAJ2000'] = stars['ra'] + stars['pmra'] * (self.time_mid.jyear - 2016) / 1000 / 60 / 60 / np.cos(
            np.deg2rad(stars['dec']))
        stars['DEJ2000'] = stars['dec'] + stars['pmdec'] * (self.time_mid.jyear - 2016) / 1000 / 60 / 60
        stars = stars[['UCAC4', 'RAJ2000', 'DEJ2000', 'Vmag']]
        stars = stars.rename(columns={'UCAC4': 'ID', 'RAJ2000': 'RA', 'DEJ2000': 'DEC', 'Vmag': 'FLUX_V'})
        stars = stars.dropna()

        # Add the target to source list.
        stars = stars.append([{'FLUX_V': self.target.magnitude, 'ID': 'TARGET'}])
        self.stars = Stars(stars)

    def add_photon(self):
        """
        Add photons from sources to scene.
        :return:
        """
        print('add photons ...')
        self.station2target = self.target - self.station

        # The air-mass at mid-exposure.
        alt_mid, az_mid, distance_mid = self.station2target.at(ts.from_astropy(self.time_mid)).altaz()
        airmass = 1 / np.cos(0.5 * np.pi - alt_mid.radians)

        # The number of photons survive to the image plane within exposure.
        self.stars['FLUX'] = 10 ** (-0.4 * (self.stars['FLUX_V']
                                            - 25 + self.telescope.zero_point
                                            + self.telescope.k * airmass))
        self.stars['FLUX_exp'] = np.random.poisson(self.stars['FLUX'] * self.exposure.value)

        # genarate the arrive time of photons and store photons' information.
        photons = pd.DataFrame()
        photons['ID'] = self.stars['ID']
        photons['FLUX_exp'] = self.stars['FLUX_exp']
        photons['T_IJ'] = photons['FLUX_exp'].map(lambda x: np.random.uniform(0, self.exposure.value, x))
        self.photons = pd.DataFrame({'ID': photons['ID'].repeat(photons['T_IJ'].str.len()),
                                     'T_IJ': np.concatenate(photons['T_IJ'].values)})

    def tracking(self, inter=True):
        """
        Photon tracing.
        :param inter: If interpolate the tracking path or not.
        :return:
        """
        print('tracking ... ')

        # Sort the photons by arrive time.
        photons = self.photons.sort_values(by='T_IJ')

        # Generate 2 2D Brownian motion.
        photons['DT'] = photons['T_IJ'].diff(1).fillna(0)
        photons['DT_CUM'] = photons['DT'].cumsum()
        photons['DRA'] = photons['DT'].map(lambda x: self.telescope.ra_sigma * np.random.normal(0, x ** 0.5, 1).item())
        photons['DDEC'] = photons['DT'].map(
            lambda x: self.telescope.dec_sigma * np.random.normal(0, x ** 0.5, 1).item())
        photons['DRA'].iloc[0] = self.eA
        photons['DDEC'].iloc[0] = self.eD
        photons['DRA_CUM'] = photons['DRA'].cumsum()
        photons['DDEC_CUM'] = photons['DDEC'].cumsum()

        # trace the photons according the tracking mode of telescope.
        if self.tracking_mode == 'target':
            if inter:
                t = ts.from_astropy(
                    self.time_init + TimeDelta(np.arange(0, self.exposure.value + 0.1, 0.1), format='sec'))
                RA_, Dec_, dis = self.station2target.at(t).radec()

                RA_tck = interpolate.splrep(np.arange(0, self.exposure.value + 0.1, 0.1), RA_._degrees, s=0)
                RA_Tar = interpolate.splev(photons['T_IJ'], RA_tck, der=0)

                Dec_tck = interpolate.splrep(np.arange(0, self.exposure.value + 0.1, 0.1), Dec_._degrees, s=0)
                Dec_Tar = interpolate.splev(photons['T_IJ'], Dec_tck, der=0)
            else:
                t = ts.from_astropy(self.time_init + TimeDelta(photons['T_IJ'].values, format='sec'))
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

            t = ts.from_astropy(
                self.time_init + TimeDelta(photons['T_IJ'].loc[photons['ID'] == 'TARGET'].values, format='sec'))
            RA_, Dec_, dis = self.station2target.at(t).radec()
            photons['RA_TAR'] = 0
            photons['DEC_TAR'] = 0
            photons['RA_TAR'].loc[photons['ID'] == 'TARGET'] = RA_._degrees
            photons['DEC_TAR'].loc[photons['ID'] == 'TARGET'] = Dec_._degrees
            photons['A'] = self.A0 + photons['DRA_CUM']
            photons['D'] = self.D0 + photons['DDEC_CUM']

        elif self.tracking_mode == 'parking':
            t = ts.from_astropy(
                self.time_init + TimeDelta(photons['T_IJ'].loc[photons['ID'] == 'TARGET'].values, format='sec'))
            RA_, Dec_, dis = self.station2target.at(t).radec()
            photons['RA_TAR'] = 0
            photons['DEC_TAR'] = 0
            photons['RA_TAR'].loc[photons['ID'] == 'TARGET'] = RA_._degrees
            photons['DEC_TAR'].loc[photons['ID'] == 'TARGET'] = Dec_._degrees
            photons['A'] = self.A0 + 15 * photons['T_IJ'] / 3600 + photons['DRA_CUM']
            photons['D'] = self.D0 + photons['DDEC_CUM']
        else:
            print('The tracking mode can not be recognized!')

        # Initialize a photon map.
        photon_map = pd.merge(self.stars[['RA', 'DEC', 'ID']], photons, on='ID')

        photon_map['RA'] = np.where(photon_map['ID'] == 'TARGET', photon_map['RA_TAR'], photon_map['RA'])
        photon_map['DEC'] = np.where(photon_map['ID'] == 'TARGET', photon_map['DEC_TAR'], photon_map['DEC'])

        # Calculate the standard coordinates of photons.
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

        # Calculate the measured coordinates of photons.
        xy = np.dot(np.linalg.inv(self.ccd.cd), np.rad2deg(photon_map[['XI', 'ETA']].values.T))
        photon_map['X'] = xy[0, :]
        photon_map['Y'] = xy[1, :]

        photon_map['X_SE'] = photon_map['X'] + rd.normal(0, self.seeing, len(photon_map['X']))
        photon_map['Y_SE'] = photon_map['Y'] + rd.normal(0, self.seeing, len(photon_map['Y']))

        self.photon_map = photon_map

    def rendering(self, bias=0, flat=1,  background=0, dx=0, dy=0):
        """
        Image rendering.
        :param bias: the bias.
        :param flat: the flat.
        :param background: the background.
        :param dx:
        :param dy:
        :return:
        """
        print('rendering ...')

        img_photon = np.zeros(self.ccd.shape)
        x_size = self.ccd.shape[0]
        y_size = self.ccd.shape[1]

        # Image plane sampling.
        for x, y in zip(self.photon_map['X_SE'], self.photon_map['Y_SE']):
            x = x + dx
            y = y + dy
            if (-0.5 < x < x_size) and (-0.5 < y < y_size):
                img_photon[int(y), int(x)] = img_photon[int(y), int(x)] + 1

        self.img_photon = img_photon
        self.image_adu = (img_photon / self.ccd.gain + background)*flat + bias

    def writeto(self, img_photon='photon.fits', img_adu='adu.fits', header=None, **kwargs):
        """
        write the image ti fits file.
        :param img_photon: the file name of image in photon-electrons.
        :param img_adu: the file name of image in ADU.
        :param header: the header of fits file.
        :param kwargs: keywords of fits.PrimaryHDU.writeto.
        :return:
        """
        print('write to ...')

        hdu = fits.PrimaryHDU(self.img_photon.astype(np.float32), header=header)
        hdu.writeto(img_photon, **kwargs)

        hdu = fits.PrimaryHDU(self.image_adu.astype(np.float32), header=header)
        hdu.writeto(img_adu, **kwargs)


if __name__ == '__main__':
    pass
