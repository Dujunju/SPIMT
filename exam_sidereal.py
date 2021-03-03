import spimt as psim
from astropy.stats import SigmaClip
import numpy as np
import astropy.io.fits as fits
from photutils import Background2D, SExtractorBackground
import numpy.random as rd


# input file, a file to mimic.
fits_file = 'exam_sidereal.fits'
img_photon = 'eSIMU_'+fits_file
img_adu = 'SIMU_'+fits_file

# the parameter of input file
data = fits.getdata(fits_file)
header = fits.getheader(fits_file)
background = Background2D(
    data, box_size=(32, 32),
    filter_size=(3, 3),
    sigma_clip=SigmaClip(sigma=3.0),
    bkg_estimator=SExtractorBackground())

# the simulation parameters.
Ps = {
    'A0': 35.66464205398764,
    'D0': 43.04244009345635,
    'latitude': 37.53591667,
    'longitude': 122.04961111,
    'elevation': 100,
    'exposure': 4,
    'line1': '1 00000U 15083A   21028.91508324 -.00000329  00000-0  00000-0 0  9996',
    'line2': '2 00000   0.1408  93.4305 0000435  96.7329  13.2856  1.00270407 18786',
    'magnitude': 99,
    'mag_range': (-10, 20),
    'tracking_mode': 'sidereal',
    'ra_sigma': 0,
    'dec_sigma': 0,
    'shape': (1024, 1024),
    'fov': 0.4,
    'cd': np.array([[-6.191165792747975E-06, 0.0001963007872101023], [-0.0001963564892993812, -6.169550777235073E-06]]),
    'dx': 512-0.5,
    'dy': 512-0.5,
    'eA': 0,
    'eD': 0,
    'zero_point': 2.185,
    'K': 0,
    'gain': 1.8,
    'inter': True,
    'plate_const': [],
    'pixel_size': 0.7066,
    'pixel_scale': 1.0,
    'seeing': 1.2,
    'time_init': '2021-01-11T11:28:15',
    'bias': 0,
    'flat': 1,
    'back': np.random.normal(background.background, background.background_rms, (1024, 1024)),
    'img_photo': img_photon,
    'img_adu': img_adu}


# Initialize a scene.
a_fits = psim.SynFits()

# Add a ground station.
a_fits.add_station(
    latitude_degrees=Ps['latitude'],
    longitude_degrees=Ps['longitude'],
    elevation_m=Ps['elevation'])

# Add a moving target.
a_fits.add_target(
    Ps['line1'],
    Ps['line2'],
    magnitude=Ps['magnitude'])

# Add a telescope.
a_fits.add_telescope(
    ra_sigma=Ps['ra_sigma'],
    dec_sigma=Ps['dec_sigma'],
    fov=Ps['fov'],
    zero_point=Ps['zero_point'],
    k=Ps['K'])

# Add a CCD.
a_fits.add_ccd(
    shape=Ps['shape'],
    gain=Ps['gain'],
    plate_const=Ps['plate_const'],
    cd=Ps['cd'],
    pixel_size=Ps['pixel_size'],
    pixel_scale=Ps['pixel_scale'])

# Set the observation setting.
a_fits.set_setting(
    seeing=Ps['seeing'],
    time_init=Ps['time_init'],
    exposure=Ps['exposure'],
    tracking_mode=Ps['tracking_mode'],
    A0=Ps['A0'],
    D0=Ps['D0'],
    eA=Ps['eA'],
    eD=Ps['eD'])

# Add field stars.
a_fits.add_stars(mag_range=Ps['mag_range'])

# Add photons.
a_fits.add_photon()

# Photon tracing.
a_fits.tracking(inter=Ps['inter'])

# Image rendering.
a_fits.rendering(
    bias=Ps['bias'],
    background=Ps['back'],
    dx=Ps['dx'],
    dy=Ps['dy'])

# Write image to fits file.
a_fits.writeto(
    img_photon=Ps['img_photo'],
    img_adu=Ps['img_adu'],
    overwrite=True)
