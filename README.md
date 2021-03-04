# SPIMT

## SPIMT: Simulating Photometric Images of Moving Targets with Photon-mapping

We present a novel, easy-to-use method based on the photon-mapping technique to simulate photometric images of moving targets. 
Realistic images can be created in two passes: photon tracing and image rendering. 
In the first pass, a photon map is established by emitting discrete photons from the light sources and tracing them through the scene. 
In the second pass, the images are rendered using the information stored in the photon map. 
The nature of light sources, tracking mode of the telescope, point spread function (PSF), and specifications of the CCD are taken into account in the imaging process. 
Photometric images in a variety of observation scenarios can be generated flexibly. 
We compared the simulated images with the observed ones. 
The residuals between them are negligible, and the correlation coefficients between them are high, which means a high fidelity and similarity. 
The method is versatile and can be used to plan future observations of moving targets, interpret existing observations, and provide test images for image processing algorithms.

The simulation architecture was implemented based on some off-the-shelf libraries, such as Astropy, Skyfield, and Astroquery. 
A demo script written in Python code is available here.

## Some packages are required：

+ numpy
+ pandas
+ scipy
+ astropy
+ skyfield
+ astroquery
+ photutils

## Contents

+ spimt.py: the main module of SPIMT.
+ exam_target.py:  a demo script to generate a fits image in target tracking mode.
  + target.fits:  the fits image to mimic in target tracking mode.
  + SIMU_target.fits: the simulated fits image in ADU corresponding to target.fits
  + eSIMU_target.fits: the simulated fits image in photon-electron correspongding to target.fits.
  
+ exam_sidereal.py:  a demo script to generate a fits image in sidereal tracking mode.
  + sidereal.fits:  the fits image to mimic in sidereal tracking mode.
  + SIMU_sidereal.fits: the simulated fits image in ADU corresponding to sidereal.fits
  + eSIMU_sidereal.fits: the simulated fits image in photon-electron correspongding to sidereal.fits.

+ exam_parking.py:  a demo script to generate a fits image in parking tracking mode.
  + parking.fits:  the fits image to mimic in parking tracking mode.
  + SIMU_parking.fits: the simulated fits image in ADU corresponding to parking.fits
  + eSIMU_parking.fits: the simulated fits image in photon-electron correspongding to parking.fits.

## Feedback

Any question can be sent to dujunju@126.com.
