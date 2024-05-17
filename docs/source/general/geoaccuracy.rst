Geolocation Accuracy
====================

Item 4.3 of the CARD4L NRB specification requires, as minimum, an estimate of the absolute location error (ALE) "as
bias and standard deviation, provided in slant range/azimuth, or Northing/Easting" :cite:`ceos_2021`. As desired target the accuracy
is less or equal 0.1 pixels radial root mean square error (rRMSE), which can be defined as:

.. math::
   RMSE_{planar} = \sqrt{RMSE_{SLC,Az}^2 + (\frac{RMSE_{SLC,Rg}}{sin(\theta_{i,min})})^2 + RMSE_{DEM,planar}^2 + RMSE_{proc}^2}

The error induced by the DEM can be described as:

.. math::
   RMSE_{DEM,planar} = \frac{\sigma_{DEM}}{tan(\theta_{i,min})}

where

:math:`\theta_{i,min}` = The minimum possible angle of incidence

:math:`RMSE_{SLC,Az/Rg}` = Error induced by SLC source data in azimuth/range

:math:`RMSE_{DEM,planar}` = Error induced by DEM inaccuracy

:math:`RMSE_{proc}` = Error induced by other processing steps

:math:`\sigma_{DEM}` = DEM accuracy at :math:`1\sigma` (LE68)


Limitations
-----------
Currently, the following simplifications need to be considered for the calculation of rRMSE values found in the metadata
of each S1-NRB product:

* Processing induced errors (:math:`RMSE_{proc}`) and the error term related to DEM interpolation are not further considered and assumed to be 0.
* The DEM accuracy (:math:`\sigma_{DEM}`) is estimated on the global mean accuracy LE90 reported for the COP-DEM :cite:`airbus_2022` under the assumption of gaussian distribution:

    * Global: LE90 = 2.57; LE68 :math:`\approx` 1.56

* rRMSE is only calculated if a COP-DEM was used for processing, otherwise the value is set to ``None``

Development Status
------------------
The development status is tracked and discussed in the following Github issue: https://github.com/SAR-ARD/s1ard/issues/33
