Geolocation Accuracy
====================

Item 4.3 of the CARD4L NRB specification [1] requires, as minimum, an estimate of the absolute location error (ALE) “as
bias and standard deviation, provided in slant range/azimuth, or Northing/Easting” [RD-2]. As desired target theaccuracy
is less or equal 0.1 pixels radial root mean [square error (rRMSE), which can be defined as:

.. math::
   RMSE_planar = \sqrt{RMSE_{SLC,Az}^2 + (\frac{RMSE_{SLC,Rg}}{sin(\theta_{i,min})})^2 + RMSE_{DEM,planar}^2 + RMSE_{proc}^2}

The error induced by the DEM can be described as:

.. math::
   RMSE_DEM,planar = \frac{\sigma_DEM}{tan(\theta_i,min)}

where

:math:`\theta_{i,min}` = The minimum possible angle of incidence

:math:`RMSE_{SLC,Az/Rg}` = Error induced by SLC source data in azimuth/range

:math:`RMSE_{DEM,planar}` = Error induced by DEM inaccuracy

:math:`RMSE_proc` = Error induced by other processing steps

:math:`\sigma_DEM` = DEM accuracy at :math:`1\sigma` (LE68)


Currently, the following simplifications need to be considered for the calculation of rRMSE values for the S1-NRB
product:

- Processing induced errors (:math:`RMSE_proc`) and the error term related to DEM interpolation are not further considered and assumed to be 0.
- The DEM accuracy (:math:`\sigma_DEM`) is estimated on the global mean accuracy LE90 reported for the COP-DEM [2] under the assumption of gaussian distribution. The actual S1-NRB products make use of the per-tile LE68 accuracy, which is provided in the metadata of each DEM tile.


[1] https://ceos.org/ard/files/PFS/NRB/v5.5/CARD4L-PFS_NRB_v5.5.pdf
[2] https://spacedata.copernicus.eu/documents/20126/0/GEO1988-CopernicusDEM-SPE-002_ProductHandbook_I1.00.pdf/082dd479-f908-bf42-51bf-4c0053129f7c?t=1586526993604
