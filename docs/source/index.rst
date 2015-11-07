pyAtmosThermo Documentation
===========================
This library allows calculation of moist atmospheric thermodynamic quantities.

For example, the following snippet will calculate moist static energy from 
mixing ratio of 0.001 kg/kg, 298 Kelvin, and 200 meters height. Thus, MSE has 
the value 304161.6 J/kg.

.. code:: python

	import atmos_thermo
	h = atmos_thermo.h_from_r_T_z(0.001, 298, 200)

This code calculates the saturation vapor pressure and produces 3533.7 Pa.

.. code:: python

	e_saturation = atmos_thermo.esat_from_T(300)

All functions are designed around using numpy arrays transparently, and using
numpy arrays is vastly faster than passing in values one by one. 

Contents:

.. toctree::

   atmos_thermo


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

