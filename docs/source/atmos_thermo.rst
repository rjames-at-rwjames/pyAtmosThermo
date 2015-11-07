atmos_thermo package
====================

Submodules
----------

atmos_thermo.atmoslib module
----------------------------

.. automodule:: atmos_thermo.atmoslib
    :members:
    :undoc-members:
    :show-inheritance:

atmos_thermo.constant module
----------------------------

The following constants are accessible as ``atmos_thermo.catmos``. For 
example, ``atmos_thermo.catmos.R_d`` will return the gas constant for 
dry air. All quantities are in SI units.

.. code:: python

    @readonly
    def R_d():
        'gas constant of dry air'
        return 287.04  # , ureg['J / (kg * degK)']

    @readonly
    def R_v():
        'gas constant of water vapor'
        return 461.50  # , ureg['J / (kg * degK)']

    @readonly
    def epsilon():
        'ratio of R_d / R_v'
        return 0.6220  # , ureg['dimensionless']

    @readonly
    def C_vd():
        'heat capacity of dry air at constant volume'
        return 719.0  # , ureg['J / (kg * degK)']

    @readonly
    def C_pd():
        'heat capacity of dry air at constant pressure'
        return 1005.7  # , ureg['J / (kg * degK)']

    @readonly
    def C_vv():
        'heat capacity of water vapor at constant volume'
        return 1410.0  # , ureg['J / (kg * degK)']

    @readonly
    def C_pv():
        'heat capacity of water vapor at constant pressure'
        return 1870.0  # , ureg['J / (kg * degK)']

    @readonly
    def C_l():
        'heat capacity of liquid water'
        return 4190.0  # , ureg['J / (kg * degK)']

    @readonly
    def Lv0():
        'Latent heat of vaporization for water at 0 degrees C'
        return 2.501E6  # , ureg['J / (kg * degK)']

    @readonly
    def L_v():
        'Latent heat of vaporization for water at 100 degrees C'
        return 2.26E6

Similarly, these two are earth constants. For example ``atmos_thermo.cearth.g``.

.. code:: python

    @readonly
    def R_Earth():
        return 6.371E6  # , ureg['m']

    @readonly
    def g():
        return 9.81  # , ureg['m/s/s']
