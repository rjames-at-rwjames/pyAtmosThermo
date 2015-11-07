# This block declares constants that cannot be changed. Useful for things
# like acceleration due to gravity or Earth's pressure.


def readonly(f):
    def fset(self, value):
        raise SyntaxError

    def fget(self):
        return f()
    return property(fget, fset)


class _atmos_const(object):

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


class _earth_const(object):

    @readonly
    def R_Earth():
        return 6.371E6  # , ureg['m']

    @readonly
    def g():
        return 9.81  # , ureg['m/s/s']

atmos = _atmos_const()
earth = _earth_const()
