#!/usr/bin/env python

from .constant import atmos as catmos
from .constant import earth as cearth
import numpy


def esat_from_T(T):
    '''Calculates saturation vapor pressure.

    Using the reference Emanuel 1994 (Atmospheric Convection) p.116, equation
    4.4.13 , this function calculates the saturation vapor pressure (esat)
    at the given temperature.

    Args:
        T: Air temperature in Kelvin. Can be numpy array.

    Returns:
        Saturation vapor pressure (esat) in units of Pa.'''

    log_esat = 53.67957 - (6743.769 / T) - (4.8451 * numpy.log(T))
    return numpy.exp(log_esat) * 100.0


def q_from_r(r):
    '''Calculates specific humidity.

    Calculates specific humidity given mixing ratio. Typically mixing ratio
    is r and specific humidity is q in the literature.

    r = [mass of H20 vapor] / [mass of dry air]
    q = [mass of H20 vapor] / [mass of total air]

    relationship is : q = r / (1+r)

    Args:
        r: Mixing ratio in dimensionless ([kg H2O]/[kg Air]) units.

    Returns:
        Specific humidity in dimensionless ([kg H2O]/[kg Air]) units.'''

    return r / (1. + r)


def q_from_e_p(e, p):
    '''Calculates specific humidity.

    Calculates specific humidity given vapor pressure (e) and total pressure
    (p). Typically vapor pressure is e and specific humidity is q in the
    literature.

    q = [mass of H20 vapor] / [mass of total air]

    relationship is : q = epsilon * e / (p - e*(1-epsilon))

    Args:
        e: Vapor pressure in Pascals
        p: Total pressure in Pascals

    Returns:
        Specific humidity in dimensionless ([kg H2O]/[kg Air]) units.'''

    epsilon = catmos.epsilon

    return epsilon * e / (p - e * (1. - epsilon))


def e_from_r_p(r, p):
    '''Calculate partial pressure of water vapor.

    Calculate partial pressure of water vapor (e) from mixing ratio (r) and
    total pressure (p).

    Args:
        r: Mixing ratio in dimensionless ([kg H2O]/[kg Air]) units.
        p: Total pressure in Pascals.
    Returns:
        Vapor pressure (e) in Pascals.
    '''

    return p * r / (r + catmos.epsilon)


def r_from_e_p(e, p):
    '''Calculate mixing ratio of water vapor.

    Calculate mixing ratio of water vapor(r) from vapor pressure (e) and
    total pressure (p).

    Args:
        r: Mixing ratio in dimensionless ([kg H2O]/[kg Air]) units.
        p: Total pressure in Pascals.
    Returns:
        Vapor pressure (e) in Pascals.
    '''
    return catmos.epsilon * e / (p - e)


def relhum_from_e_esat(e, esat):
    '''Calculate the relative humidity of air.

    Calculate the relative humidity of air (H) given the vapor pressure (e)
    and saturation vapor pressure (esat).

    Args:
        e: Vapor pressure in Pascals.
        esat: Saturation vapor pressure in Pascals.
    Returns:
        relative humidity in dimensionless units.'''

    return e / esat


def relhum_from_r_p_T(r, p, T):
    '''Calculate the relative humidity of air.

    Calculate the relative humidity of air (H) given the mixing ratio (r),
    pressure (p) and temperature (T).

    Args:
        r: mixing ratio (dimensionless)
        p: total pressure (Pascals)
        T: temperature (Kelvin)
    Returns:
        relative humidity in dimensionless units.'''

    esat = esat_from_T(T)
    e = e_from_r_p(r, p)

    return e / esat


def theta_from_p_T(p, T, p0=100000.):
    '''Calculate the potential temperature of air.

    Calculate the potential temperature (theta) given pressure (p) and
    temperature (T).

    Formula is : theta = T * (p0 / p) ** (R_d / C_pd)

    Args:
        p: total pressure (Pascals)
        T: temperature (Kelvin)
        p0: reference pressure used in calculating potential temperature
    Returns:
        potential temperature in Kelvin'''

    return T * (p0 / p) ** (catmos.R_d / catmos.C_pd)


def T_from_p_theta(p, theta, p0=100000.):
    '''Calculate the temperature of air.

    Calculate the temperature (T) given pressure (p) and
    potential temperature (theta).

    Formula is : T = theta * (p / p0) ** (R_d / C_pd)

    Args:
        p: total pressure (Pascals)
        theta: potential temperature (Kelvin)
        p0: reference pressure used in calculating potential temperature
    Returns:
        temperature in Kelvin'''

    return theta * (p / p0) ** (catmos.R_d / catmos.C_pd)


def Tv_from_r_T(r, T):
    '''Calculate the virtual temperature.

    Calculate the virtual temperature (T_v) given temperature (T) and mixing
    ratio (r).

    Args:
        r: mixing ratio (dimensionless)
        T: temperature (Kelvin)
    Returns:
        Virtual Temperature in Kelvin
    '''

    return T * (1. + r / catmos.epsilon) / (1. + r)


def thetav_from_r_p_T(r, p, T):
    '''Calculate the virtual potential temperature.

    Calculate the virtual potential temperature (thetav) given temperature (T)
    and mixing ratio (r).

    Args:
        r: mixing ratio (dimensionless)
        p: total pressure (Pascals)
        T: temperature (Kelvin)
    Returns:
        Virtual Temperature in Kelvin
    '''

    return theta_from_p_T(p, Tv_from_r_T(r, T))


def pd_from_r_p(r, p):
    '''Calculate dry pressure without water vapor.

    Calculate dry pressure (pd) from mixing ratio (r) and total pressure (p)

    Args:
        r: Mixing ratio in dimensionless ([kg H2O]/[kg Air]) units.
        p: Total pressure in Pascals.
    Returns:
        Dry pressure (pd) in Pascals.
    '''
    return p / (1 + r / catmos.epsilon)


def moist_specific_entropy_from_r_p_T(r, p, T):
    '''Calculate the moist specific entropy.

    Calculate moist specific entropy for the mixing ratio (r), total pressure
    (p), and temperature (T).

    Uses the definition given in Benedict et. al 2014 - Gross Moist Stability
    and MJO Simulation Skill ... Equation 1

    Args:
        r: mixing ratio (dimensionless)
        p: total pressure (Pascals)
        T: temperature (Kelvin)
    Returns:
        Moist specific entropy in J / (kg K)
    '''
    Lv = catmos.Lv0
    Cpd = catmos.C_pd
    Cpv = catmos.C_pv
    R_d = catmos.R_d
    R_v = catmos.R_v

    T_R = 273.1
    e_sf = 611.
    P_R = 100000.

    e = e_from_r_p(r, p) #vapor pressure of water vapor
    p_dry = p - e # dry vapor pressure

    term0 = (Cpd + r * Cpv) * numpy.log(T/T_R)
    term1 = R_d * numpy.log(p_dry/P_R)
    term2 = r * R_v * numpy.log(e/e_sf)
    term3 = Lv * r/T_R

    return (term0 - term1 - term2 + term3)


def moist_specific_entropy_Emanuel_from_r_p_T(r, p, T,
                                              r_rain=0, r_cloud=0, r_ice=0):
    '''Calculate the moist specific entropy.

    Using the definition of Emanuel 1994 (Atmospheric Convection).
    Equation 4.5.9

    s =   (C_pd + C_pv * r_t) * ln(T)
        - R_d ln (pd)
        + Lv r / T
        - r R_v ln (H)

    Args:
        r: mixing ratio (dimensionless)
        p: total pressure (Pascals)
        T: temperature (Kelvin)
    Returns:
        Equivalent potential temperature in Kelvin
    '''
    r_t = r + r_rain + r_cloud + r_ice  # total mixing ratio
    Lv = Lv_from_T(T)

    esat = esat_from_T(T)  # saturation vapor pressure
    e = e_from_r_p(r, p)  # vapor pressure
    pd = p - e  # dry pressure
    C_p_total = (catmos.C_pd + catmos.C_pv * r_t)

    term0 = C_p_total * numpy.log(T)
    term1 = catmos.R_d * numpy.log(pd)
    term2 = Lv * r / T
    term3 = r * catmos.R_v * numpy.log(e / esat)

    return (term0 - term1 + term2 - term3)


def moist_specific_entropy_Emanuel_sat_from_p_T(p, T, r_rain=0,
                                                r_cloud=0, r_ice=0):
    '''Calculate the moist specific entropy.

    Using the definition of Emanuel 1994 (Atmospheric Convection).
    Equation 4.5.9

    s =   (C_pd + C_pv * r_t) * ln(T)
        - R_d ln (pd)
        + Lv r / T
        - r R_v ln (H)

    Args:
        p: total pressure (Pascals)
        T: temperature (Kelvin)
    Returns:
        Equivalent potential temperature in Kelvin
    '''
    esat = esat_from_T(T)  # saturation vapor pressure
    r = r_from_e_p(esat, p)

    r_t = r + r_rain + r_cloud + r_ice  # total mixing ratio
    Lv = Lv_from_T(T)

    pd = p - esat  # dry pressure
    C_p_total = (catmos.C_pd + catmos.C_pv * r_t)

    term0 = C_p_total * numpy.log(T)
    term1 = catmos.R_d * numpy.log(pd)
    term2 = Lv * r / T

    return (term0 - term1 + term2)


def Lv_from_T(T):
    '''Calculate latent heat of vaporization for water, Lv.

    Args:
        T: temperature (Kelvin)
    Returns:
        Latent heat of vaporization [J / kg K]
    '''
    return catmos.Lv0 + (catmos.C_pv - catmos.C_pv) * (T - 273.15)


def thetae_Bolton_from_r_p_T(r, p, T, p0=100000, r_rain=0., r_cloud=0.
                             ):
    '''Calculate the pseudoadiabatic equivalent potential temperature.

    Using the Bolton 1980 definition.
    
    Args:
        r: mixing ratio (dimensionless)
        p: total pressure (Pascals)
        T: temperature (Kelvin)
    Returns:
        Pseudoadibatic equivalent potential temperature in Kelvin'''
    e = e_from_r_p(r, p)
    T_star = 2840. / (3.5 * numpy.log(T) - numpy.log(e) - 4.805) + 55.

    theta_ep = T * (p0 / p) ** (0.2854 * (1 - 0.28 * r))
    theta_ep *= numpy.exp(r * (1. + 0.81 * r) * (3376. / T_star - 2.54))

    return theta_ep


def thetae_from_r_p_T(r, p, T, p0=100000, r_cloud=0.):
    '''Calculate the reversible equivalent potential temperature.

    Using the definition of Emanuel 1994 (Atmospheric Convection).
    Equation 4.5.11

    theta_e = T * [(p0 / pd) ** (R_d / (C_pd + C_pv * r_t))]
                * [H ** (-r R_v / (C_pd + C_pv * r_t)) ]
                * exp[L_v * r / ((C_pd + C_pv * r_t) * T)]

    Args:
        r: mixing ratio (dimensionless)
        p: total pressure (Pascals)
        T: temperature (Kelvin)
    Returns:
        Reversible equivalent potential temperature in Kelvin
    '''
    r_t = r + r_cloud
    Lv = catmos.Lv0

    esat = esat_from_T(T)  # saturation vapor pressure
    e = e_from_r_p(r, p)  # vapor pressure
    pd = p - e  # dry pressure
    C_p_total = (catmos.C_pd + catmos.C_l * r_t)

    pot_temp = T * (p0 / pd) ** (catmos.R_d / C_p_total)
    humid_term = (e / esat) ** (- r * (catmos.R_v / C_p_total))
    latent_term = numpy.exp(Lv * r / (C_p_total * T))

    return pot_temp * humid_term * latent_term


def thetae_star_from_p_T(p, T):
    '''Calculate the saturation equivalent potential temperature.

    Using the definition of Emanuel 1994 (Atmospheric Convection).
    Equation 4.5.11

    thetaesat = T * [(p0 / pd) ** (R_d / (C_pd + C_pv * r_t))]
                * exp[L_v * r / ((C_pd + C_pv * r_t) * T)]

    Args:
        r: mixing ratio (dimensionless)
        p: total pressure (Pascals)
        T: temperature (Kelvin)
    Returns:
        Saturation equivalent potential temperature in Kelvin
    '''
    esat = esat_from_T(T)
    r = r_from_e_p(esat, p)

    return thetae_from_r_p_T(r, p, T)


def h_from_r_T_z(r, T, z):
    '''Calculate the moist static energy (h).

    Using the classical definiton of moist static energy
    h = C_pd T + Lv r + g z

    Args:
        r: mixing ratio (dimensionless, or kg/kg)
        T: temperature (Kelvin)
        z: height (meters)
    Returns:
        Moist static energy in J/kg (m^2 / s^2)
    '''

    h = (catmos.C_pd) * T + catmos.Lv0 * r + cearth.g * z
    return h

def h_Em1994_from_r_T_z(r, T, z):
    '''Calculate the moist static energy (h).

    Using the definition of Emanuel 1994 (Atmospheric Convection).
    Equation 4.5.23

    h = (C_pd + r_t C_pv) T + Lv r + g z

    Args:
        r: mixing ratio (dimensionless)
        T: temperature (Kelvin)
        z: height (meters)
    Returns:
        Moist static energy in J/kg (m^2 / s^2)
    '''
    r_t = r
    Lv = Lv_from_T(T)
    # Lv = catmos.Lv0

    h = (catmos.C_pd + r_t * catmos.C_pv) * T + Lv * r + (1+r_t) * cearth.g * z
    return h


def hstar_Em1994_from_p_T_z(p, T, z):
    '''Calculate the saturation moist static energy (h).

    Using the definition of Emanuel 1994 (Atmospheric Convection).
    Equation 4.5.23

    h = (C_pd + r_sat C_pv) T + Lv (r_sat) + g z

    Args:
        p: total pressure (Pascals)
        T: temperature (Kelvin)
        z: height (meters)
    Returns:
        Saturation moist static energy in J/kg (m^2 / s^2)
    '''

    esat = esat_from_T(T)
    rsat = r_from_e_p(esat, p)
    r = rsat

    return h_Em1994_from_r_T_z(r, T, z)


def hstar_classical_from_p_T_z(p, T, z):
    '''Calculate the saturation moist static energy (h).

    Using the definition of Emanuel 1994 (Atmospheric Convection).
    Equation 4.5.23

    h = (C_pd + r_sat C_pv) T + Lv (r_sat) + g z

    Args:
        p: total pressure (Pascals)
        T: temperature (Kelvin)
        z: height (meters)
    Returns:
        Saturation moist static energy in J/kg (m^2 / s^2)
    '''

    esat = esat_from_T(T)
    rsat = r_from_e_p(esat, p)
    r = rsat

    return h_from_r_T_z(r, T, z)


def h_KF_from_r_T_z_rliq_rfroz(r, T, z, rliq, rfroz):
    """Calculate moist static energy according to the Kain-Fritsch definition.

    The standard frozen MSE formulation is:
        h = cp * T + Lv * r + g * z - Lf * r_froz

    Kain-Fritsch uses:
        cp_moist = cp_dry * (1 + 0.887 * r) # this is only the vapor r

        Lv = XLV0 - T * XLV1
        Lf = XLS0 - T * XLS1

        XLV0 = 3.15E6 J/kg
        XLV1 = 2370 J/(kg K)
        XLS0 = 2.905E6 J/kg
        XLS1 = 259.532 J/(kg K)
    """
    XLV0 = 3.15E6  # J/kg
    XLV1 = 2370  # J/(kg K)
    XLS0 = 2.905E6  # J/kg
    XLS1 = 259.532  # J/(kg K)

    cp_moist = (catmos.C_pd + (catmos.C_pd * 0.887) * r)
    Lv = XLV0 - T * XLV1
    Ls = XLS0 - T * XLS1

    h = cp_moist * T
    h += Lv * r
    h += (1.0 + rliq + rfroz) * cearth.g * z
    h += - Ls * rfroz

    return h


def hd_from_r_T_z(r, T, z):
    '''Calculate the classical dry static energy (hd).

    Args:
        p: total pressure (Pascals)
        T: temperature (Kelvin)
        z: height (meters).
    Returns:
        Dry static energy in J/kg (m^2 / s^2)
    '''

    hd = (catmos.C_pd) * T + cearth.g * z
    return hd


def hd_Em1994_from_r_T_z(r, T, z):
    '''Calculate the dry static energy (hd).

    Using the definition of Emanuel 1994 (Atmospheric Convection).
    Equation 4.5.23. Notice that the mixing ratio is required to
    calculate the heat capacity correctly.

    hd = (C_pd + r_t C_pv) T + g z

    Args:
        r: mixing ratio (dimensionless)
        p: total pressure (Pascals)
        T: temperature (Kelvin)
        z: height (meters).
    Returns:
        Dry static energy in J/kg (m^2 / s^2)
    '''

    hd = (catmos.C_pd + r * catmos.C_pv) * T + (1 + r) * cearth.g * z
    # hd = catmos.C_pd * T + cearth.g * z
    return hd


def streamfunction_from_v_p(v, p):
    "Calculate streamfunction from v, p"
    sfn = numpy.zeros(p.shape, dtype=p.dtype)

    dp = numpy.diff(p[..., ::-1, :, :], axis=-3)
    v_avg = 0.5 * \
        (v[..., 1:, :, :] + v[..., :-1, :, :])[..., ::-1, :, :]

    sfn[..., :-1, :, :] = numpy.cumsum(dp * v_avg, axis=-3)[..., ::-1, :, :]

    return sfn / cearth.g
