"""
Canopy Production Module
"""
from __future__ import division
from math import exp
from utils import get_days_in_month
from constants import molPAR_MJ, gDM_mol

def calc_modifier_temp(T_av, T_min, T_max, T_opt):
    if (T_av <= T_min) or (T_av >= T_max):
        res = 0
    else:
        res = ((T_av - T_min) / (T_opt - T_min)) * \
                ((T_max - T_av) / (T_max - T_opt)) ** \
                ((T_max - T_opt) / (T_opt - T_min))
    return res


def calc_modifier_VPD(VPD, CoeffCond):
    res = exp(-1 * CoeffCond * VPD)
    return res


def calc_modifier_soilwater(ASW, MaxASW, SWconst, SWpower):
    moist_ratio = ASW / MaxASW
    res = 1 / (1 + ((1 - moist_ratio) / SWconst) ** SWpower)
    return res


def calc_modifier_soilnutrition(FR, fN0):
    res = fN0 + (1 - fN0) * FR
    return res


def calc_modifier_frost(frost_days, kF):
    res = 1 - kF * (frost_days / 30)
    return res


def calc_modifier_age(stand_age, MaxAge, rAge, nAge):
    rel_age = stand_age / MaxAge
    res = (1 / (1 + (rel_age / rAge) ** nAge))
    return res


def calc_physiological_modifier(modifier_VPD,
        modifier_soilwater, modifier_age):
    # calculate physiological modifier applied to conductance and APARu.
    res = min(modifier_VPD, modifier_soilwater) * modifier_age
    return res


def calc_canopy_cover(stand_age, LAI, fullCanAge, canpower, k):
    # calc canopy cover and light interception.
    canopy_cover = 1.0
    if (fullCanAge > 0) and (stand_age < fullCanAge):
        canopy_cover = (stand_age / fullCanAge) ** canpower
    light_interception = (1 - (exp(-1 * k * LAI)))

    return canopy_cover, light_interception

#added modifier_frost here -Danielle
def calc_canopy_conductance(T_av, LAI, modifier_frost, modifier_physiology,
        TK2, TK3, MaxCond, LAIgcx):
    # calculate canopy conductance from stomatal conductance
    # with added temperature modifier_ Liang Wei
	#added modifier_frost here -Danielle
    canopy_conductance = max(0, min(1, TK2 + TK3 * T_av)) * \
        MaxCond * modifier_frost * modifier_physiology * min(1, LAI / LAIgcx)
    if canopy_conductance == 0:
        canopy_conductance = 0.0001
    return canopy_conductance


def calc_canopy_production(solar_rad, month,
        light_interception, canopy_cover,
        modifier_physiology, modifier_nutrition,
        modifier_temperature, modifier_frost,
        alpha, y):
    # Determine gross and net biomass production
    # Calculate PAR, APAR, APARu and GPP

    RAD = solar_rad * get_days_in_month(month)        # MJ/m^2
    PAR = RAD * molPAR_MJ                      # mol/m^2
    APAR = PAR * light_interception * canopy_cover
    APARu = APAR * modifier_physiology
    alphaC = alpha * modifier_nutrition * modifier_temperature * modifier_frost
    GPPmolc = APARu * alphaC                   # mol/m^2
    GPPdm = (GPPmolc * gDM_mol) / 100          # tDM/ha
    NPP = GPPdm * y                     # tDM/ha - assumes constant respiratory rate

    return PAR, APAR, APARu, GPPmolc, GPPdm, NPP


def canopy_production(T_av, VPD, ASW, frost_days, stand_age,
        LAI, solar_rad, month, CounterforShrub, config):
    
    config_canopy = config.CanopyProduction
    config_shrub = config.ShrubEffect
    config_bio = config.BiomassPartition

    T_min = float(config_canopy.t_min)
    T_max = float(config_canopy.t_max)
    T_opt = float(config_canopy.t_opt)

    CoeffCond = float(config_canopy.coeffcond)

    MaxASW = float(config_canopy.maxasw)
    SWconst = float(config_canopy.swconst0)
    SWpower = float(config_canopy.swpower0)

    FR = float(config_canopy.fr)
    fN0 = float(config_canopy.fn0)

    kF = float(config_canopy.kf)

    MaxAge = float(config_canopy.maxage)
    rAge = float(config_canopy.rage)
    nAge = float(config_canopy.nage)
    
    TK2 = float(config_bio.tk2)
    TK3 = float(config_bio.tk3)
    MaxCond = float(config_bio.maxcond)
    LAIgcx = float(config_bio.laigcx)

    fullCanAge = float(config_canopy.fullcanage)
    canpower = float(config_canopy.canpower)
    k = float(config_canopy.k)

    alpha = float(config_canopy.alpha)
    y = float(config_canopy.y)

    if CounterforShrub is None:
        CounterforShrub = float(config_shrub.counterforshrub)
    KL = float(config_shrub.kl)
    Lsx = float(config_shrub.lsx)

    modifier_temperature = calc_modifier_temp(T_av, T_min, T_max, T_opt)
    modifier_VPD = calc_modifier_VPD(VPD, CoeffCond)
    modifier_soilwater = calc_modifier_soilwater(ASW, MaxASW, SWconst, SWpower)
    modifier_nutrition = calc_modifier_soilnutrition(FR, fN0)
    modifier_frost = calc_modifier_frost(frost_days, kF)
    modifier_age = calc_modifier_age(stand_age, MaxAge, rAge, nAge)
    modifier_physiology = calc_physiological_modifier(modifier_VPD,
            modifier_soilwater, modifier_age)

    canopy_cover, light_interception = calc_canopy_cover(stand_age, LAI, fullCanAge, canpower, k)
    #added modifier_frost here -Danielle
    canopy_conductance = calc_canopy_conductance(T_av, LAI, modifier_frost, modifier_physiology,
        TK2, TK3, MaxCond, LAIgcx)
    PAR, APAR, APARu, GPPmolc, GPPdm, NPP = calc_canopy_production(solar_rad, month,
        light_interception, canopy_cover,
        modifier_physiology, modifier_nutrition,
        modifier_temperature, modifier_frost,
        alpha, y)
        
    modifiers = [modifier_temperature, modifier_VPD,
            modifier_soilwater, modifier_nutrition,
            modifier_frost, modifier_age, modifier_physiology]

    if CounterforShrub == 0:
        LsOpen = LAI * KL
        LsClosed = Lsx * exp(-k * LAI)
        LAIShrub = min(LsOpen, LsClosed)
    elif CounterforShrub == 1:
        LAIShrub = Lsx * exp(-k * LAI)

    if LsClosed <= LsOpen:
        CounterforShrub = 1


    return PAR, APAR, APARu, GPPmolc, GPPdm, NPP, modifiers, LAIShrub, CounterforShrub, canopy_conductance
