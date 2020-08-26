from __future__ import division
from math import exp, log, pi

"""
Stem Mortality Module
"""

# def doThinning(n As Integer, table As Variant)
# function for manually thinning
# not implemented yet


# def doDefoliation(n As Integer, table As Variant)
# function for defoliation
# not implemented yet



# Update tree and stand data at the end of this time period,
# taking mortality, thinning or defoliation into account

"""
    # Perform any thinning or defoliation events for this time period
    if thinEventNo <= nThinning:
        doThinning(thinEventNo, Thinning)
    if defoltnEventNo <= nDefoliation:
        doDefoliation(defoltnEventNo, Defoliation)
"""

def getMortality(oldN, oldW,
        mS, wSx1000, thinPower):
    """
    Input:
        oldN, Double
        oldW, Double
    Output:
        mortality rate, Double
    Description:
        This function determines the number of stems to remove to ensure the
        self-thinning rule is satisfied. It applies the Newton-Rhapson method
        to solve for N to an accuracy of 1 stem or less. To change this,
        change the value of "accuracy".

        This was the old mortality function:
          getMortality = oldN - 1000 * (wSx1000 * oldN / oldW / 1000) ^ (1 / thinPower)
        which has been superceded by the following ...
    """
    accuracy = 1 / 1000
    n = oldN / 1000
    x1 = 1000 * mS * oldW / oldN
    i = 0
    while True:
        i = i + 1
        x2 = wSx1000 * (n ** (1 - thinPower))
        fN = x2 - x1 * n - (1 - mS) * oldW
        dfN = (1 - thinPower) * x2 / n - x1
        dN = -1 * fN / dfN
        n = n + dN
        if (abs(dN) <= accuracy) or (i >= 5):
            break
    res = oldN - 1000 * n
    return int(res)


def calc_mortality(WF, WR, WS, StemNo, delStemNo,
        wSx1000, thinPower, mF, mR, mS):
    # Calculate mortality
    wSmax = wSx1000 * (1000 / StemNo) ** thinPower
    AvStemMass = WS * 1000 / StemNo
    delStems = 0
    if wSmax < AvStemMass:
        delStems = getMortality(StemNo, WS,
              mS, wSx1000, thinPower)
        WF = WF - mF * delStems * (WF / StemNo)
        WR = WR - mR * delStems * (WR / StemNo)
        WS = WS - mS * delStems * (WS / StemNo)
        # wSmax = wSx1000 * (1000 / StemNo) ** thinPower
    StemNo = StemNo - delStems
    AvStemMass = WS * 1000 / StemNo
    delStemNo = delStemNo + delStems
    return WF, WR, WS, AvStemMass, StemNo, delStemNo


def calc_factors_age(stand_age, SLA0, SLA1, tSLA,
        fracBB0, fracBB1, tBB):
    # update age-dependent factors
    SLA = SLA1 + (SLA0 - SLA1) * exp(-log(2) * (stand_age / tSLA) ** 2)
    fracBB = fracBB1 + (fracBB0 - fracBB1) * exp(-log(2) * (stand_age / tBB))
    return SLA, fracBB


def update_stands(stand_age, WF, WS, AvStemMass, StemNo,
        SLA, fracBB, StemConst, StemPower, Density, HtC0, HtC1):
    # update stsand characteristics
    LAI = WF * SLA * 0.1
    avDBH = (AvStemMass / StemConst) ** (1 / StemPower) #CJS note must be in [cm]
    BasArea = (((avDBH / 200) ** 2) * pi) * StemNo #CJS must be in [m]
    StandVol = WS * (1 - fracBB) / Density
    if stand_age > 0:
        MAI = StandVol / stand_age
    else:
        MAI = 0

    # Height equation (Wykoff 1982) is in English unit,
    # DBH is first convert to inch.
    # Finally Ht is convert form feet to meters._Liang
    Height = (exp(HtC0 + HtC1 / (avDBH / 2.54 + 1)) + 4.5) * 0.3048

    return LAI, MAI, avDBH, BasArea, Height, StandVol


def stem_mortality(WF, WR, WS,
        StemNo, delStemNo, stand_age, config, doThinning=None, doDefoliation=None):

    config_stem = config.StemMortality

    wSx1000 = float(config_stem.wsx1000)
    thinPower = float(config_stem.thinpower)
    mF = float(config_stem.mf)
    mR = float(config_stem.mr)
    mS = float(config_stem.ms)

    SLA0 = float(config_stem.sla0)
    SLA1 = float(config_stem.sla1)
    tSLA = float(config_stem.tsla)
    fracBB0 = float(config_stem.fracbb0)
    fracBB1 = float(config_stem.fracbb1)
    tBB = float(config_stem.tbb)

    StemConst = float(config_stem.stemconst)
    StemPower = float(config_stem.stempower)
    Density = float(config_stem.density)
    HtC0 = float(config_stem.htc0)
    HtC1 = float(config_stem.htc1)

    if doThinning is not None:
        doThinning()
    if doDefoliation is not None:
        doDefoliation()

    stand_age = stand_age + 1.0 / 12

    WF, WR, WS, AvStemMass, StemNo, delStemNo = calc_mortality(WF, WR, WS, StemNo, delStemNo,
        wSx1000, thinPower, mF, mR, mS)
    SLA, fracBB = calc_factors_age(stand_age, SLA0, SLA1, tSLA,
        fracBB0, fracBB1, tBB)
    LAI, MAI, avDBH, BasArea, Height, StandVol = update_stands(stand_age, WF, WS, AvStemMass, StemNo,
        SLA, fracBB, StemConst, StemPower, Density, HtC0, HtC1)
    return stand_age, LAI, MAI, avDBH, BasArea, Height, StemNo, delStemNo, StandVol, WF, WR, WS, AvStemMass

