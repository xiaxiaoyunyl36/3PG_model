from __future__ import division

import re
from math import pi

import numpy as np

from framework import Model, BookKepper
from utils import get_stand_age, get_day_length

from CanopyProduction import canopy_production
from BiomassPartition import biomass_partition
from WaterBalance import water_balance
from StemMortality import stem_mortality, calc_factors_age


mapper = {"stand_age": "stand_age",
        "lai": "LAI",
        "mai": "MAI",
        "basarea": "BasArea",
        "height": "Height",
        "d13ctissue": "D13CTissue",
        "modifier_physiology": "modifier_physiology",
        "npp": "NPP",
        "asw": "ASW",
        "transp": "transp",
        "loss_water": "loss_water",
        "standvol": "StandVol",
        "stemno": "StemNo",
        "par": "PAR",
        "intercippm": "InterCiPPM",
        "wf": "WF",
        "ws": "WS",
        "wr": "WR",
        "avstemmass": "AvStemMass",
        "delwf": "delWF", 
        "delwr": "delWR", 
        "delws": "delWS",
        "d18oleaf": "d18Oleaf", 
        "d18ocell": "d18Ocell",
        "d18ocell_peclet": "d18Ocell_peclet",
        "avdbh": "avDBH",
        "canopy_conductance": "canopy_conductance",
        "canopy_transpiration_sec": "canopy_transpiration_sec",
        "l": "l",
        "gppdm": "GPPdm",
        "totallitter": "TotalLitter"}

class Model3PG(Model):
    def __init__(self, fpath_setting):
        super(Model3PG, self).__init__(fpath_setting)
        self.initialize()

    def initialize(self):
        fpath_input = self.config.IO.input
        fpath_output = self.config.IO.output

        self.data = np.loadtxt(fpath_input, skiprows=1)
        self.keeper = BookKepper(fpath_output)
        self.keeper.initialize(self.config.Output)

    def teardown(self):
        self.data = None
        self.keeper.shutdown()

    def run(self):
        config_time = self.config.TimeRange
        config_site = self.config.SiteCharacteristics
        config_initial = self.config.InitialState
        config_stem = self.config.StemMortality

        lat = float(config_site.lat)
        EndYear = int(config_time. endyear)
        InitialYear = int(config_time.initialyear)
        InitialMonth = int(config_time.initialmonth)
        YearPlanted = int(config_time.yearplanted)
        MonthPlanted = int(config_time.monthplanted)
        EndAge = int(config_time.endage)

        nYears = EndYear - InitialYear + 1
        
        # Assign initial state of stand
        stand_age, StartAge, \
            InitialYear, InitialMonth, MonthPlanted = get_stand_age(config_site.lat,
                        InitialYear, InitialMonth,
                        YearPlanted, MonthPlanted, EndAge)

        # do annual calculation
        metMonth = InitialMonth
        for year in range(StartAge, EndAge + 1):
            print('year', year)

            # do monthly calculations
            month = InitialMonth
            for month_counter in range(1, 12 + 1):
                if (year == 0) and (month == InitialMonth):
                    WS = float(config_initial.initialws)
                    WF = float(config_initial.initialwf)
                    WR = float(config_initial.initialwr)
                    StemNo = float(config_initial.initialstocking)
                    ASW = float(config_initial.initialasw)
                    TotalLitter = 0
                    # thinEventNo = 1
                    # defoltnEventNo = 1
                    irrig = 0 # TODO

                    SLA0 = float(config_stem.sla0)
                    SLA1 = float(config_stem.sla1)
                    tSLA = float(config_stem.tsla)
                    fracBB0 = float(config_stem.fracbb0)
                    fracBB1 = float(config_stem.fracbb1)
                    tBB = float(config_stem.tbb)
                    StemConst = float(config_stem.stemconst)
                    StemPower = float(config_stem.stempower)
                    Density = float(config_stem.density)

                    SLA, fracBB = calc_factors_age(stand_age, SLA0,
                            SLA1, tSLA, fracBB0, fracBB1, tBB)
                    AvStemMass = WS * 1000 / StemNo                 # kg/tree
                    avDBH = (AvStemMass / StemConst) ** (1 / StemPower)
                    BasArea = (((avDBH / 200) ** 2) * pi) * StemNo
                    LAI = WF * SLA * 0.1
                    StandVol = WS * (1 - fracBB) / Density

                    if stand_age > 0:
                        MAI = StandVol / stand_age
                    else:
                        MAI = 0

                    Height = D13CTissue = NPP = InterCiPPM = delWF = delWR = delWS = 0
                    d18Oleaf = d18Ocell = d18Ocell_peclet = canopy_conductance = GPPdm = 0
                    transp = loss_water = canopy_transpiration_sec = l = TotalLitter = 0
                    modifiers = 7 * [0]
                    delStemNo = 0
                    modifier_physiology = 0
                    PAR = 0
                else:
                    #print 'month', month
                    config_site = self.config.SiteCharacteristics

                    lat = float(config_site.lat)
                    elev = float(config_site.elev)
                    # assign meteorological data at this month
                    if month >= 12:
                        month = 1
                    # if metMonth > 12 * mYears:
                    #     metMonth = 1

                    # T_max = self.data[metMonth, 0] #CJS note: does not need Tmax met. data: VPD and SRAD are already in inputs
                    # T_min = self.data[metMonth, 1] #CJS note: does not need Tmax met. data: VPD and SRAD are already in inputs
                    T_av = self.data[metMonth, 2]
                    # VPD = get_VPD(T_min, T_max) #CJS note: does not need VPD met. data: VPD data are already in inputs
                    VPD = self.data[metMonth, 3]
                    rain = self.data[metMonth, 4]
                    solar_rad = self.data[metMonth, 5]
                    # rain_days = int(self.data[metMonth, 6])
                    day_length = get_day_length(lat, month)
                    frost_days = int(self.data[metMonth, 7])
                    CaMonthly = self.data[metMonth, 8]
                    D13Catm = self.data[metMonth, 9]
                    d18Osrc = self.data[metMonth, 10]

                    CounterforShrub = None

                    # Canopy Production Module
                    PAR, APAR, APARu, \
                        GPPmolc, GPPdm, NPP, \
                        modifiers, LAIShrub, \
                        CounterforShrub, canopy_conductance = canopy_production(T_av, VPD,
                                    ASW, frost_days, stand_age,
                                    LAI, solar_rad, month, CounterforShrub, self.config)

					# Water Balance Module
                    transpall, transp, transpshrub, loss_water, ASW, \
                        monthlyIrrig, canopy_transpiration_sec = water_balance(solar_rad, VPD,
                                day_length, LAI, rain, irrig,
                                month, ASW, canopy_conductance, LAIShrub, self.config)

                    # Biomass Partion Module
                    modifier_physiology = modifiers[-1]
                    WF, WR, WS, TotalW, TotalLitter, \
                        D13CTissue, InterCiPPM, \
                        delWF, delWR, delWS, d18Oleaf, d18Ocell, \
                        d18Ocell_peclet = biomass_partition(T_av, LAI,
                                elev, CaMonthly, D13Catm,
                                WF, WR, WS, TotalLitter,
                                NPP, GPPmolc, stand_age, month, avDBH,
                                modifier_physiology, VPD, d18Osrc,  
                                canopy_conductance, canopy_transpiration_sec, self.config)

                    # Stem Mortality Module
                    stand_age, LAI, MAI, \
                        avDBH, BasArea, Height, \
                        StemNo, delStemNo, StandVol, \
                        WF, WR, WS, AvStemMass = stem_mortality(WF, WR, WS, StemNo, delStemNo,
                                    stand_age, self.config)

                self.keeper.keep(mapper, locals())

                metMonth = metMonth + 1
                month = month + 1
            # break


if __name__ == '__main__':
    # fpath_test = r'../test/Test_config.cfg'
    fpath_test = r'/Users/admin/workspace/3PG_python/test/Test_config.cfg'
    model = Model3PG(fpath_test)
    print(model.data)

    model.run()
