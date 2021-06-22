# -*- coding: utf-8 -*-
"""
    config of radioactive sources
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright:  2021 Xu Hangkun <xuhangkun@163.com>
    :license: MIT, see LICENSE for more details.

"""
radioactive_sources_info = {
    "Ge68":{
        "total_gamma_e":0.510998910*2,
        "mean_gamma_e":0.510998910,
        "activity":500,
        "nonlin_calib_time":100,
        "n_gamma":2
    },
    "Mn54":{
        "total_gamma_e":0.835,
        "mean_gamma_e":0.835,
        "activity":50,
        "nonlin_calib_time":3600*10,
        "n_gamma":1
    },
    "Cs137":{
        "total_gamma_e":0.6617,
        "mean_gamma_e":0.6617,
        "activity":50,
        "nonlin_calib_time":3600*10,
        "n_gamma":1
    },
    "K40":{
        "total_gamma_e":1.461,
        "mean_gamma_e":1.461,
        "activity":10,
        "nonlin_calib_time":3600*10,
        "n_gamma":1
    },
    "Co60":{
        "total_gamma_e":1.173237+1.332501,
        "mean_gamma_e":(1.173237+1.332501)/2,
        "activity":10,
        "nonlin_calib_time":3600*10,
        "n_gamma":2
    },
    "amc_n_prompt":{
        "total_gamma_e":6.13,
        "mean_gamma_e":6.13,
        "activity":0.16,
        "nonlin_calib_time":3600*10,
        "n_gamma":1
    },
    "amc_n_delay":{
        "total_gamma_e":2.22,
        "mean_gamma_e":2.22,
        "activity":0.2,
        "nonlin_calib_time":3600*10,
        "n_gamma":1
    },
    "nc_gamma":{
        "total_gamma_e":3.746,
        "mean_gamma_e":3.746,
        "activity":10,
        "nonlin_calib_time":3600,
        "n_gamma":1
    }
}
