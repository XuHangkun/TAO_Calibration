"""
    Evaluate model
"""
import numpy as np
from math import acos

def diff_mean(nonuni_func_1,nonuni_func_2,radius_cut=650):
    """
    args:
        nonuni_func : scipy function describe nonuniformity
    return:
        mean(V) of the map with radius_cut
    """
    radius = np.power(np.linspace(1,pow(radius_cut,3),300),1.0/3)
    theta = np.array([acos(x) for x in np.linspace(-0.999,0.999,200)])*180./np.pi
    radius,theta = np.meshgrid(radius,theta)
    value = (nonuni_func_1(radius,theta) - nonuni_func_2(radius,theta))/nonuni_func_2(radius,theta)
    return np.mean(value)

def diff_rms(nonuni_func_1,nonuni_func_2,radius_cut=650):
    """
    args:
        nonuni_func : scipy function describe nonuniformity
    return:
        mean(V) of the map with radius_cut
    """
    radius = np.power(np.linspace(1,pow(radius_cut,3),300),1.0/3)
    theta = np.array([acos(x) for x in np.linspace(-0.999,0.999,200)])*180./np.pi
    radius,theta = np.meshgrid(radius,theta)
    value = (nonuni_func_1(radius,theta) - nonuni_func_2(radius,theta))/nonuni_func_2(radius,theta)
    return np.std(value)
