# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 15:03:05 2019

@author: MSI
"""



# basic moudle
import numpy as np 
import pandas as pd
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import math as mt
import random
import time

# seaborn moudle
import seaborn as sns

# cmatch moudle
import cmath as cmath



def ReflectionPP(a1, a2, b1, b2, p1, p2, a_min=0.0, a_max=0.0, N=1):
    # Takes the p-wave velocity, s-wave velocity and density for two layers and plots
    # the normalized energy coefficients of the reflected and transmitted p and s waves
    
    # Compute the aucoustic Impedance
    
    # a1: vp1
	# a2: vp2
	# b1: vs1
	# b2: vs2
	# p1: dense1
	# p2: dense2
    # a_min: 初始角度
	# a_max: 终止角度
    # N: 初始角度和终止角度之间采样点的个数
    
    # Create array of incidence angles between 0 and 90 degrees
    # N = 15
    theta1deg = np.linspace(a_min, a_max, N)
    theta1 = np.radians(theta1deg)
    
    #make zero arrays to be filled in later
    Rpp = np.zeros(N)
    
    for i in range(0,N):
        t1 = theta1[i]
        t11= np.arcsin(a2/a1*np.sin(t1))
        d1 = np.arcsin(b1/a1*np.sin(t1))
        d11= np.arcsin(b2/a1*np.sin(t1))
        
        p1p2 = a1 / a2
        p1s1 = a1 / b1
        p1s2 = a1 / b2
        p2p1 = a2 / a1
        p2s1 = a2 / b1
        p2s2 = a2 / b2
        s2p1 = b2 / a1
        s2s1 = b2 / b1
        s2p2 = b2 / a2
        s1p1 = b1 / a1
        s1s2 = b1 / b2
        s1p2 = b1 / a2
        d1d2 = p1 / p2
        d2d1 = p2 / p1
        
        #Zoepritz equations
        A = np.array([[np.sin(t1), -np.cos(d1), -np.sin(t11), -np.cos(d11)],
                     [np.cos(t1), -np.sin(d1), np.cos(t11), np.sin(d11)],
                     [np.sin(2 * t1), p1s1 * np.cos(2 * d1), p1p2 * s2s1 * s2s1 * d2d1 * np.sin(2 * t11), -d2d1 * p1s1 * s2s1 * np.cos(2 * d11)],
                     [np.cos(2 * d1), -s1p1 * np.sin(2 * d1), -d2d1 * p2p1 * np.cos(2 * d11), -d2d1 * s2p1 * np.sin(2 * d11)]])
        
        B = np.array([[-np.sin(t1)],
                      [np.cos(t1)],
                      [np.sin(2 * t1)],
                      [-np.cos(2 * d1)]])


# =============================================================================
#         print(A)
#         print(B)
#         print(p1s1, np.cos(2 * d1), p1p2, s2s1, s2s1, d2d1, np.sin(2 * t11), -d2d1, p1s1, s2s1, np.cos(2 * d11))
#         print()
# =============================================================================
        
        #Solve the matrix equations
        c = np.dot(linalg.pinv(A), B)
        
        Rpp[i] = c[0]
    
    return Rpp
	
	
	
	
	
	