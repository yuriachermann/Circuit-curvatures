#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 21:41:03 2018

@author: yuri
"""

import numpy as np
import spl3
import curve


x , y = np.loadtxt( "Circuit of the Americas.txt" , unpack=True )

n = len( x )

tl = np.arange( n , dtype=float )

spl3.spl3Plot( x , y , 'k' , 4.5 , 'off' )

Rt = curve.curveRadius( x , y )

spl3.spl3Plot( tl , Rt , 'g' , 0.7 )
