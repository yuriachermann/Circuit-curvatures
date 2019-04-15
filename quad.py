#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


#%%
"""
    Nome: qgauss

    Descrição: Calcula aproximadamente a integral de uma função no intervalo [−1,+1] utilizando a quadratura de Gauss-Legendre.

    Argumentos:

        f	-> A função a ser integrada 
        n	-> Número de nós utilizado na quadratura

    Retorno: O valor da integral
"""

def qgauss( f , n ):
    
    x , c = np.polynomial.legendre.leggauss( n )

    sum=0

    for i in range ( n ):

        sum += c[i] * f( x[i] )

    return sum


#%%
"""
    Nome: qgaussab

    Descrição: Calcula aproximadamente a integral de uma função no intervalo [a,b] utilizando a quadratura de Gauss-Legendre.

    Argumentos:

        f	-> A função a ser integrada 
        n	-> Número de nós utilizado na quadratura
        a   -> Extremo inferior do intervalo de integração
        b   -> Extremo superior do intervalo de integração

    Retorno: O valor da integral
"""

def qgaussab( f , n , a , b ):
    
    g = lambda u: f( a + ( u + 1 ) * ( b - a ) * 0.5 )

    x , c = np.polynomial.legendre.leggauss( n )

    sum=0

    for i in range ( n ):

        sum += c[i] * g( x[i] )

    return ( b - a ) * 0.5 * sum


#%%