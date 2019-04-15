#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


#%%
"""
    Nome: spl3System

    Descrição: Monta o sistema linear para a determinação dos coeficientes ci da spline.

    Argumentos:

        x	-> Array com o valor da  variável independente dos pontos a serem interpolados
        y	-> Array com o valor da  variável dependente dos pontos a serem interpolados


    Retorno: A função deve retornar a matriz de coeficientes do sistema e o vetor com o lado direito do sistema.
"""

def spl3System( x , y ):

    assert ( x.shape == y.shape ) , "spl3System: x e y não possuem o mesmo número de elementos"

    n = len( x ) - 1
    
    M = np.zeros( [n+1,n+1] , dtype=float )
    r = np.zeros( n+1 , dtype=float )
    
    for i in range ( 1 , n ):
        r[i] = 3 * ( ( y[i+1] - y[i] ) / ( x[i+1] - x[i] ) - ( y[i] - y[i-1] ) / ( x[i] - x[i-1] ) )
        M[i,i+1] = x[i+1] - x[i]
        M[i,i-1] = x[i] - x[i-1]
        M[i,i] = 2 * ( M[i,i+1] + M[i,i-1] )
                
    M[0,0] = 1
    M[-1,-1] = 1

    return M , r


#%%
"""
    Nome: spl3Coeff

    Descrição: Calcula os coeficientes da spline.

    Argumentos:

        x	-> Array com o valor da  variável independente dos pontos a serem interpolados
        y	-> Array com o valor da  variável dependente dos pontos a serem interpolados

    Retorno: A função de retornar quatro arrays contendo, respectivamente, os coeficientes ai, bi, ci e di.
"""

def spl3Coeff( x , y ):
    
    assert( x.shape == y.shape ) , "spl3Coeff: x e y não possuem o mesmo número de elementos"

    n = len( x ) - 1

    a = np.zeros( [n] , dtype=float )
    b = np.zeros( [n] , dtype=float )
    c = np.zeros( [n] , dtype=float )
    d = np.zeros( [n] , dtype=float )

    M , r  = spl3System( x , y )

    vet = np.linalg.solve( M , r )

    for i in range( n ):
        h = ( x[i+1] - x[i] )
        a[i] = y[i]
        c[i] = vet[i]
        d[i] = ( vet[i+1] - vet[i] ) / ( 3 * h )
        b[i] = ( ( y[i+1] - y[i] ) / h ) - vet[i] * h - d[i] * h**2

    return a , b , c , d


#%%
"""
    Nome: spl3EvalS

    Descrição: Calcula o valor da spline cúbica em um dado ponto.

    Argumentos:

        a   -> Array com os coeficientes ai da spline cúbica
        b	-> Array com os coeficientes bi da spline cúbica
        c	-> Array com os coeficientes ci da spline cúbica
        d   -> Array com os coeficientes di da spline cúbica
        x   -> Array com o valor da  variável independente dos pontos a serem interpolados
        z   -> Ponto onde a interpolação será avaliada

    Retorno: O valor da spline no ponto z.
"""

def spl3EvalS( a , b , c , d , x , z ): #ESSA FOI PORRA!

    assert ( z >= x[0] and z <= x[-1] ) , "spl3EvalS: O ponto está fora do intervalo de interpolação"

    n = len( x ) - 1

    p = 0

    for i in range( n ):
    
        if( z>= x[i] and z<= x[i+1] ):
        
            p = i

    return a[p] + b[p] * ( z - x[p] ) + c[p] * ( z - x[p] )**2 + d[p] * ( z - x[p] )**3


#%%
"""
    Nome: spl3EvalD

    Descrição: Calcula o valor da derivada da spline cúbica em um dado ponto.

    Argumentos:

        a   -> Array com os coeficientes ai da spline cúbica
        b	-> Array com os coeficientes bi da spline cúbica
        c	-> Array com os coeficientes ci da spline cúbica
        d   -> Array com os coeficientes di da spline cúbica
        x   -> Array com o valor da  variável independente dos pontos a serem interpolados
        z   -> Ponto onde a interpolação será avaliada

    Retorno: O valor da derivada da spline no ponto z.
"""

def spl3EvalD( b , c , d , x , z ):
    
    assert ( z >= x[0] and z<= x[-1] ) , "spl3EvalD: O ponto está fora do intervalo de interpolação"

    n = len( x ) - 1

    p = 0

    for i in range( n ):
    
        if( z >= x[i] and z <= x[i+1] ):
        
            p = i

    return b[p] + 2 * c[p] * ( z - x[p] ) + 3 * d[p] * ( z - x[p] )**2


#%%
"""
    Nome: spl3Eval2D

    Descrição: Calcula o valor da segunda derivada da spline cúbica em um dado ponto.

    Argumentos:

        a   -> Array com os coeficientes ai da spline cúbica
        b	-> Array com os coeficientes bi da spline cúbica
        c	-> Array com os coeficientes ci da spline cúbica
        d   -> Array com os coeficientes di da spline cúbica
        x   -> Array com o valor da  variável independente dos pontos a serem interpolados
        z   -> Ponto onde a interpolação será avaliada

    Retorno: O valor da segunda derivada da spline no ponto z.
"""

def spl3Eval2D( c , d , x , z ):
    
    assert ( z >= x[0] and z<= x[-1] ) , "spl3Eval2D: O ponto está fora do intervalo de interpolação"

    n = len( x ) - 1

    p = 0

    for i in range( n ):
    
        if( z >= x[i] and z <= x[i+1] ):
        
            p = i

    return 2 * c[p] + 6 * d[p] * ( z - x[p] )


#%%
"""
    Nome: spl3Plot

    Descrição: Plota uma imagem das coordenadas interpoladas por splines cúbicas

    Argumentos:

        name    -> Nome do arquivo de coordenadas a ser plotado

    Retorno: -
"""

def spl3Plot( x , y , c='k' , ln=1 , ao='on' ):

    n = len( x )

    base = np.arange( n )

    ax , bx , cx , dx = spl3Coeff( base , x )
    ay , by , cy , dy = spl3Coeff( base , y )

    BASE = np.linspace( base[0] , base[-1] , 1000 )

    X = np.linspace( x[0] , x[-1] , 1000 )
    Y = np.linspace( y[0] , y[-1] , 1000 )

    for i in range( 1000 ):

        X[i] = spl3EvalS( ax , bx , cx , dx , base , BASE[i] )   
        Y[i] = spl3EvalS( ay , by , cy , dy , base , BASE[i] )   

    plt.plot( X , Y , color=c , linewidth=ln )
    plt.axis( ao )
    plt.show( )

    return


#%%
