#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from quad import qgaussab
from spl3 import spl3Coeff , spl3EvalS , spl3EvalD , spl3Eval2D


#%%
"""
    Nome: curveLength

    Descrição:  Retorna um array com a distância do i-ésimo até o ponto inicial percorrendo a curva.

    Argumentos:

        x	-> Array com o valor da  variável independente dos pontos a serem interpolados
        y	-> Array com o valor da  variável dependente dos pontos a serem interpolados


    Retorno: Array com a distância de cada ponto ao ponto inicial percorrendo a curva.
"""

def curveLength( x , y ):
    
    assert( x.shape == y.shape ) , "curveLength: x e y não possuem o mesmo número de elementos"

    n = len( x )
    
    a , b , c , d = spl3Coeff( x , y )
    
    l = lambda z: spl3EvalD( a , b , c , d , x , z )
    
    f = lambda z: np.sqrt( 1 + l(z)**2 )
    
    s = np.zeros( [n] , dtype=float )
    
    for i in range( n ):

        for j in range( i ):

            s[i] += qgaussab( f , 3 , x[j] , x[j+1] )

    return s


#%%
"""
    Nome: curveParametric

    Descrição:  Monta as funções paramétricas x(t) e y(t) usando spline cúbica

    Argumentos:

        xa	-> Array com o valor da  variável independente dos pontos a serem interpolados
        ya	-> Array com o valor da  variável dependente dos pontos a serem interpolados

    Retorno: As funções x(t) e y(t)
"""

def curveParametric( xa , ya ):
    
    n = len( xa )
    
    tl = np.arange( n , dtype=float )

    def xt( t ):
    
        ax , bx , cx , dx = spl3Coeff( tl , xa )
        
        return spl3EvalS( ax , bx , cx , dx , tl , t )

    def yt( t ):
    
        ay , by , cy , dy = spl3Coeff( tl , ya )
        
        return spl3EvalS( ay , by , cy , dy , tl , t )
        
    return xt , yt


#%%
"""
    Nome: curveSParametric

    Descrição:  Monta as funções paramétricas x(s) e y(s) usando spline cúbica. 

    Argumentos:

        xa	-> Array com o valor da  variável independente dos pontos a serem interpolados
        ya	-> Array com o valor da  variável dependente dos pontos a serem interpolados

    Retorno: As funções x(s) e y(s)
"""

def curveSParametric( xa , ya ):
    
    sl = curveLength( xa , ya )
    
    a , b , c , d = spl3Coeff( xa , ya )

    def xs( s ):
    
        ax , bx , cx , dx = spl3Coeff( sl , xa )
        
        return spl3EvalS( ax , bx , cx , dx , sl , s )

    def ys( s ):
    
        ay , by , cy , dy = spl3Coeff( sl , ya )
        
        return spl3EvalS( ay , by , cy , dy , sl , s )
        
    return xs , ys


#%%
"""
    Nome: curveTangent

    Descrição:  Calcula a tangente em um dado ponto da curva.

    Argumentos:

        xa	-> Array com o valor da  variável independente dos pontos a serem interpolados
        ya	-> Array com o valor da  variável dependente dos pontos a serem interpolados
        s   -> Distância do primeiro ponto ao longo da curva para o qual desejemos conhecer o vetor tangente. 

    Retorno: As componentes x e y do vetor tangente no ponto especificado. Ou seja, o comando return será algo na forma 'return tx,ty' onde tx,ty são dois floats.
"""

def curveTangent( xa , ya , s ):
    
    sl = curveLength( xa , ya )
    
    a , b , c , d = spl3Coeff( xa , ya )

    def tx( u ):
    
        ax , bx , cx , dx = spl3Coeff( sl , xa )
        
        return spl3EvalD( ax , bx , cx , dx , sl , u )

    def ty( u ):
    
        ay , by , cy , dy = spl3Coeff( sl , ya )
        
        return spl3EvalD( ay , by , cy , dy , sl , u )
        
    return tx( s ) , ty( s )


#%%
"""
    Nome: curveRadius

    Descrição: Monta a função R(t) usando spline cúbica.

    Argumentos:

        x	-> Array com o valor da  variável independente dos pontos a serem interpolados
        y	-> Array com o valor da  variável dependente dos pontos a serem interpolados

    Retorno: Função R(t) que devolve o raio de curvatura na i-ésima posição do vetor t.
"""

def curveRadius( x , y ):
    
    assert( x.shape == y.shape ) , "curveRadius: x e y não possuem o mesmo número de elementos"

    n = len( x )

    tl = np.arange( n , dtype=float )
    
    def R( t ):        
                
        
        ax , bx , cx , dx = spl3Coeff( tl , x )
        
        ay , by , cy , dy = spl3Coeff( tl , y )
        
        DX = spl3EvalD( bx , cx , dx , tl , t )
        
        D2X = spl3Eval2D( cx , dx , tl , t )
        
        DY = spl3EvalD( by , cy , dy , tl , t )
        
        D2Y = spl3Eval2D( cy , dy , tl , t )
        
        return ( ( DX**2 + DY**2 )**(1.5) ) / abs( DX * D2Y - DY * D2X )
    
    Rt = np.zeros( n , dtype=float )
    
    for i in range( n ):
        
        Rt[i] = R( i )
        
    Rt[0] = Rt[1]
    Rt[-1] = Rt[-2]
    
    return Rt


#%%
