from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

l = 3
h = 1
b = 0
a = sqrt( ( ( l - b )**2 )/4  - h**2  )

N = 11
dzeta = 1/( N - 1 )

tot_size = (N - 1)*( N - 1 )

A = np.zeros( (tot_size, tot_size) )
f = np.ones( ( tot_size, 1 ) )

def computePoint( idx, N ):

    return idx/( N - 1 )

def computeX( zeta, eta, a, b ):

    return zeta*( 2*a*eta + b )/( 2 )

def computeY( eta, h ):

    return eta*h

Nstride = N - 1
Ny = N - 1

X = np.zeros( ( N - 1, N - 1 ) )
Y = np.zeros( ( N - 1, N - 1 ) )

for j in range( Ny - 1 ):
    eta = computePoint( j + 1, N )  

    for i in range( 1, Nstride ):
        zeta = computePoint( i, N )

        Y[j][i] = computeY( eta, h )
        X[j][i] = computeX( zeta, eta, a, b )

        rowIdx = j*Nstride + i
        colIdx = rowIdx
      
        c1 = 8*(a**2)*zeta/( ( h*( 2*a*eta + b ) )**2 )
        
        c2 = 4*( a*a*zeta*zeta + h*h )/( ( h*( 2*a*eta + b ) )**2 )

        c3 = -4*a*zeta/( h*h*( 2*a*eta + b ) )

        c4 = 1/( h*h )

        A[ rowIdx ][ colIdx ] = ( 2*c2/( dzeta*dzeta ) + 2*c4/( dzeta*dzeta ) )

        if colIdx + 1 < tot_size:
            A[ rowIdx ][ colIdx + 1 ] = -c1/( 2*dzeta ) -c2/( dzeta*dzeta )

        if colIdx - 1 >= 0:
            A[ rowIdx ][ colIdx - 1 ] = c1/( 2*dzeta ) - c2/( dzeta*dzeta )

        A[ rowIdx ][ colIdx + Nstride ] = -c4/( dzeta*dzeta )

        if colIdx - Nstride >= 0:
            A[ rowIdx ][ colIdx - Nstride ] = -c4/( dzeta*dzeta )

        if colIdx + Nstride +  1 < tot_size:
            A[ rowIdx ][ colIdx + Nstride + 1 ] = -c3/( 4*dzeta*dzeta )

        if colIdx - Nstride + 1 >= 0:
            A[ rowIdx ][ colIdx - Nstride + 1 ] = c3/( 4*dzeta*dzeta )

        A[ rowIdx ][ colIdx + Nstride - 1 ] = c3/( 4*dzeta*dzeta )

        if colIdx - Nstride - 1 >= 0:
            A[ rowIdx ][ colIdx - Nstride - 1 ] = -c3/( 4*dzeta*dzeta )

for j in range( Ny ):

    i = 0

    eta = computePoint( j + 1, N )  
    zeta = computePoint( i, N )
    Y[j][i] = computeY( eta, h )
    X[j][i] = computeX( zeta, eta, a, b )

    rowIdx = j*Nstride
    colIdx = rowIdx
    f[rowIdx] = 0


    A[ rowIdx ][ colIdx ] = 3
    A[ rowIdx ][ colIdx + 2 ] = 1
    A[ rowIdx ][ colIdx + 1 ] = -4

for i in range( 1, Nstride ):

    rowIdx = j*Nstride + i
    colIdx = rowIdx
    j = Ny - 1
    f[rowIdx] = 0

    eta = computePoint( j + 1, N )  
    zeta = computePoint( i, N )

    Y[j][i] = computeY( eta, h )
    X[j][i] = computeX( zeta, eta, a, b )

    k = ( 2*a*zeta )/( 2*a*eta + b )

    A[ rowIdx ][ colIdx ] = 3*k + 3

    if colIdx + 1 < tot_size:
        A[ rowIdx ][ colIdx + 1 ] = -4*k

    if colIdx + 2 < tot_size:
        A[ rowIdx ][ colIdx + 2 ] = k

    A[ rowIdx ][ colIdx - Nstride ] = -4
    A[ rowIdx ][ colIdx - 2*Nstride ] = 1


# deleteRowIdx = np.arange( N - 1, N*N, N )
# deleteRowIdx = np.append( deleteRowIdx, np.arange( 0, N - 1 ) )

# A = np.delete( A, deleteRowIdx, 0 )
# A = np.delete( A, deleteRowIdx, 1 )

u = np.linalg.solve( A, f )

#u = spsolve( A, f )

u = u.reshape( ( N - 1, N - 1 ) )

plt.pcolormesh(X, Y, u)
plt.colorbar()
plt.show()