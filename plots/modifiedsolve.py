from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import utils

l = 3
h = 1
b = 0.5
a = sqrt( ( ( l - b )**2 )/4  - h**2  )

Nvals = [321, 161, 81, 41, 21]

iter = 0
l2error = []
sumref = 0
uref = np.zeros( ( Nvals[0], Nvals[0] ) )
sumerror = []
linferror = []
deltax = []

def computeUpperBCForwardDiff( A, Nstride, Ny, N, a, b ):

    for i in range( 1, Nstride ):

        rowIdx = j*Nstride + i
        colIdx = rowIdx
        j = Ny - 1

        eta = utils.computePoint( j + 1, N )  
        zeta = utils.computePoint( i, N )

        k = ( 2*a*zeta )/( 2*a*eta + b )

        A[ rowIdx , colIdx ] = 3*k + 3

        if colIdx + 1 < tot_size:
            A[ rowIdx , colIdx + 1 ] = -4*k

        if colIdx + 2 < tot_size:
            A[ rowIdx , colIdx + 2 ] = k

        A[ rowIdx , colIdx - Nstride ] = -4
        A[ rowIdx , colIdx - 2*Nstride ] = 1


def computeUpperBCCentralDiff( A, Nstride, Ny, N, a, b ):

    for i in range( 1, Nstride ):

        rowIdx = j*Nstride + i
        colIdx = rowIdx
        j = Ny - 1

        eta = utils.computePoint( j + 1, N )  
        zeta = utils.computePoint( i, N )

        k = ( 2*a*zeta )/( 2*a*eta + b )

        A[ rowIdx , colIdx ] = 0

        if colIdx + 1 < tot_size:
            A[ rowIdx , colIdx + 1 ] = k

        

        A[ rowIdx , colIdx - Nstride ] = 1
        A[ rowIdx , colIdx + Nstride ] = -1

for N in Nvals:
    dzeta = 1/( N - 1 )

    tot_size = (N - 1)*( N - 1 )

    A = lil_matrix( (tot_size, tot_size) )
    f = np.zeros( ( tot_size, 1 ) )

    Nstride = N - 1
    Ny = N - 1

    X, Y = utils.computeXY( N, h, a, b )

    for j in range( Ny - 1 ):
        eta = utils.computePoint( j + 1, N )  

        for i in range( 1, Nstride ):
            zeta = utils.computePoint( i, N )

            rowIdx = j*Nstride + i
            colIdx = rowIdx
        
            c1 = 8*(a**2)*zeta
            
            c2 = 4*( a*a*zeta*zeta + h*h )

            c3 = -4*a*zeta*( 2*a*eta + b )

            c4 = ( 2*a*eta + b )**2

            A[ rowIdx , colIdx ] = ( 2*c2 + 2*c4 )

            if colIdx + 1 < tot_size:
                A[ rowIdx , colIdx + 1 ] = -c1*dzeta/( 2 ) -c2

            if colIdx - 1 >= 0:
                A[ rowIdx , colIdx - 1 ] = c1*dzeta/( 2) - c2

            A[ rowIdx , colIdx + Nstride ] = -c4

            if colIdx - Nstride >= 0:
                A[ rowIdx , colIdx - Nstride ] = -c4

            if colIdx + Nstride +  1 < tot_size:
                A[ rowIdx , colIdx + Nstride + 1 ] = -c3/( 4 )

            if colIdx - Nstride + 1 >= 0:
                A[ rowIdx , colIdx - Nstride + 1 ] = c3/( 4 )

            A[ rowIdx , colIdx + Nstride - 1 ] = c3/( 4 )

            if colIdx - Nstride - 1 >= 0:
                A[ rowIdx , colIdx - Nstride - 1 ] = -c3/( 4 )

            f[rowIdx] = dzeta*dzeta*( ( 2*a*eta + b )**2 )

    for j in range( Ny ):

        i = 0

        eta = utils.computePoint( j + 1, N )  
        zeta = utils.computePoint( i, N )

        rowIdx = j*Nstride
        colIdx = rowIdx

        A[ rowIdx , colIdx ] = 3
        A[ rowIdx , colIdx + 2 ] = 1
        A[ rowIdx , colIdx + 1 ] = -4

    for i in range( 1, Nstride ):

        rowIdx = j*Nstride + i
        colIdx = rowIdx
        j = Ny - 1

        eta = utils.computePoint( j + 1, N )  
        zeta = utils.computePoint( i, N )

        k = ( 2*a*zeta )/( 2*a*eta + b )

        A[ rowIdx , colIdx ] = 3*k + 3

        if colIdx + 1 < tot_size:
            A[ rowIdx , colIdx + 1 ] = -4*k

        if colIdx + 2 < tot_size:
            A[ rowIdx , colIdx + 2 ] = k

        A[ rowIdx , colIdx - Nstride ] = -4
        A[ rowIdx , colIdx - 2*Nstride ] = 1

    A = A.tocsr()
    u = spsolve( A, f )

    u = u.reshape( ( N - 1, N - 1 ) )

    u = np.append( u, np.zeros( (N - 1, 1) ), axis=1 )
    lowervals = np.zeros( ( 1, N ) )
    u = np.append( lowervals, u, axis=0 )

    sum = utils.computeSimpsons( u, N, a, h, b )
    print(sum)

    plt.figure()
    plt.pcolormesh( X, Y, u )
    plt.colorbar()
    plt.show()

    if iter > 0:

        currl2error = utils.computeL2Error( u, X, Y, N, uref, Xref, Yref, Nvals[0] )
        l2error.append( currl2error )

        currlinferror = utils.computeLinfError( u, X, Y, N, uref, Xref, Yref, Nvals[0] )
        linferror.append( currlinferror )

        currsumerror = abs( sum - sumref )
        sumerror.append( currsumerror )

        deltax.append( utils.computeSpacing( N ) )

    else:

        uref = u
        Xref = X
        Yref = Y
        sumref = sum

    iter += 1

l2error = np.array( l2error )
linferror = np.array( linferror )
sumerror = np.array( sumerror )
deltax = np.array( deltax )

plt.figure()
plt.yscale('log')
plt.xscale('log')
plt.plot( deltax, l2error, "bo-", label = "L2 Error Convergence" )
plt.plot( deltax, deltax**2, label = " ${\Delta x}^2$ " )
plt.legend()

plt.figure()
plt.yscale('log')
plt.xscale('log')
plt.plot( deltax, linferror/np.max(linferror), "bo-", label = "Linf Error Convergence" )
plt.plot( deltax, deltax**2,label = " ${\Delta x}^2$ ")
plt.legend()

print( deltax )
print( l2error )
print( linferror )
print( Nvals[0]*linferror )

plt.figure()
plt.yscale('log')
plt.xscale('log')
plt.plot( deltax, sumerror, "bo-", label = "Flowrate Error Convergence" )
plt.plot( deltax, deltax**2, label = " ${\Delta x}^2$ ")

plt.legend()
plt.show()