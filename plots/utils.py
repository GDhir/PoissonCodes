from math import sqrt, isnan
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from scipy import interpolate

def computePoint( idx, N ):

    return idx/( N - 1 )

def computeX( zeta, eta, a, b ):

    return zeta*( 2*a*eta + b )/( 2 )

def computeY( eta, h ):

    return eta*h

def computeSpacing( N ):

    return 1/(N - 1)

def computeSimpsons( u, N, a, h, b ):

    f = np.array( u )

    for j in range( N ):

        eta = computePoint( j, N )

        for i in range( N ):

            J = (a*eta + b/2)*h

            f[j, i] = f[j, i]*J

    dzeta = computeSpacing( N )

    sumvals = np.zeros( N )

    # for j in range(N):

    #     sumvals[j] += f[j, 0] + f[j, N - 1 ]

    #     for i in range( 1, N - 1 ):

    #         if( i%2 != 0 ):

    #             sumvals[j] += 4*f[j, i]

    #         else:

    #             sumvals[j] += 2*f[j, i]

    for j in range(N):

        sumvals[j] += 2*f[j, N - 1 ] + 2*f[j, 0]

        for i in range( 1, N - 1 ):

            if( i%2 != 0 ):

                sumvals[j] += 8*f[j, i]

            else:

                sumvals[j] += 4*f[j, i]

    
    sum = sumvals[0] + sumvals[-1]

    for j in range( 1, N - 1 ):

        if j%2 != 0:
            sum += 4*sumvals[j]

        else:
            sum += 2*sumvals[j]

    return 2*sum*dzeta*dzeta/9

def computeXY( N, h, a, b ):

    X = np.zeros( ( N, N ) )
    Y = np.zeros( ( N, N ) )

    for j in range(N):
        eta = computePoint( j, N )  

        for i in range(N):
            zeta = computePoint( i, N )

            Y[j][i] = computeY( eta, h )
            X[j][i] = computeX( zeta, eta, a, b )

    return [ X, Y ]

def computeL2Error( u, X, Y, N, uref, Xref, Yref, Nref ):

    points = []
    values = []
    for j in range(N):
        for i in range(N):

            points.append( [ X[j, i], Y[j, i] ] )
            values.append( u[j, i] )

    f = interpolate.griddata(points, values, (Xref, Yref), method='cubic')

    error = 0
    for j in range(  Nref):

        for i in range( Nref ):

            if isnan( f[j, i] ):
                m = 1
                #print( f[j, i], j, i, sep=" " )
            else:
                error += ( f[j, i] - uref[j, i] )**2

    return sqrt( error )

def computeLinfError( u, X, Y, N, uref, Xref, Yref, Nref ):

    points = []
    values = []
    for j in range(N):
        for i in range(N):

            points.append( [ X[j, i], Y[j, i] ] )
            values.append( u[j, i] )

    f = interpolate.griddata(points, values, (Xref, Yref), method='cubic')

    error = 0
    for j in range(  Nref ):
        for i in range( Nref - 1 ):

            # if isnan( f[j][i] ):
            #     print( f[j][i], j, i, sep=" " )
            # else:
            error = max( error, abs( f[j, i] - uref[j, i] ) )

    return error

def computek( eta, zeta, a, b ):

    return ( 2*a*zeta )/( 2*a*eta + b )

def computeUpperBCForwardDiff( A, Nstride, Ny, N, a, b, tot_size ):

    for i in range( 1, Nstride ):

        j = Ny - 1
        rowIdx = j*Nstride + i
        colIdx = rowIdx

        eta = computePoint( j + 1, N )  
        zeta = computePoint( i, N )

        k = computek( eta, zeta, a, b )

        A[ rowIdx , colIdx ] = 3*k + 3

        if colIdx + 1 < tot_size:
            A[ rowIdx , colIdx + 1 ] = -4*k

        if colIdx + 2 < tot_size:
            A[ rowIdx , colIdx + 2 ] = k

        A[ rowIdx , colIdx - Nstride ] = -4
        A[ rowIdx , colIdx - 2*Nstride ] = 1

    return A


def computeUpperBCCentralDiff( A, Nstride, Ny, N, a, b, tot_size ):

    for i in range( 1, Nstride ):

        rowIdx = j*Nstride + i
        colIdx = rowIdx
        j = Ny - 1

        eta = computePoint( j + 1, N )  
        zeta = computePoint( i, N )

        k = ( 2*a*zeta )/( 2*a*eta + b )

        A[ rowIdx , colIdx ] = 0

        if colIdx + 1 < tot_size:
            A[ rowIdx , colIdx + 1 ] = k

        

        A[ rowIdx , colIdx - Nstride ] = 1
        A[ rowIdx , colIdx + Nstride ] = -1


def computeUpperBCThirdOrder( A, Nstride, Ny, N, a, b, tot_size ):

    for i in range( 1, Nstride ):

        j = Ny - 1
        rowIdx = j*Nstride + i
        colIdx = rowIdx

        eta = computePoint( j + 1, N )  
        zeta = computePoint( i, N )

        k = ( 2*a*zeta )/( 2*a*eta + b )

        A[ rowIdx , colIdx ] = 11*k + 11

        if colIdx + 1 < tot_size:
            A[ rowIdx , colIdx + 1 ] = -18*k

        if colIdx + 2 < tot_size:
            A[ rowIdx , colIdx + 2 ] = 9*k

        if colIdx + 3 < tot_size:
            A[ rowIdx , colIdx + 3 ] = -2*k

        A[ rowIdx , colIdx - Nstride ] = -18
        A[ rowIdx , colIdx - 2*Nstride ] = 9
        A[ rowIdx , colIdx - 3*Nstride ] = -2



    return A

def computeLeftBoundarySecondOrder( A, Ny, Nstride  ):

    for j in range( Ny ):

        i = 0

        rowIdx = j*Nstride
        colIdx = rowIdx

        A[ rowIdx , colIdx ] = 3
        A[ rowIdx , colIdx + 2 ] = 1
        A[ rowIdx , colIdx + 1 ] = -4

    return A


def computeLeftBoundaryFourthOrder( A, Ny, Nstride ):

    for j in range( Ny ):

        i = 0

        rowIdx = j*Nstride + i
        colIdx = rowIdx

        A[ rowIdx , colIdx ] = -11
        A[ rowIdx , colIdx + 1 ] = 18
        A[ rowIdx , colIdx + 2 ] = -9
        A[ rowIdx, colIdx + 3 ] = 2

    return A


def computea( h, b, l ):
    return sqrt( ( ( l - b )**2 )/4  - h**2  )


def computeMomentofInertia( h, b, l ):

    d = ( l - b )/2
    y = h*d/( 2*d + b )
    t = 0.05

    I = (y**2)*b*t + 2*d*t*( h**2 - 3*h*y + 3*y**2 )/3

    return I

def saveSolutionData( filename, u ):
    
    np.savetxt( filename, u )

def savedata( filename, values ):

    np.savetxt( filename, values, delimiter="," )

def computeInnerSecondOrder( N, Nstride, Ny, tot_size, a, h, b, dzeta ):

    A = lil_matrix( (tot_size, tot_size) )
    f = np.zeros( ( tot_size, 1 ) )

    for j in range( Ny - 1 ):
        eta = computePoint( j + 1, N )  

        for i in range( 1, Nstride ):
            zeta = computePoint( i, N )

            rowIdx = j*Nstride + i
            colIdx = rowIdx
        
            c1 = 8*(a**2)*zeta/( ( h*( 2*a*eta + b ) )**2 )
            
            c2 = 4*( a*a*zeta*zeta + h*h )/( ( h*( 2*a*eta + b ) )**2 )

            c3 = -4*a*zeta/( h*h*( 2*a*eta + b ) )

            c4 = 1/( h*h )

            A[ rowIdx , colIdx ] = ( 2*c2 + 2*c4 )/(dzeta*dzeta)

            if colIdx + 1 < tot_size:
                A[ rowIdx , colIdx + 1 ] = -c1*dzeta/( 2*(dzeta*dzeta) ) -c2/(dzeta*dzeta)

            if colIdx - 1 >= 0:
                A[ rowIdx , colIdx - 1 ] = c1*dzeta/( 2*(dzeta*dzeta)  ) - c2/(dzeta*dzeta)

            A[ rowIdx , colIdx + Nstride ] = -c4/(dzeta*dzeta)

            if colIdx - Nstride >= 0:
                A[ rowIdx , colIdx - Nstride ] = -c4/(dzeta*dzeta)

            if colIdx + Nstride +  1 < tot_size:
                A[ rowIdx , colIdx + Nstride + 1 ] = -c3/( 4*(dzeta*dzeta) )

            if colIdx - Nstride + 1 >= 0:
                A[ rowIdx , colIdx - Nstride + 1 ] = c3/( 4*(dzeta*dzeta) )

            A[ rowIdx , colIdx + Nstride - 1 ] = c3/( 4*(dzeta*dzeta) )

            if colIdx - Nstride - 1 >= 0:
                A[ rowIdx , colIdx - Nstride - 1 ] = -c3/( 4*(dzeta*dzeta) )

            f[rowIdx] = 1

    return [A, f]


def getData( filename ):
    with open( filename ) as f:
        lines = f.readlines()

    return lines

def convertToFloat( arr ):

    return [ float(elem) for elem in arr ]

def getSolutionData( filename, N ):

    sol = getData( filename )
    
    xdata = []
    ydata = []
    udata = []

    for line in sol:

        values = line.split( "\t" )
        values = [ float(value) for value in values ]

        xdata.append( values[0] )
        ydata.append( values[1] )

        if len(values) == 3:
            udata.append( values[2] )

    udata = np.array( udata )
    xdata = np.array( xdata )
    ydata = np.array( ydata )

    return [ xdata, ydata, udata ]
            

def plotSolution( filename, xdata, ydata, udata ):

    plt.figure()
    plt.pcolormesh( xdata, ydata, udata )
    plt.colorbar()
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.savefig( filename +".png" )

    plt.figure()
    a = plt.contour( xdata, ydata, udata )
    plt.colorbar(a)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.savefig( filename + "contour.png" )

    plt.pause(3)
    plt.close("all")