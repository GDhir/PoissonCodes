
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
import utils

def computeSolution( N, h, a, b, tol ):

    dzeta = 1/( N - 1 )

    u = np.zeros( ( N, N ) )

    error = 10**5

    prevval = 0

    iter = 0
    while( iter < 100 ):
        iter += 1

        error = 0
        for j in range( 1, N - 1 ):

            eta = utils.computePoint( j, N )

            prevval = u[ j, 0 ]
            u[ j, 0 ] = ( 18*u[ j, 1 ] -9*u[ j, 2 ] + 2*u[ j, 3 ] )/11
            error += ( u[ j, 0 ] - prevval )**2

            for i in range( 1, N - 1 ):

                zeta = utils.computePoint( i, N )

                c1 = 8*(a**2)*zeta
                
                c2 = 4*( a*a*zeta*zeta + h*h )

                c3 = -4*a*zeta*( 2*a*eta + b ) 

                c4 = ( 2*a*eta + b )**2

                c5 = ( h*( 2*a*eta + b ) )**2

                prevval = u[ j, i ]
                u[ j, i ] = ( c1*dzeta*( u[ j, i + 1 ] - u[ j, i - 1 ] )/2 + c2*( u[ j, i + 1 ] + u[ j, i - 1 ] ) + \
                c3*( u[ j + 1, i + 1 ] - u[ j - 1, i + 1 ] - u[ j + 1, i - 1 ] + u[ j - 1, i - 1 ] )/4 + \
                c4*( u[ j + 1, i ] + u[ j - 1, i ] ) + c5*dzeta*dzeta )/( 2*c2 + 2*c4 )

                error += ( u[ j, i ] - prevval )**2

        for i in range( N - 1 ):

            j = N - 1
            eta = utils.computePoint( j, N )
            zeta = utils.computePoint( i, N )
            k = utils.computek( eta, zeta, a, b )

            prevval = u[ j, i ]

            if i + 3 < N:
                u[ j, i ] = ( 18*u[ j - 1, i ] - 9*u[ j - 2, i ] + 2*u[ j - 3, i ] + 18*k*u[ j, i + 1 ] \
                - 9*k*u[ j, i + 2 ] + 2*k*u[ j, i + 3 ] )/( 11*k + 11 )
            elif i + 2 < N:
                u[ j, i ] = ( 18*u[ j - 1, i ] - 9*u[ j - 2, i ] + 2*u[ j - 3, i ] + 18*k*u[ j, i + 1 ] \
                - 9*k*u[ j, i + 2 ] )/( 11*k + 11 )
            else:
                u[ j, i ] = ( 18*u[ j - 1, i ] - 9*u[ j - 2, i ] + 2*u[ j - 3, i ] + 18*k*u[ j, i + 1 ] )/( 11*k + 11 )

            error += ( u[ j, i ] - prevval )**2

        error = sqrt(error)
        print( error )

    return u        


prefix = "/home/gaurav/plots/PythonGaussSiedel_Solution/"

l = 3
h = 1
b = 0.5
a = utils.computea( h, b, l)
tol = 10**(-4)
Nvals = [321, 161, 81, 41, 21]
sumvals = []

iter = 0
l2error = []
sumref = 0
uref = np.zeros( ( Nvals[0], Nvals[0] ) )
sumerror = []
linferror = []
deltax = []

prefix += "plotsb=" + str(b) + "/"

for N in Nvals:

    u = computeSolution( N, h, a, b, tol )
    X, Y = utils.computeXY( N, h, a, b )

    sum = utils.computeSimpsons( u, N, a, h, b )
    print(sum)
    sumvals.append(sum)

    filename = prefix + "udataN=" + str(N)
    utils.saveSolutionData( filename, u )


    plt.figure()
    plt.pcolormesh( X, Y, u )
    plt.colorbar()
    titledata = "Solution plot for b=" + str(b) + " and N=" + str(N)
    plt.title( label=titledata )
    filename = prefix + "umeshN=" + str(N)
    plt.savefig( filename + ".png" )
    plt.show()

    plt.figure()
    plt.contour( X,Y,u )
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
titledata = "L2 Error Convergence plot for b=" + str(b)
plt.title( label=titledata )
filename = prefix + "l2error"
plt.savefig( filename + ".png" )
# plt.close()
utils.savedata( filename + ".txt", l2error )

plt.figure()
plt.yscale('log')
plt.xscale('log')
plt.plot( deltax, linferror/np.max(linferror), "bo-", label = "Linf Error Convergence" )
plt.plot( deltax, deltax**2,label = " ${\Delta x}^2$ ")
plt.legend()
titledata = "Linf Error Convergence plot for b=" + str(b)
plt.title( label=titledata )
filename = prefix + "linferror"
plt.savefig( filename + ".png" )
# plt.close()
utils.savedata( filename + ".txt" , linferror)


print( deltax )
print( l2error )
print( linferror )
print( Nvals[0]*linferror )

plt.figure()
plt.yscale('log')
plt.xscale('log')
plt.plot( deltax, sumerror, "bo-", label = "Flowrate Error Convergence" )
plt.plot( deltax, deltax**2, label = " ${\Delta x}^2$ ")
titledata = "Flowrate Error Convergence plot for b=" + str(b)
plt.title( label=titledata )
plt.legend()
filename = prefix + "flowrateerror"
plt.savefig( filename + ".png" )
# plt.close()
utils.savedata( filename + ".txt", sumerror )

filename = prefix + "flowrate"
utils.savedata( filename + ".txt", sumvals )

plt.show()
plt.pause(3)
#plt.close("all")










