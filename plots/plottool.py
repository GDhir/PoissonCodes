from email.mime import base
import numpy as np
import matplotlib.pyplot as plt
import os
from math import sqrt

import utils

def computeGrid( N, l, h, b, a ):

    X = np.zeros( (N, N) )
    Y = np.zeros( (N, N) )

    for j in range(N):

        eta = utils.computePoint( j + 1, N )  

        for i in range(N):

            zeta = utils.computePoint( i, N )

            Y[j, i] = utils.computeY( eta, h )
            X[j, i] = utils.computeX( zeta, eta, a, b )

    return [X, Y]

def computeErrorCpp( prefixdir, a, h, b ):

    Nvals = [641, 321, 161, 81, 41]
    l2error = []
    uref = np.zeros( ( Nvals[0], Nvals[0] ) )
    linferror = []
    deltax = []
    sumref = 0
    sum = 0
    sumerror = []

    iter = 0
    for N in Nvals:
        xdata, ydata, udata = utils.getSolutionData( prefixdir + "solnsN=" + str(N) + ".txt", N )
        udata = udata.reshape( ( N, N ))
        xdata = xdata.reshape( (N, N) )
        ydata = ydata.reshape( (N, N) )
        
        sum = utils.computeSimpsons( udata, N, a, h, b )
        filename = prefixdir + "N=" + str(N)
        utils.plotSolution( filename, xdata, ydata, udata )


        if iter > 0:

            currl2error = utils.computeL2Error( udata, xdata, ydata, N, uref, Xref, Yref, Nvals[0] )
            l2error.append( currl2error )

            currlinferror = utils.computeLinfError( udata, xdata, ydata, N, uref, Xref, Yref, Nvals[0] )
            linferror.append( currlinferror )

            currsumerror = abs( sum - sumref )
            sumerror.append( currsumerror )

            deltax.append( utils.computeSpacing( N ) )

        else:

            uref = udata
            Xref = xdata
            Yref = ydata
            sumref = sum

        
        iter += 1

    return [deltax, l2error, linferror, sumerror]

def plotSol(prefixdir, a, h, b):

    Nvals = [161]

    for N in Nvals:
        xdata, ydata, udata = utils.getSolutionData( prefixdir + "solnsN=" + str(N) + ".txt", N )
        udata = udata.reshape( ( N, N ))
        xdata = xdata.reshape( (N, N) )
        ydata = ydata.reshape( (N, N) )
        
        sum = utils.computeSimpsons( udata, N, a, h, b )
        filename = prefixdir + "N=" + str(N)
        utils.plotSolution( filename, xdata, ydata, udata )

def computePareto( ):

    from scipy.spatial import ConvexHull, convex_hull_plot_2d

    l = 3
    hvals = np.arange( 0.1, 1.1, 0.1 )
    bvals = np.arange( 0, 1.1, 0.1 )
    N = 41
    prefixstr = "/home/gaurav/CS6220/HW1/plots/pareto/"

    QIinfo = []

    filename = "qvalpareto.txt"

    [hvals, bvals, sumvals] = utils.getSolutionData( prefixstr + filename, N )

    for idx in range( len(hvals) ):
        h = hvals[idx]
        b = bvals[idx]
        I = utils.computeMomentofInertia( h, b, l )
        QIinfo.append( [ I, sumvals[idx] ] )

    QIinfo = np.array( QIinfo )
    hull = ConvexHull( QIinfo )

    plt.plot( QIinfo[ hull.vertices, 0 ], QIinfo[ hull.vertices, 1 ], 'r--', lw=2, label="Convex Hull" )
    plt.plot( QIinfo[ :, 0 ], QIinfo[ :, 1 ], 'bo', label="($I, Q$)")
    plt.xlabel("Moment of Inertia")
    plt.ylabel("Flowrate")
    plt.legend()
    plt.savefig( prefixstr + "paretoplot.png" )
    plt.show()

if __name__ == "__main__":

    # os.system("/usr/bin/cmake --build /home/gaurav/CS6220/HW1/build --config Debug --target all -- -j 10")
    # os.system( "/home/gaurav/CS6220/HW1/build/bin/solver" )

    # bvals = [1.0]
    # b = 0.7
    # l = 3
    # h = 1
    # #a = sqrt( ( ( l - b )**2 )/4  - h**2  )
    # prefix = "/home/gaurav/CS6220/HW1/plots/pareto" + str(b) + "/"
    # plotSol( prefix, a, h, b )

    computePareto()

    # for b in bvals:

    #     a = sqrt( ( ( l - b )**2 )/4  - h**2  )
    #     prefix = "/home/gaurav/CS6220/HW1/plots/plotdir_b=" + str(b) + "/"
    #     deltax, l2Error, linfError, qvalError = computeErrorCpp( prefix, a, h, b)

    #     l2Error = np.array( l2Error )
    #     linfError = np.array( linfError )
    #     deltax = np.array( deltax )

    #     plt.figure()
    #     plt.yscale('log')
    #     plt.xscale('log')
    #     plt.plot( deltax, l2Error, "bo-", label = "L2 Error Convergence" )
    #     plt.plot( deltax, 2*deltax**2, label = " ${\Delta x}^2$ " )
    #     titledata = "L2 Error Convergence plot for b=" + str(b)
    #     #plt.title( label=titledata )
    #     plt.legend()
    #     plt.xlabel( "$\Delta x$" )
    #     plt.ylabel( "$L_2$ Error" )
    #     filename = prefix + "l2error"
    #     plt.savefig( filename + "1.png" )

    #     plt.figure()
    #     plt.yscale('log')
    #     plt.xscale('log')
    #     plt.plot( deltax, linfError, "bo-", label = "Linf Error Convergence" )
    #     plt.plot( deltax, 2*deltax**2, label = " ${\Delta x}^2$ " )
    #     titledata = "Linf Error Convergence plot for b=" + str(b)
    #     plt.xlabel( "$\Delta x$" )
    #     plt.ylabel( "$L_{\infty}$ Error" )
    #     #plt.title( label=titledata )
    #     plt.legend()
    #     filename = prefix + "linferror"
    #     plt.savefig( filename + "1.png" )

    #     print( deltax )
    #     print( l2Error )
    #     print( linfError )
    #     print( linfError*321 )
    #     #plt.plot( dx )

    #     plt.figure()
    #     plt.yscale('log')
    #     plt.xscale('log')
    #     plt.plot( deltax, qvalError, "bo-", label = "Flowrate Error Convergence" )
    #     plt.plot( deltax, 2*deltax**2 )
    #     titledata = "Flowrate Error Convergence plot for b=" + str(b)
    #     #plt.title( label=titledata )
    #     plt.xlabel( "$\Delta x$" )
    #     plt.ylabel( "Flowrate Error" )
    #     filename = prefix + "flowrateerror"
    #     plt.savefig( filename + "1.png" )
    #     plt.show()
