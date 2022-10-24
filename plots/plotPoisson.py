from sys import prefix
import numpy as np
import matplotlib.pyplot as plt
import os
from math import sqrt

import utils    

def getBestW():

    prefixstr = "/home/gaurav/plots/PoissonMG/GS/"
    N = 25

    totalIters = []
    wvals = np.arange( 0.10, 2, 0.01 )
    bestiter = 10**6
    bestw = 0

    for w in wvals:

        residual, niter, uval = utils.getSolutionData( prefixstr + "residualN=" + str(N) + "w=" + str( round(w, 2) ) + ".txt", N )

        totalIters.append( niter[-1] )

        if niter[-1] < bestiter:
            bestiter = niter[-1]
            bestw = w


        # plt.figure()
        # plt.yscale('log')
        # plt.xscale('log')
        # plt.xlabel( "Number of Iterations" )
        # plt.ylabel( "Residual" )
        # plt.plot( niter, residual, "-" )
        # plt.show()

    print( bestw )
    print( bestiter )
    plt.figure()
    plt.plot( wvals, totalIters, "-" )
    plt.yscale("log")
    plt.xlabel( "$\omega$" )
    plt.ylabel( "Number of Iterations" )
    plt.savefig( prefixstr + "NIterVSw.png" )
    plt.show()

def computeMG():
    prefixstr = "/home/gaurav/plots/PoissonMG/MG/"
    N = 25
    w = 0.50

    os.system("/usr/bin/cmake --build /home/gaurav/build --config Debug --target all -- -j 10")
    os.system( "/home/gaurav/CS6220/HW1/build/bin/solver" )

    xdata, ydata, udata = utils.getSolutionData( prefixstr + "solnsN=" + str(N) + "w=" + "{:.2f}".format(w) + ".txt", N )

    residual, niter, uval = utils.getSolutionData( prefixstr + "residualN=" + str(N) + "w=" + "{:.2f}".format(w) + ".txt", N )

    udata = udata.reshape( ( N, N ))
    xdata = xdata.reshape( (N, N) )
    ydata = ydata.reshape( (N, N) )

    utils.plotSolution( prefixstr + "N=" + str(N) + "w=" + str( round(w, 2) ), xdata, ydata, udata )

    plt.figure()
    plt.yscale('log')
    #plt.xscale('log')
    plt.xlabel( "Number of Iterations" )
    plt.ylabel( "Residual" )
    plt.plot( niter, residual, "-" )
    plt.xlabel( "Number of Gauss Siedel Iterations" )
    plt.ylabel( "Residual" )
    plt.savefig( prefixstr + "NIterVSResidual.png" )
    plt.show()

def computeIterative( method ):

    prefixstr = "/home/gaurav/plots/PoissonMG/" + method + "/"
    N = 25
    w = 0.8

    # os.system("/usr/bin/cmake --build /home/gaurav/CS6220/HW1/build --config Debug --target all -- -j 10")
    # os.system( "/home/gaurav/CS6220/HW1/build/bin/solver" )

    xdata, ydata, udata = utils.getSolutionData( prefixstr + "solnsN=" + str(N) + "w=" + "{:.2f}".format(w) + ".txt", N )

    residual, niter, uval = utils.getSolutionData( prefixstr + "residualN=" + str(N) + "w=" + "{:.2f}".format(w) + ".txt", N )

    udata = udata.reshape( ( N, N ))
    xdata = xdata.reshape( (N, N) )
    ydata = ydata.reshape( (N, N) )

    utils.plotSolution( prefixstr + "N=" + str(N) + "w=" + str( round(w, 2) ), xdata, ydata, udata )

    plt.figure()
    # plt.yscale('log')
    #plt.xscale('log')
    plt.xlabel( "Number of Iterations" )
    plt.ylabel( "Residual" )
    plt.plot( niter, residual, "-" )
    plt.savefig( prefixstr + "NIterVSResidual.png" )
    plt.show()

def compareConvergence():

    prefixstr = "/home/gaurav/plots/PoissonMG/"
    N = 241
    wvals = [0.8]
    # dirvals = [ "GS", "Jacobi", "MG" ]
    # labels = [ "Gauss Siedel", "Jacobi", "Multigrid" ]
    dirvals = [ "MG" ]
    labels = [ "Gauss Siedel"]

    # os.system("/usr/bin/cmake --build /home/gaurav/CS6220/HW1/build --config Debug --target all -- -j 10")
    # os.system( "/home/gaurav/CS6220/HW1/build/bin/solver" )

    # xdata, ydata, udata = utils.getSolutionData( prefixstr + "solnsN=" + str(N) + "w=" + "{:.2f}".format(w) + ".txt", N )

    for w in wvals:

        plt.figure()
        idx = 0
        for dirval in dirvals:
            residual, niter, uval = utils.getSolutionData( prefixstr + dirval + "/residualN=" + str(N) + "w=" + "{:.2f}".format(w) + ".txt", N )

            # udata = udata.reshape( ( N, N ))
            # xdata = xdata.reshape( (N, N) )
            # ydata = ydata.reshape( (N, N) )

            # utils.plotSolution( prefixstr + "N=" + str(N) + "w=" + str( round(w, 2) ), xdata, ydata, udata )

            # plt.yscale('log')
            #plt.xscale('log')
            plt.xlabel( "Number of Iterations" )
            plt.ylabel( "Residual" )
            plt.plot( niter, residual, "-o" , label = labels[idx])
            #plt.show()
            idx += 1

        plt.legend()
        plt.savefig( prefixstr + "N=" + str(N) + "w=" + "{:.2f}".format(w) + "NIterVSResidualIterative.png" )
        plt.show()

def compareMGConvergence():

    prefixstr = "/home/gaurav/plots/PoissonMG/MG/"
    N = 25
    wvals = [0.50, 0.8, 1.0]

    # os.system("/usr/bin/cmake --build /home/gaurav/CS6220/HW1/build --config Debug --target all -- -j 10")
    # os.system( "/home/gaurav/CS6220/HW1/build/bin/solver" )

    # xdata, ydata, udata = utils.getSolutionData( prefixstr + "solnsN=" + str(N) + "w=" + "{:.2f}".format(w) + ".txt", N )

    for w in wvals:

        residual, niter, uval = utils.getSolutionData( prefixstr + "residualN=" + str(N) + "w=" + "{:.2f}".format(w) + ".txt", N )

        # udata = udata.reshape( ( N, N ))
        # xdata = xdata.reshape( (N, N) )
        # ydata = ydata.reshape( (N, N) )

        # utils.plotSolution( prefixstr + "N=" + str(N) + "w=" + str( round(w, 2) ), xdata, ydata, udata )

        # plt.yscale('log')
        #plt.xscale('log')
        plt.xlabel( "Number of Iterations" )
        plt.ylabel( "Residual" )
        plt.plot( niter, residual, "-o" , label = "w = {:.2f}".format(w) )
        #plt.show()

        plt.legend()

    plt.savefig( prefixstr + "N=" + str(N) + "NIterVSResidualForW.png" )

def compareMGGrids():

    prefixstr = "/home/gaurav/plots/PoissonMG/MG/"
    Nvals = [ 49, 97 ]
    wvals = [0.50, 0.8, 1.0]
    labels = [ "(1, 2, 1)", "(2, 4, 2)" ]

    for w in wvals:

        idx = 0
        plt.figure()
        for N in Nvals:
            residual, niter, uval = utils.getSolutionData( prefixstr + "residualN=" + str(N) + "w=" + "{:.2f}".format(w) + ".txt", N )

            xdata, ydata, udata = utils.getSolutionData( prefixstr + "solnsN=" + str(N) + "w=" + "{:.2f}".format(w) + ".txt", N )
            udata = udata.reshape( ( N, N ))
            xdata = xdata.reshape( (N, N) )
            ydata = ydata.reshape( (N, N) )

            #utils.plotSolution( prefixstr + "N=" + str(N) + "w=" + str( round(w, 2) ), xdata, ydata, udata )

            plt.yscale('log')
            #plt.xscale('log')
            plt.xlabel( "Number of Iterations" )
            plt.ylabel( "Residual" )
            plt.plot( niter, residual, "-o" , label = "MG at" + labels[idx] )
            #plt.show()

            plt.legend()
            idx += 1

        plt.savefig( prefixstr + "w = {:.2f}".format(w) + "NIterVSResidualForMGGrids.png" )

def compareMGLevels():

    prefixstr = "/home/gaurav/plots/PoissonMG/MG/"
    N = 97
    nlevels = [ 2, 3, 4 ] 
    wvals = [0.50, 0.8, 1.0]

    for w in wvals:

        idx = 0
        plt.figure()
        for maxlevel in nlevels:
            residual, niter, uval = utils.getSolutionData( prefixstr + "residualN=" + str(N) + "w=" + "{:.2f}".format(w) + "maxlevels=" + str(maxlevel) + ".txt", N )

            xdata, ydata, udata = utils.getSolutionData( prefixstr + "solnsN=" + str(N) + "w=" + "{:.2f}".format(w) + "maxlevels=" + str(maxlevel) + ".txt", N )
            udata = udata.reshape( ( N, N ))
            xdata = xdata.reshape( (N, N) )
            ydata = ydata.reshape( (N, N) )

            #utils.plotSolution( prefixstr + "N=" + str(N) + "w=" + str( round(w, 2) ), xdata, ydata, udata )

            plt.yscale('log')
            #plt.xscale('log')
            plt.xlabel( "Number of Iterations" )
            plt.ylabel( "Residual" )
            plt.plot( niter, residual, "-o" , label = "cycles = " + str(maxlevel) )
            #plt.show()

            plt.legend()
            idx += 1

        plt.savefig( prefixstr + "w = {:.2f}".format(w) + "NIterVSResidualForMGLevels.png" )

    Nvals = [49, 97]

    for w in wvals:

        idx = 0

        plt.figure()

        for N in Nvals:
            totaliter = []
            for maxlevel in nlevels:
                residual, niter, uval = utils.getSolutionData( prefixstr + "residualN=" + str(N) + "w=" + "{:.2f}".format(w) + "maxlevels=" + str(maxlevel) + ".txt", N )
                totaliter.append( niter[-1] )


            plt.yscale('log')
            plt.xlabel( "Number of Levels" )
            plt.ylabel( "Iterations until Convergence" )
            plt.plot( nlevels, totaliter, "-o" , label = "N = " + str(N) )

            plt.legend()
            plt.savefig( prefixstr + "w = {:.2f}".format(w) + "NlevelsVSNiterationsMG.png" )

def compareMGLevelsWithIterative():


    prefixstr = "/home/gaurav/plots/PoissonMG/"
    N = 49
    nlevels = [ 2, 3, 4 ] 
    wvals = [0.50, 0.8, 1.0]
    dirvals = [ "GS", "Jacobi", "MG" ]
    labels = [ "Gauss Siedel", "Jacobi", "Multigrid" ]

    for w in wvals:

        idx = 0
        plt.figure()

        for dirval in dirvals:

            if dirval == "MG":
                for maxlevel in nlevels:
                    residual, niter, uval = utils.getSolutionData( prefixstr + dirval + "/residualN=" + str(N) + "w=" + "{:.2f}".format(w) + "maxlevels=" + str(maxlevel) + ".txt", N )

                    # plt.yscale('log')
                    plt.xlabel( "Number of Iterations" )
                    plt.ylabel( "Residual" )
                    plt.plot( niter, residual, "-o" , label = "Multigrid, Cycles = " + str(maxlevel) )

            else:
                residual, niter, uval = utils.getSolutionData( prefixstr + dirval + "/residualN=" + str(N) + "w=" + "{:.2f}".format(w) + ".txt", N )

                # plt.yscale('log')
                plt.xlabel( "Number of Iterations" )
                plt.ylabel( "Residual" )
                plt.plot( niter, residual, "-o" , label = labels[idx] )

            idx += 1


        plt.legend()
        plt.savefig( prefixstr + "w = {:.2f}".format(w) + "N=" + str(N) + "NIterVSResidualForMGLevelsAndIterative.png" )


if __name__=="__main__":

    # compareMGLevelsWithIterative()
    # getBestW()
    compareConvergence()