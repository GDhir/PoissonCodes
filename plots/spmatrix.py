from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve, splu
import utils

def computeSolution( N, h, a, b ):

    dzeta = 1/( N - 1 )

    tot_size = (N - 1)*( N - 1 )

    Nstride = N - 1
    Ny = N - 1

    A, f = utils.computeInnerSecondOrder( N, Nstride, Ny, tot_size, a, h, b, dzeta )

    A = utils.computeLeftBoundaryFourthOrder( A, Ny, Nstride )

    A = utils.computeUpperBCThirdOrder( A, Nstride, Ny, N, a, b, tot_size )

    A = A.tocsr()
    u = spsolve( A, f )

    u = u.reshape( ( N - 1, N - 1 ) )

    u = np.append( u, np.zeros( (N - 1, 1) ), axis=1 )
    lowervals = np.zeros( ( 1, N ) )
    u = np.append( lowervals, u, axis=0 )

    return u


def computePareto(  ):

    from scipy.spatial import ConvexHull, convex_hull_plot_2d

    l = 3
    hvals = np.arange( 0.1, 1.1, 0.1 )
    bvals = np.arange( 0, 1.1, 0.1 )
    N = 81

    sumvals = []
    Ivals = []

    QIinfo = []

    for h in hvals:
        for b in bvals:

            a = sqrt( ( ( l - b )**2 )/4  - h**2  )

            u = computeSolution( N, h, a, b )
            #X, Y = utils.computeXY( N, h, a, b )

            sum = utils.computeSimpsons( u, N, a, h, b )
            sumvals.append( sum )

            I = utils.computeMomentofInertia( h, b, l )
            Ivals.append(I)

            QIinfo.append( [ I, sum ] )

    QIinfo = np.array( QIinfo )
    hull = ConvexHull( QIinfo )

    plt.plot( QIinfo[ hull.vertices, 0 ], QIinfo[ hull.vertices, 1 ], 'r--', lw=2 )
    plt.plot( QIinfo[ :, 0 ], QIinfo[ :, 1 ], 'ro')
    plt.show()


def computeConvergence():
    l = 3
    h = 1
    bvals = [0.0, 0.5, 1.0]

    for b in bvals:

        a = utils.computea( h, b, l)

        Nvals = [1001]

        iter = 0
        l2error = []
        sumref = 0
        uref = np.zeros( ( Nvals[0], Nvals[0] ) )
        sumerror = []
        linferror = []
        deltax = []

        prefixstr = "/home/gaurav/plots/MatsolveData/" + "plotb=" + str(b) + "/"

        for N in Nvals:
            
            u = computeSolution( N, h, a, b )
            X, Y = utils.computeXY( N, h, a, b )

            sum = utils.computeSimpsons( u, N, a, h, b )
            print(sum)

            filename = prefixstr + "udataN=" + str(N)
            utils.saveSolutionData( filename, u )


            plt.figure()
            plt.pcolormesh( X, Y, u )
            plt.colorbar()
            titledata = "Solution plot for b=" + str(b) + " and N=" + str(N)
            plt.title( label=titledata )
            filename = prefixstr + "umeshN=" + str(N)
            plt.savefig( filename + ".png" )
            # plt.close()
            #plt.show()

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

        # l2error = np.array( l2error )
        # linferror = np.array( linferror )
        # sumerror = np.array( sumerror )
        # deltax = np.array( deltax )

        # plt.figure()
        # plt.yscale('log')
        # plt.xscale('log')
        # plt.plot( deltax, l2error, "bo-", label = "L2 Error Convergence" )
        # plt.plot( deltax, deltax**2, label = " ${\Delta x}^2$ " )
        # plt.legend()
        # titledata = "L2 Error Convergence plot for b=" + str(b)
        # plt.title( label=titledata )
        # filename = prefixstr + "l2error"
        # plt.savefig( filename + ".png" )
        # # plt.close()
        # utils.savedata( filename + ".txt", l2error )

        # plt.figure()
        # plt.yscale('log')
        # plt.xscale('log')
        # plt.plot( deltax, linferror/np.max(linferror), "bo-", label = "Linf Error Convergence" )
        # plt.plot( deltax, deltax**2,label = " ${\Delta x}^2$ ")
        # plt.legend()
        # titledata = "Linf Error Convergence plot for b=" + str(b)
        # plt.title( label=titledata )
        # filename = prefixstr + "linferror"
        # plt.savefig( filename + ".png" )
        # # plt.close()
        # utils.savedata( filename + ".txt" , linferror)


        # print( deltax )
        # print( l2error )
        # print( linferror )
        # print( Nvals[0]*linferror )

        # plt.figure()
        # plt.yscale('log')
        # plt.xscale('log')
        # plt.plot( deltax, sumerror, "bo-", label = "Flowrate Error Convergence" )
        # plt.plot( deltax, deltax**2, label = " ${\Delta x}^2$ ")
        # titledata = "Flowrate Error Convergence plot for b=" + str(b)
        # plt.title( label=titledata )
        # plt.legend()
        # filename = prefixstr + "flowrateerror"
        # plt.savefig( filename + ".png" )
        # # plt.close()
        # utils.savedata( filename + ".txt", sumerror )
        # #plt.show()
        # plt.pause(3)
        # plt.close("all")


def LUSolver():

    l = 3
    prefixstr = "/home/gaurav/plots/"
    h = 1
    b = 0.5
    a = utils.computea( h, b, l )
    N = 7

    dzeta = 1/( N - 1 )

    tot_size = (N - 1)*( N - 1 )

    Nstride = N - 1
    Ny = N - 1

    A, f = utils.computeInnerSecondOrder( N, Nstride, Ny, tot_size, a, h, b, dzeta )

    A = utils.computeLeftBoundarySecondOrder( A, Ny, Nstride )

    A = utils.computeUpperBCForwardDiff( A, Nstride, Ny, N, a, b, tot_size )
    X, Y = utils.computeXY( N, h, a, b )

    from scipy.sparse.linalg import eigs
    from scipy.linalg import lu

    A = A.toarray()

    # u = np.triu( -A )

    # np.fill_diagonal(u, 0)
    # R = np.matmul(np.linalg.inv(A), u) + np.eye( Nstride*Nstride )

    # print(R)
    # vals, vecs = eigs(R, k=20)
    #print(vals)

    # plt.figure()
    # plt.plot(vals)
    # plt.show()

    plt.figure()

    # A = A.toarray()

    plt.imshow(np.where(A != 0, 1., 0), interpolation='None')
    plt.xlabel("Coefficients")
    plt.ylabel("Equations")
    plt.savefig(prefixstr + "Amatrix")
    plt.show()

    # print(A)

    # plt.show()
    #B = splu( A )
    #u = B.solve( f )

    # u = u.reshape( ( N - 1, N - 1 ) )

    # u = np.append( u, np.zeros( (N - 1, 1) ), axis=1 )
    # lowervals = np.zeros( ( 1, N ) )
    # u = np.append( lowervals, u, axis=0 )

    # plt.figure()
    # plt.pcolormesh( X, Y, u )
    # plt.colorbar()
    # plt.show()





# computePareto()
#computeConvergence()
LUSolver()