#ifndef solverval
#define solverval

#include<Grid.hpp>
#include<vector>
#include<cmath>
#include<ConstantData.hpp>
#include<utils.hpp>
#include<unordered_map>
#include<functional>
#include<string>
#include<utility>

using std::pow;
using std::vector;
using std::sqrt;
using std::pair;

using vecpair = vector<pair<double, int>>;
using cmptype = std::function<bool( double, double, int, int )>;
class Solver;

typedef std::vector<std::pair<double, int> > (Solver::*mytype)(Grid&, std::vector<double>&, std::vector<double>&, double, double, int, cmptype);

class Solver {

    public:

        std::unordered_map<std::string, mytype>  funcMap;

        Solver();

        void solveChannelFlow( Grid& grid, std::vector<double>& dataValue, const ConstantData<double>& cvals, const double tol );

        vector<pair<double, int>> solveJacobiPoisson(Grid &grid, std::vector<double> &dataValue, vector<double>& fval, const double w, const double tol, int maxiter, cmptype comparator);

        vector<pair<double, int>> computeMultigrid( Grid refgrid, std::vector<double> &dataValue, vector<double> &fval, const int maxlevels,
                                      const double w, const double tol, string method );

        void multigridVCycle(Grid refgrid, std::vector<double> &dataValue, vector<double> &fval, int currlevel,
                                      const double w, const double tol, string method, int& totalIter, int niter );

        vector<pair<double, int>> solveGSPoisson(Grid &grid, std::vector<double> &dataValue, vector<double> &fval, const double w,
                     const double tol, int maxiter, cmptype comparator );

        ~Solver(){};

};

#endif