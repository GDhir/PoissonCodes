#ifndef Gridval
#define Gridval

#include<utility>
#include<ConstantData.hpp>
#include<vector>

using std::vector;

class Grid {

    int neta;
    int nzeta;

    public:

        Grid( int, int );

        int getNeta();
        double getH();

        int getNzeta();
        int getNx();
        int getNy();

        std::pair<double, double> getXY( int i, int j, ConstantData<double>& cvals );
        std::pair<double, double> getXYRegular( int i, int j );

        double calcXFromIndex( int idx, int npoints );
        vector<double> interpolateVal( vector<double>& coarseVals, Grid& coarsegrid );
        vector<double> restrictVal( vector<double>& fineVals, Grid& coarsegrid );

        ~Grid() = default;

};

#endif