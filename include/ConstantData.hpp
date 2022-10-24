#ifndef cdata
#define cdata

#include<cmath>

using std::pow;
using std::sqrt;

template<typename T>
struct ConstantData {

    T l;
    T h;
    T b;
    T a;

    ConstantData(T _l, T _h, T _b ):l(_l),h(_h),b(_b){

        a = sqrt( pow( ( l - b ), 2 )/4 - h*h );

    };

};

#endif