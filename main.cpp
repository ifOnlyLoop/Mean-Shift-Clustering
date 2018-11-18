#include <bits/stdc++.h>
#include "/usr/include/eigen3/Eigen/Eigen"
#include "/usr/include/eigen3/Eigen/MatrixFunctions"
#include "MSC.h"
//#include "meanShift.h"
using namespace std;
using namespace Eigen;
//using namespace MatrixFunctions
double EXP(double D){
    return exp(D);
}
int main()
{
    MSC cluster;
    freopen("in.txt" , "r", stdin);
    freopen("out.txt", "w", stdout);
        Matrix<double,1,Dynamic> DATAx,DATAy;
    int N,n; cin>>N; n=N-1;
        DATAx.resize(NoChange,N);
        DATAy.resize(NoChange,N);
    double x,y;
    char c;
    while(N--){
        cin>>x>>c>>y;
     
        DATAx(0,n-N)=x;
        DATAy(0,n-N)=y;
    }

    cluster.cluster(3,DATAx,DATAy);
    fclose (stdin);  // stop reading data from in.txt
    fclose (stdout); // stop writing data from out.txt
    return 0;
}