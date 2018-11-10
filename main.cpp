#include <bits/stdc++.h>
#include "/usr/include/eigen3/Eigen/Eigen"
#include "/usr/include/eigen3/Eigen/MatrixFunctions"
#include "meanShift.h"
//#include "meanShift.h"

using namespace std;
using namespace Eigen;
//using namespace MatrixFunctions

double EXP(double D){
    return exp(D);
}

int main()
{
    meanShift cluster;

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
        //cout<<x<<" "<<y<<endl;
        DATAx(0,n-N)=x;
        DATAy(0,n-N)=y;
    }
//cout<<DATAx<<endl<<DATAy<<endl;
  /*cout<<*/  cluster.Cluster(3,DATAx,DATAy);
    //Matrix<double,Dynamic,Dynamic> MATT; MATT.resize(n,n);
    //MATT=(((DATAx.replicate(n,1))-(DATAx.replicate(n,1)).transpose()).cwiseAbs()).triangularView<Upper>();
    //cout<<MATT<<endl;
    //MATT=(MATT.unaryExpr(&EXP)).triangularView<Upper>();
    //cout<<MATT.MATT.unaryExpr(&EXP)<<endl<<MATT.CwiseUnaryOp(exp(),);
    //cwiseInverse

    fclose (stdin);  // stop reading data from in.txt
    fclose (stdout); // stop writing data from out.txt
    return 0;
}
