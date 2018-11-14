#include <bits/stdc++.h>
#include "/usr/include/eigen3/Eigen/Eigen"
#include "/usr/include/eigen3/Eigen/MatrixFunctions"

using namespace std;
using namespace Eigen;

typedef Matrix<double,1,Dynamic> MAT1D;

class MSC{

private:

    static double BW; 
    double Wsum, Dmax=0.00001, x, y, xi, yi,xo, yo,K;
    MAT1D TEMP1Dx,TEMP1Dy, TEMP1D,X,Y, W;
    
protected:

    static double SQ2(double D);
    static double PO2(double D);
    static double KF(double D);
    double  SHIFT(MAT1D &DATA);
public:

    void cluster(double bw, MAT1D &DATAx, MAT1D &DATAy);

};

void MSC::cluster(double bw, MAT1D &DATAx, MAT1D &DATAy){

    BW = bw; int N = DATAx.cols(), n=N-1;
    
    TEMP1D.resize(NoChange,N); 
    TEMP1D.setOnes();
    
    TEMP1Dx.resize(NoChange,N);
    TEMP1Dy.resize(NoChange,N);
    vector<bool> isClustered;
    
    for(int i=0;i<N;i++){
        isClustered.push_back(false);
    }

    W.resize(NoChange,N);
    X.resize(NoChange,N);
    Y.resize(NoChange,N);
    
//do
    while(N--){
        xo=DATAx(0,n-N);
        yo=DATAy(0,n-N);
        K=1;
        do{
            //if((n-N)!=467) continue;
            // exit condition for one point
            if(isClustered[n-N]) continue;
            // find differance squared
            if(K){xi=xo;yi=yo;K=0;}else{xi=x;yi=y;}
            //cout<<xi<<" "<<yi<<"\n";
            TEMP1Dx=(DATAx-xi*TEMP1D).unaryExpr(&PO2);
            TEMP1Dy=(DATAy-yi*TEMP1D).unaryExpr(&PO2);
            // find weights
            W = (TEMP1Dx+TEMP1Dy).unaryExpr(&SQ2).unaryExpr(&KF);
            Wsum = W.sum(); 
            // shift points
            x=SHIFT(DATAx); X(0,n-N)=x;
            y=SHIFT(DATAy); Y(0,n-N)=y;
        }while(PO2(x-xi)+PO2(y-yi>Dmax));
        
        printf("%f %f -> %f %f\n",xo,yo,x,y);
    }
//while

    return;
}

double MSC::BW=0;

double MSC::PO2(double D){
    return D*D;
}

double MSC::SQ2(double D){
    return std::sqrt(D);
}

double MSC::KF(double D){
    return std::exp(-0.5*D/BW*D/BW);
}

double MSC::SHIFT(MAT1D &DATA){
    return (W*DATA.transpose()).sum()/Wsum;
}
