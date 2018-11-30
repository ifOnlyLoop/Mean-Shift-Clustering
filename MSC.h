#include <bits/stdc++.h>
#include "/usr/include/eigen3/Eigen/Eigen"
#include "/usr/include/eigen3/Eigen/MatrixFunctions"

// using eigen is a must here

using namespace std;
using namespace Eigen;

// since the matrix width is dynamic
typedef Matrix<double,1,Dynamic> MAT1D;
typedef vector<pair<double,double>> VP;
class MSC{

private:

        static double BW; 
        double Wsum, Dmax=0.00001, x, y, xi, yi,xo, yo,K;
        VP OUT;
        MAT1D TEMP1Dx,TEMP1Dy, TEMP1D,X,Y, W;
    
        /*
         * BW : band width
         * Wsum : sum of weights for after one shift iteration
         * Dmax : exit condition (point reached cluster centeroid)
         * TEMP1D_ : temporary vector
         * X : point shift in x axis
         * Y : point shift in y axis
         * W : weights matrix 
         */
    
protected:

        static double SQ2(double D); // return x^0.5 for input x
        static double PO2(double D); // return x^2 for input x
        static double KF(double D);  // return kernel wieght for input x
        double  SHIFT(MAT1D &DATA);  // shift point to cluster
    
public:

        VP cluster(double bw, MAT1D &DATAx, MAT1D &DATAy); // cluster data

};

VP MSC::cluster(double bw, MAT1D &DATAx, MAT1D &DATAy){
    
    // extract bandwidth and data size
    BW = bw; int N = DATAx.cols(), n=N-1;
    
    // create temp vectors and resize matrices 
    TEMP1D.resize(NoChange,N); 
    TEMP1D.setOnes(); // eigen3 can't sub a number (a neat way to do it) 
    
    TEMP1Dx.resize(NoChange,N);
    TEMP1Dy.resize(NoChange,N);

    W.resize(NoChange,N);
    X.resize(NoChange,N);
    Y.resize(NoChange,N);
    

    while(N--){ //
        
        // select point to shift
        xo=DATAx(0,n-N);
        yo=DATAy(0,n-N);
        K=1;
        
        do{
            // find differance squared
            if(K){xi=xo;yi=yo;K=0;}else{xi=x;yi=y;}
            
            // find axes diff^2
            TEMP1Dx=(DATAx-xi*TEMP1D).unaryExpr(&PO2);
            TEMP1Dy=(DATAy-yi*TEMP1D).unaryExpr(&PO2);
            
            // find weights
            W = (TEMP1Dx+TEMP1Dy).unaryExpr(&SQ2).unaryExpr(&KF);
            Wsum = W.sum();
            
            // shift points
            x=SHIFT(DATAx); //X(0,n-N)=x;
            y=SHIFT(DATAy); //Y(0,n-N)=y;
            
        }while(PO2(x-xi)+PO2(y-yi)>Dmax);
        
        OUT.push_back({x,y});
        //printf("%f %f -> %f %f\n",xo,yo,x,y);
    }

    return OUT;
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
