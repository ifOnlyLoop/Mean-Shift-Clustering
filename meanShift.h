#include<bits/stdc++.h>
#include"/usr/include/eigen3/Eigen/Eigen"
#include"/usr/include/eigen3/Eigen/MatrixFunctions"

using namespace std;
using namespace Eigen;

typedef std::pair<double,double> PofD;
typedef std::vector<PofD> VofP;
typedef VofP::iterator VofPit;
typedef Matrix<double,Dynamic,Dynamic> MATDD;
typedef Matrix<double,1,Dynamic> MAT1D;


class meanShift {

private:

    int N;
    
    MAT1D xTEMP1D, yTEMP1D, Wsum;
    MATDD xTEMPDD, yTEMPDD, TEMPDD, W;
    
    double PI = tan(-1), 
           e=2.718281828459;

    static double h/*bandwidth*/;
    
    double Xshift = 0,Yshift = 0, scaleFactor =0, // para
           XshOld = 0,YshOld = 0, Dshift = 0, Eps = 0.5; // para
    
    Matrix<double,Dynamic,Dynamic> diffMATx,diffMATy,temp;

protected:

    static double euclideanDist_cwise(double D);
    static double KF_cwise(double D);
    static double Inv(double D);

    MAT1D Shift(MAT1D &DATA);   
    MAT1D W_cwise();
    MATDD diffMAT(MAT1D &DATA);
    MATDD euclideanDist(); // distance

public:

    MATDD Cluster(double bandwidth, MAT1D &DATAx,MAT1D &DATAy);



};
/*
 *
 * 
 * 
 * 
 * 
 * 
 * 
 */

double meanShift::h = 0;

MATDD meanShift::diffMAT(MAT1D &DATA){  /*
                                         * this function return the difference matrix 
                                         * between the elements of an Eigen vector
                                         */
    
    return ((DATA.replicate(N,1))-(DATA.replicate(N,1)).transpose()).cwiseAbs(); //.triangularView<Upper>()
    
}

MATDD meanShift::euclideanDist(){ /*
                                   * this function return the distance matrix 
                                   * between the elements of an Eigen vector
                                   */
    
    return (diffMATx.unaryExpr(&euclideanDist_cwise)+diffMATy.unaryExpr(&euclideanDist_cwise)).cwiseSqrt();
    
}

double meanShift::euclideanDist_cwise(double D){ /*
                                                  * return the square of a number
                                                  */

    return (D*D);

}

double meanShift::KF_cwise(double D){ /*
                                       * return the kernel function value cwise
                                       */
    return std::exp(-0.5*D/h*D/h);
    
} 

MAT1D meanShift::W_cwise(){ /*
                             * return the kernel function value cwise
                             */
    W = (euclideanDist()).unaryExpr(&KF_cwise);
    return W.colwise().sum(); //Wsum
    
}

double meanShift::Inv(double D){ /*
                                       * return the kernel function value cwise
                                       */
    
    return 1/D;
    
} 

MAT1D meanShift::Shift(MAT1D &DATA){

    TEMPDD = W.cwiseProduct(DATA.replicate(N,1));
    return   (TEMPDD.rowwise().sum()).cwiseProduct(Wsum.cwiseInverse().transpose());
}

MATDD meanShift::Cluster(double bandwidth, MAT1D &DATAx,MAT1D &DATAy){

    /*
     * find difffence matrix
     * for x and y
     */

    h=bandwidth;
    N=min(DATAx.cols(),DATAx.cols()); // add a part later to clip extra elements

//  resize
    xTEMP1D.resize(NoChange,N); TEMPDD.resize(N,N);    //diffMATx.resize(N,N); 
    yTEMP1D.resize(NoChange,N); //yTEMPDD.resize(N,N);    //diffMATy.resize(N,N); 
//  differance    
    diffMATx  = diffMAT(DATAx);
    diffMATy  = diffMAT(DATAy);//cout<<diffMATx;
//  weights    
    Wsum = W_cwise();
//  shift

    xTEMP1D=Shift(DATAx);
    yTEMP1D=Shift(DATAy);
    
    cout<<xTEMP1D<<endl<<yTEMP1D;
    
    
    return diffMATy;
}



/*
 *
 * 
 */

/*
 *
 * 
 */
