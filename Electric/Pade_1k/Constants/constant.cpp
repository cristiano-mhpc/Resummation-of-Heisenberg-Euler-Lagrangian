#include <fstream>
#include <iostream>
#include "unsupported/Eigen/MPRealSupport"
#include "Eigen/LU"
#include "mpreal.h"
#include <chrono>

using namespace std::chrono;
using namespace mpfr;
using namespace Eigen;
using namespace std;


int main()
{
    auto start = high_resolution_clock::now();
    using mpfr::mpreal;


    int i,j,k,m,n;
    const int N = 999;
    const int M = N + 1;
    const int d = N + M ;
    const int digits = 3000; //for setting the precision


    /*Setup default precision for all subsequent computations
    MPFR accepts precision in bits - so we do the conversion*/
    const int bits = mpfr::digits2bits(digits);
    mpreal::set_default_prec(bits);

    //some mp literals
    const mpreal pi = mpfr::const_pi();


    //Declare matrix and vector types with multiprecision scalar types
    typedef Matrix<mpreal,Dynamic,Dynamic> MatrixXmp;
    typedef Matrix<mpreal,Dynamic,1>      VectorXmp;
    typedef Matrix<mpreal,1,Dynamic>   RowVectorXmp;


    //allocate and read-in the berno
    VectorXmp moments(d+1);

    std::string st;
    std::ifstream infile;
    infile.open("../../moments/moments.txt", std::ios_base::in);
    for (int i = 0; i < d+1; ++i )
    {
      std::getline(infile, st);
      moments(i) = mpreal(st);
    }

    infile.close();

    //Allocate matrices vetors
    MatrixXmp P(M,M);

    VectorXmp coeff(M);

    //populate the matrix and vector
    for (i = 1; i < M + 1; i++)
    {
        coeff(i-1) = -moments(N + i);

        for (j = 1; j < M + 1; j++)
        {
            P(i-1, j-1) = moments(N + i - j);
        }
    }

    //solve the system P constant = moments
    VectorXmp constant = P.partialPivLu().solve(coeff);

    //print out the constant
    ofstream outfile1;
    outfile1.open("Constant.txt",ios_base::out);
    outfile1.precision(digits);
    for (i = 0; i < constant.size(); i++)
    {
        outfile1 << scientific << constant(i) << endl;
    }
    outfile1.close();

    P.resize(0,0);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<minutes>(stop - start);
    std::cout << duration.count() << "minutes" << std::endl;

    return 0;

}
