#include <string>
#include <boost/lexical_cast.hpp>

//Boost.Multiprecision mpfr_float
#include <boost/multiprecision/mpfr.hpp>

//for parallel implementation
#include <boost/mpi.hpp>
#include <boost/thread.hpp>

//Boost.Math headers
#include <boost/math/special_functions/log1p.hpp> //computes the ln(1+x)
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>

//Boost.Serialization headers
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/version.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/array.hpp>
#include <stdio.h> // used in creating File object

//Boost.Ublas headers
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

//for serialization and miscellaneous functionalities
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include<stdlib.h>
#include<time.h>

#define digits 1550

using namespace boost::archive;
using namespace boost::serialization;
using namespace boost::multiprecision;
using namespace boost::numeric::ublas;
using namespace std::chrono;

//serialization code for mpfr_float

using realtype = number<mpfr_float_backend<digits, allocate_stack>>;
#define MPFR_BUFFER_SIZE 1560

    namespace boost {
        namespace serialization {
            template<class Archive>
            void save(Archive& ar, const realtype& x, const boost::serialization::version_type&) {
                static char buffer[MPFR_BUFFER_SIZE];
                FILE* fid = fmemopen(buffer, MPFR_BUFFER_SIZE, "wb+");
                mpfr_fpif_export(fid, const_cast<mpfr_ptr>(x.backend().data()));
                fseek(fid, 0L, SEEK_END);
                long length = ftell(fid);
                ar& length;
                ar& boost::serialization::make_array(buffer, length);
                fclose(fid);
            }

            template<class Archive>
            void load(Archive& ar, realtype& x, const boost::serialization::version_type&) {
                static char buffer[MPFR_BUFFER_SIZE];
                long length = 0;

                ar& length;
                ar& boost::serialization::make_array(buffer, length);

                FILE* fid = fmemopen(buffer, length, "r");
                mpfr_fpif_import(x.backend().data(), fid);
                fclose(fid);
            }

            template<class Archive>
            inline void serialize(Archive& ar, realtype& t, const unsigned int file_version) {
                split_free(ar, t, file_version);
            }
        }
    }

namespace mpi = boost::mpi;

int i,j,k,m,n,l;
int betas = 30; // number of betas to test
const int s = 150, b = 9; // l*(b+1)  = d+1
const int d = s*(b+1) - 1;
const realtype Pi = boost::math::constants::pi<realtype>();
const realtype loga = boost::math::log1p( -realtype(1)/realtype(2) );

vector<realtype> factor(d+1);
vector<realtype> factor1(2*floor((d-1)/2) + 1);
vector<realtype> factor2(2*floor((d-1)/2) + 1);

std::string st1; //for reading- in constant

int main()
{
    namespace math = boost::math;

    vector<realtype> beta(betas);

    for (i = 0; i < betas - 2; ++i)
    {
        beta(i) = pow( realtype(10), realtype( i - 5 ) );
    }
    beta(betas-2) = realtype(1)/realtype(5);
    beta(betas-1) = realtype(4);

    //--------------------factors--------------------------------------------
    for (i=0; i < factor.size (); ++i)
    {
        factor(i) = realtype(1)/ boost::math::factorial < realtype > ( d - i  ); //factorial

    }


    for (i = 0; i < factor1.size (); ++i)
    {

        factor1(i) = ( pow( realtype(1)/realtype(2), realtype( i  ) ) * ( loga - math::digamma(realtype( i + 1  ) ) ) )/
                        ( boost::math::factorial <realtype> (  i  ) );

        factor2(i) =  realtype(1) /  pow ( math::factorial < realtype > ( i ), realtype(2) );

    }


    //-------------read-in the constants------------------------------------
    std::ifstream infile;
    vector<realtype> constants(d+1);

    infile.open("../Constants/Constant.txt");
    for (i=0; i < constants.size(); ++i)
    {
        std::getline(infile, st1);
        constants(i) = realtype(st1);
    }

    infile.close();
    //----------------------------------------------------------------------------------
    // allocate the ublas::vectors containing the terms to be summed over k in the summation

    vector<realtype> second(floor((d-1)/2) + 1), secondter(betas), term1(d+1), total(d+1), term(2*floor((d-1)/2) + 1);

    for (int k = 0; k < floor((d-1)/2) + 1; ++k)
    {
        for(int l = 0; l < 2*k+1; ++l)
        {
            term(l) = factor2(l)*factor1(2*k - l ); //factor1(2k+1)...factor1(0), for any k
        }

        for(int m =2*k+1; m < d+1; ++m)
        {
            term1(m) = inner_prod( project(term, range(0,2*k+1) ), project( factor, range(d-m , d - m + 2*k + 1) ) ); //1/m!, 1/(m+1)!,...1/(m - 2k)
        }

        second(k) = inner_prod( project ( constants, range (2*k + 1, d+1) ), project ( term1 , range (2*k+1, d+1) ) );
    }

    std::ofstream outfile1;
    outfile1.open("../results/SECOND.txt",std::ios_base::out);
    outfile1.precision(digits);

    for ( i = 0; i < betas ; ++i )
    {
        for ( k = 0; k < floor((d-1)/2) + 1; ++k )
        {
            total(k) = pow(-realtype(1)/beta(i), realtype(k) ) * second(k);
        }

        secondter(i) = -sum(total);
        outfile1 << std::scientific << secondter(i) << std::endl;
    }


    return 0;

}
