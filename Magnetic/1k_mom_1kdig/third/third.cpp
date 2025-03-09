#include <string>
#include <boost/lexical_cast.hpp>

//Boost.Multiprecision mpfr_float
#include <boost/multiprecision/mpfr.hpp>

//for parallel implementation
#include <boost/mpi.hpp>
#include <boost/thread.hpp>

//Boost.Math headers
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

#define digits 1050

using namespace boost::archive;
using namespace boost::serialization;
using namespace boost::multiprecision;
using namespace boost::numeric::ublas;
using namespace std::chrono;

//serialization code for mpfr_float

using realtype = number<mpfr_float_backend<digits, allocate_stack>>;
#define MPFR_BUFFER_SIZE 1060

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

int i,j,k,m,n, l;
int betas = 30; // number of omegas to test
const int s = 100, b = 9; // l*(b+1)  = d+1
const int d = s*(b+1) - 1;
const realtype Pi = boost::math::constants::pi<realtype>();



/*vectors to contain the factors needed in the computation. We store them
so we dont have to recompute them each time we need them.*/
vector<realtype> factor(d+1);
vector<realtype> factor1(d+1);
vector<realtype> factor2(d+1);

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

    //--------------------factor------------------------------------------------------------

    for ( i = 0; i < factor.size(); ++i )
    {
        factor(i) = realtype(1)/ math::factorial < realtype > ( d - i );

        factor1(i) = pow( realtype(2), realtype(i + 1) ) * math::factorial <realtype> ( i );

        factor2(i) = pow(-realtype(1), realtype(i))/ pow(math::factorial <realtype> ( i ), realtype(2));
    }

    //-------------read-in the constants----------------------------------------------------

    vector<realtype> constants(d+1);
    std::ifstream infile;
    std::string st1;
    infile.open( "../Constants/Constant.txt" );
    for (i=0; i < constants.size(); ++i)
    {
        std::getline(infile, st1);
        constants(i) = realtype(st1);
    }

    infile.close();
    //----------------------------------------------------------------------------------
    // allocate the ublas::vectors containing the terms to be summed over k in the summation

    vector<realtype> third(floor((d-1)/2) + 1), thirdter(betas), total(floor((d-1)/2) + 1), term(d+1), term1(d+1);

    for ( int k = 0; k < floor((d-1)/2) + 1; ++k )
    {
        for( int l = 2*k+1; l < d + 1 ; ++l )
        {
            term( l ) = factor2(l) * factor1(l-2*k-1);
        }

        for( int m = 2*k+1; m < d + 1; ++m )
        {
            term1(m) = inner_prod( project(term, range(2*k+1, m+1 ) ),  project( factor, range (  d - m + 2*k + 1,  d + 1) ) );//(m-(2k+1))!...(0)!
        }

        third(k) = inner_prod( project ( constants, range ( 2*k + 1, d+1 ) ), project( term1, range ( 2*k + 1, d+1 ) ) );
    }

    std::ofstream outfile1;
    outfile1.open("../results/THIRD.txt",std::ios_base::out);
    outfile1.precision(digits);

    for ( i = 0; i < betas ; ++i )
    {
        for ( k = 0; k < floor((d-1)/2) + 1 ; ++k)
        {
            total(k) = pow( -realtype(1)/beta(i), realtype( k ) ) * third(k);
        }

        thirdter(i) = sum(total);
        outfile1 << std::scientific << thirdter(i) << std::endl;
    }


    return 0;

}
