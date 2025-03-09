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

#define digits 100

using namespace boost::archive;
using namespace boost::serialization;
using namespace boost::multiprecision;
using namespace boost::numeric::ublas;
using namespace std::chrono;

//serialization code for mpfr_float

using realtype = number<mpfr_float_backend<digits, allocate_stack>>;
#define MPFR_BUFFER_SIZE 110

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
int betas = 30; // number of omegas to test

const int s = 1, b = 4; // l*(b+1)  = d+1
const int d = s*(b+1) - 1;
int num_terms = 500;

const realtype Pi = boost::math::constants::pi<realtype>();
const realtype loga = boost::math::log1p( -realtype(1)/realtype(2) );

vector<realtype> factor(d+1); // lifetime of factor spans the entire program
vector<realtype> constants(d+1);

vector<realtype> factor1( 2*num_terms + 2);
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

    //--------------------factor--------------------------------------------

    for (i=0; i < factor.size(); ++i)
    {
        factor(i) = realtype(1)/math::factorial < realtype > ( d - i ); //factorial

        factor2(i) =  realtype(1)/pow( math::factorial < realtype > ( i ), realtype(2) );
    }


    for (i = 0; i < factor1.size() ; ++i)
    {
        factor1(i) =  ( pow( realtype(1)/realtype(2), realtype( i ) ) * ( loga - math::digamma(realtype( i + 1  ) ) ) )  /
                        ( boost::math::factorial < realtype > (  i   ) );
    }

    //-------------read-in the constants------------------------------------

    std::ifstream infile;
    std::string st1; // we could use std::string
    infile.open("../Constants/Constant.txt");
        for ( i=0; i < constants.size(); ++i )
    {
        std::getline(infile, st1);
        constants(i) = realtype(st1);
    }

    infile.close();
    //----------------------------------------------------------------------------------
    // allocate the ublas::vectors containing the terms to be summed over k in the summation

    vector<realtype> fourth(num_terms - floor((d-1)/2) ),fourthter(betas), term(d+1), term1(d+1), total(num_terms - floor((d-1)/2) );

    for (k = 0; k < num_terms - floor((d-1)/2) ; ++k)
    {
        for( l = 0; l < d+1; ++l)
        {
            term(l) =  factor2(l)*factor1(2*(k + floor((d-1)/2) + 1)  - l );
        }

        for( m = 0; m < d+1; ++m)
        {
            term1(m) = inner_prod( project( term, range(0, m+1) ), project( factor, range (d - m , d + 1) ) );
        }

        fourth( k ) = inner_prod( constants, term1 );
    }


    std::ofstream outfile1;
    outfile1.open( "../results/FOURTH.txt",std::ios_base::out );
    outfile1.precision(digits);

    for (i = 0; i < betas ; ++i )
    {
        for ( k = 0; k < total.size(); ++k )
        {
            total( k ) = pow( -realtype(1)/ beta(i), realtype(k + floor((d-1)/2 + 1) ) ) * fourth(k);
        }

        fourthter(i) = -sum(total);

        outfile1 << std::scientific << fourthter(i) << std::endl;
    }


    return 0;

}
