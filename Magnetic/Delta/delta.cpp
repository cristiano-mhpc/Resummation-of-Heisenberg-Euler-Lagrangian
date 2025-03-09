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
#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/math/special_functions/binomial.hpp>

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
#include <stdlib.h>
#include <time.h>

#define digits 3000

using namespace boost::archive;
using namespace boost::serialization;
using namespace boost::multiprecision;
using namespace boost::numeric::ublas;
using namespace std::chrono;

//serialization code for mpfr_float

using realtype = number<mpfr_float_backend<digits, allocate_stack>>;
#define MPFR_BUFFER_SIZE 3000

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
namespace math = boost::math;

int i,j,k,m,n, l;
int betas = 30; // number of betas to test
const int s = 10, b = 9; // l*(b+1)  = d+1//
//const int d = s*(b+1) - 1;

const realtype Pi = boost::math::constants::pi<realtype>();



int main(int argc, char** argv)
{
    if(argc < 1){
      std::cout << "Usage: input the number of coefficients to use." << std::endl;
      return 1;
    }
   
    int d = std::stoi(argv[1]);

    //ublas vectors
    vector<realtype> terms(d+1);
    vector<realtype> moments(d+1);
    vector<realtype> summ(d+1); //stores the partial sum

    vector<realtype> numerator_terms(d);
    vector<realtype> denominator_terms(d);

    vector<realtype> delta_sums(betas);
    vector<realtype> beta(betas);

    for(i=0; i< betas-2; ++i)
    {
        beta(i) = pow(realtype(10), realtype(i-5));
    }
    beta(betas-2) = realtype(1)/realtype(5);
    beta(betas-1) = realtype(4);


   // std::ofstream outfile;
   // outfile.open("moments.txt",std::ios_base::out);
   // outfile.precision(digits);
    //-------------the moments ------------------------------------
    for ( int i = 0; i < d+1; ++i )
    {
        moments(i) = ( pow(realtype(-1), realtype( i )) * (realtype(2)-pow(realtype(4), realtype( i+2 ))) * math::bernoulli_b2n <realtype> ( i + 2 ) )
                        /realtype( ( 2 * ( i ) + 4 ) * ( 2 * ( i ) + 3 ) * ( 2 * ( i ) + 2 )  );

     //   outfile << std::scientific << moments(i) << std::endl;

    }

    //----------------------------------------------------------------------------------
    for (i = 0; i < beta.size();++i)
    {
        for (l = 0; l < terms.size(); ++l)
        {
            terms(l) = moments(l)*pow(-beta(i), l);//terms
        }

        for(n = 0; n < terms.size(); ++n)
        {
            summ(n) = sum(subrange(terms,0,n+1));//partial sums
        }

        for (j = 0; j < numerator_terms.size();++j)
        {
            numerator_terms(j) = pow( realtype(-1), realtype( j ) ) * math::binomial_coefficient<realtype> (d-1, j)
                                *(math::rising_factorial<realtype> (1+j,d-2)/math::rising_factorial<realtype> (d,d-2))
                                *(summ(j) / terms(j+1) );

           denominator_terms(j) = pow(realtype(-1), realtype( j ) ) * math::binomial_coefficient<realtype> (d-1, j)
                                *(math::rising_factorial<realtype> (1+j,d-2)/math::rising_factorial<realtype> (d,d-2))
                                *( realtype(1) / terms(j+1) );
        }


        delta_sums(i) = pow(beta(i),realtype(2))*sum(numerator_terms)/sum(denominator_terms);


    }
    
    std::string filename = "delta_" + std::string(argv[1]) + ".txt";
    
    std::ofstream outfile1(filename);
    //outfile1.open("delta_100.txt",std::ios_base::out);
    //std::ofstream outfile1;
    //outfile1.open("delta_100.txt",std::ios_base::out);
    outfile1.precision(digits);

    for(n = 0; n < delta_sums.size(); ++n)
    {
        outfile1 << std::scientific << delta_sums(n) << std::endl;
    }

    return 0;

}
