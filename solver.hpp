
#pragma once
using namespace std;
#include <ostream>
#include <complex>
namespace solver
{
    class RealVariable
    {
    public:
        double co1;
        double co2;
        double co3;
        double maala;
        RealVariable(double c1 =0,double c2 = 1,double c3 = 0, double m = 1 ):co1(c1),co2(c2),co3(c3),
        maala(m){
            
        }

        friend RealVariable operator == (RealVariable a , RealVariable b);
        friend RealVariable operator == (RealVariable a , int b);

        friend  RealVariable operator + ( const RealVariable a , const RealVariable b );
        friend  RealVariable operator + (  RealVariable a ,  int b );

        friend  RealVariable operator - ( const RealVariable a , const  RealVariable b );
         friend  RealVariable operator - (  RealVariable a ,   int b );

        friend  RealVariable operator * ( const RealVariable a , const RealVariable b );
        friend const RealVariable operator ^ (  RealVariable b , const RealVariable a ){
        b.maala = a.maala;
        b.co1 = 1;
        b.co2 = 0;



        return b;
        };
        friend const RealVariable operator / ( const RealVariable a , const RealVariable b );
        friend const RealVariable operator / (  RealVariable a ,  int b );

        friend RealVariable operator + (int a , solver::RealVariable b );
        friend RealVariable operator * ( const int a , solver::RealVariable b );
        friend RealVariable operator - ( int a , RealVariable b );
        friend RealVariable operator^(const RealVariable b,int a );
        friend RealVariable operator / ( int a , const RealVariable b );

        friend RealVariable operator + (double a , RealVariable b );
        friend RealVariable operator * ( double a , RealVariable b );
        friend RealVariable operator - ( double a , RealVariable b );
        friend RealVariable operator ^ ( double a , const RealVariable b );
        friend RealVariable operator / ( double a , const RealVariable b );

        friend ostream& operator<< (ostream& os, const RealVariable& c);

    };

    class ComplexVariable
    {
    public:
        double co1;
        double co2;
        double co3;
        double maala;
        double im;
        ComplexVariable(double c1 =0,double c2 =1,double c3 = 0, double m = 1,double _im =0):co1(c1),co2(c2),co3(c3),
        maala(m),im(_im){}        
        ComplexVariable(int re) :co1(0), co2(re),co3(0),maala(1), im(re) {};
        ComplexVariable(complex<double> y): co1(0),co2(0),co3(y.real()),maala(1),im(y.imag()){};
        friend ComplexVariable operator == (ComplexVariable a , ComplexVariable b);
        friend const ComplexVariable operator + ( const ComplexVariable a , const ComplexVariable b );
        friend const ComplexVariable operator - ( const ComplexVariable a , const  ComplexVariable b );
        friend const ComplexVariable operator * ( const ComplexVariable a ,const ComplexVariable b );
        friend const ComplexVariable operator / (  const ComplexVariable a , const ComplexVariable b );
        friend ComplexVariable operator + (int a , ComplexVariable b );
        friend ComplexVariable operator * ( int a , ComplexVariable b );
        friend ComplexVariable operator - ( int a , ComplexVariable b );
        friend ComplexVariable operator ^ (ComplexVariable a , ComplexVariable b );
        friend ComplexVariable operator / ( int a ,  ComplexVariable b );

        friend ComplexVariable operator + (double a , ComplexVariable b );
        friend ComplexVariable operator * ( double a , ComplexVariable b );
        friend ComplexVariable operator - ( double a , ComplexVariable b );
        friend ComplexVariable operator^(ComplexVariable b,double a );
        friend ComplexVariable operator / ( double a ,  ComplexVariable b );

        friend ostream& operator<< (ostream& os, const ComplexVariable& c);
        friend ostream& operator<< (ostream& os, const std::complex<double>& c);

    };

    double solve(RealVariable a);

    std::complex<double>  solve(ComplexVariable y);

};
