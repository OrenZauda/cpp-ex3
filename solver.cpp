
using namespace std;
#include <iostream>
#include <complex>
#include "solver.hpp"
#include <cmath>
#include <bits/exception.h>
namespace solver {

    RealVariable operator==(const RealVariable a, const RealVariable b) {
        RealVariable t(a.co1,a.co2,a.co3,a.maala);
        t.co3=a.co3-b.co1;
        cout<<"======== r r"<<endl;
        cout<<t.co1<<endl;
        cout<<t.co2<<endl;
        cout<<t.co3<<endl;
        cout<<t.maala<<endl;
        return t;
    };
  

    RealVariable operator+( RealVariable b, const RealVariable a) {
        cout<<"+++++++++++++++++a,second is "<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        cout<<endl;
        cout<<"b is"<<endl;
        b.co3+=a.co3;
        b.co2+=a.co2;
        b.co1+=a.co1;
        cout<<b.co1<<endl;
        cout<<b.co2<<endl;
        cout<<b.co3<<endl;
        return b;
    }
    RealVariable operator+( RealVariable b, double a) {
        b.co3+=a;
        return b;
    }
    RealVariable operator+( RealVariable b, int a) {
        cout<<"++++++++++++"<<endl<<endl;
        b.co3+=a;
        cout<<b.co1<<endl;
        cout<<b.co2<<endl;
        cout<<b.co3<<endl;

        return b;
    }
    
    RealVariable operator+( double a, RealVariable b) {
        b.co3+=a;
        return b;
    }
    RealVariable operator-(const RealVariable a,const  RealVariable b) {
        RealVariable t(a.co1,a.co2,a.co3,a.maala);

       t.co3=a.co3-b.co1;
        cout<<"---------------- r r"<<endl;
        cout<<t.co1<<endl;
        cout<<t.co2<<endl;
        cout<<t.co3<<endl;
        cout<<t.maala<<endl;
    return a;
    };
        RealVariable operator-(RealVariable a, int b) {
        a.co3 = a.co3- b;
        cout<<"--------------- r d"<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        cout<<a.maala<<endl;
        return a;
    }

    RealVariable operator*(const RealVariable a, const RealVariable b) {
        RealVariable t(a.co1,a.co2,a.co3,a.maala);
        t.co2 = b.co1;
        return t;
    }
    const RealVariable operator/(const RealVariable a, const RealVariable b) {
        RealVariable t(a.co1,a.co2,a.co3,a.maala);
        t.co3= t.co3/b.co1;
        return a;

    }

    RealVariable operator+(int b, RealVariable a) {
        cout<<"++++++++++++++++++"<<endl;
        
        RealVariable t(a.co1,a.co2,a.co3,a.maala);
        t.co3= t.co3+b;
        cout<<t.co1<<endl;
        cout<<t.co2<<endl;
        cout<<t.co3<<endl;
        cout<<t.maala<<endl;

        return t;
    };

    RealVariable operator-(int b, RealVariable a) {
        RealVariable t(a.co1,a.co2,a.co3,a.maala);
        t.co3-=b;
        return t;
    }


    RealVariable operator*(const int a, RealVariable &b) {
        cout<<"***********"<<endl;
        if(b.maala==1){
        b.co1=1;
        }
        
        b.co2 = a;
        cout<<b.co1<<endl;
        cout<<b.co2<<endl;
        cout<<b.co3<<endl;
        cout<<b.maala<<endl;

        return b;
    }

    RealVariable operator^( RealVariable b,int a ) {
        b.maala = a;
        b.co1 = 1;
        b.co2 = 0;

        cout<<"^^^^^^^^^^^^^"<<endl;

        cout<<b.co1<<endl;
         cout<<b.co2<<endl;
        cout<<b.co3<<endl;
         cout<<b.maala<<endl;
        return b;
    }
    

    RealVariable operator/(int b, const RealVariable a) {
        RealVariable t(a.co1,a.co2,a.co3,a.maala);
        t.co3=t.co3-b;
        return t;

    }
    

    

    RealVariable operator-(double a, RealVariable b) {
        b.co3-=a;
        return b;
    }

    RealVariable operator*(double a, RealVariable b) {
        cout<<"****************a"<<endl;
        if(b.maala==1){
            b.co2=1;
            b.co2*=a;
        }
        if(b.maala==2){
            b.co2*=a;
        }
        return b;
    }

    RealVariable operator^(double a, const RealVariable b) {

        return b;
    }

    RealVariable operator/(double a, const RealVariable b) {
        return b;
    }

    ostream &operator<<(ostream &os, const RealVariable &c) {
        return (os << "x");
    }

//*****************************

    ComplexVariable operator==(ComplexVariable a, ComplexVariable b) {
        
        cout<<"==========a"<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        cout<<"imm"<<a.im<<endl;
        cout<<a.maala<<endl;

        cout<<"bbbbbbbbbbbbbbbb?"<<endl;
        cout<<b.co1<<endl;
        cout<<b.co2<<endl;
        cout<<b.co3<<endl;
        cout<<"imm"<<b.im<<endl;
        cout<<b.maala<<endl;

        a.co3-=b.co2;
        return a;
    };


    const ComplexVariable operator+( ComplexVariable a, const ComplexVariable b) {
        

        cout<<"++++++++++++++++a"<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        cout<<a.maala<<endl;
        cout<<"b"<<endl<<endl;

        cout<<b.co1<<endl;
        cout<<b.co2<<endl;
        cout<<b.co3<<endl;
        cout<<b.maala<<endl;

        cout<<a.im<<endl;
        a.im = b.im;
        a.co1+=b.co1;
        a.co2+=b.co2;
        a.co3+=b.co3;
        cout<<"add"<<endl<<endl;;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        cout<<a.maala<<endl;
        return a;
    }

    const ComplexVariable operator-( ComplexVariable a, const ComplexVariable b) {
       a.im=b.im;
       cout<<"-------------------"<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        cout<<a.maala<<endl;
        cout<<a.im<<endl;
        return a;
    };

    const ComplexVariable operator*(const ComplexVariable a, const ComplexVariable b) {
        ComplexVariable z;
        return z;
    }


    const ComplexVariable operator/(const ComplexVariable a, const ComplexVariable b) {
       return a;

    }

    ComplexVariable operator+(int a, ComplexVariable b) {
        ComplexVariable z;
        return z;
    };

    ComplexVariable operator-(int a, ComplexVariable b) {
        ComplexVariable z;
        return z;
    }

    ComplexVariable operator*(int a, ComplexVariable b) {
      
        cout<<"****************"<<endl;
        if(b.maala==1){
            b.co2=1;
            b.co2*=a;
        }
        if(b.maala==2){
            b.co1*=a;
        }
        cout<<b.co1<<endl;
        cout<<b.co2<<endl;
        cout<<b.co3<<endl;
        cout<<b.maala<<endl;

        return b;
    }

    ComplexVariable operator^(ComplexVariable a, ComplexVariable b) {
        return a;
    }

    ComplexVariable operator/(int a, ComplexVariable b) {
        return a;
    }


    ComplexVariable operator+(double a, ComplexVariable b) {
        ComplexVariable z(b.co2 + a, b.co2);
        return z;
    };

    ComplexVariable operator-(double a, ComplexVariable b) {
        ComplexVariable z(b.co2 - a, b.im);
        return z;
    }

    ComplexVariable operator*(double a, ComplexVariable b) {
        ComplexVariable z(b.co2 * a, b.im * a);

        return z;
    }

    ComplexVariable operator^(ComplexVariable b,double a ) {
        cout<<"^^^^^^^^^^^^^^^^^^^^"<<endl;
        b.maala =a;
        if (a==2){
            b.co1=1;
        }
        cout<<b.co1<<endl;
        cout<<b.co2<<endl;
        cout<<b.co3<<endl;
        cout<<b.maala<<endl;

        return b;
    }

    ComplexVariable operator/(double a, ComplexVariable b) {
        if (b.co2 != 0) {
            ComplexVariable z((a * b.co2 + 0 * b.im) / (pow(b.co2, 2) + pow(b.im, 2)),
                              (0 * b.co2 - a * b.im) / (pow(b.co2, 2) + pow(b.im, 2)));
            return z;
        } else
            throw std::exception();
    }

    ostream &operator<<(ostream &os, const ComplexVariable &c) {
        return (os << c.co2 << '+' << c.im << 'i');
    }

    ostream &operator<<(ostream &os, const std::complex<double> &c) {
        return (os);
    }

    double solve(RealVariable a) {
        cout<<"########solve##########"<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        cout<<a.maala<<endl;

        if (a.maala ==2){
            cout<<"maala ==2"<<endl;
            double delta =a.co2*a.co2-(4*a.co1*a.co3);
            cout<<"delta "<<delta<<endl;
            try{
                if(delta<0||a.co1==0){
                    
                    throw "there is no solution";
                     
                    }
            }
            catch(const char* a){
                cout<<a<<endl;
            }
                

            
        double x = (-a.co2+sqrt(delta))/(2*a.co1);
            cout<<"answer "<<x<<endl;
            return x;
        }
        double delta =a.co2*a.co2*-a.co1*a.co3;
        if (a.co1==1){
            double answer = a.co3/(-a.co2);
            cout<<"answer "<<answer<<endl;
            return answer;
        }
       
       
   
        return 4;
    }

    std::complex<double> solve(ComplexVariable y) {
        cout<<"########solve##########"<<endl;
                    std::complex <double> x;

        double delta =y.co2*y.co2-4*y.co1*y.co3;
        cout<<"delta "<<delta<<endl;
        

        if (y.maala==1){
            double real = y.co3/(-y.co2);
            double img= -y.im/y.co2;
            std::complex <double> t(real,img);

            cout<<"real "<<real<<endl;
            cout<<"img "<<img<<endl;

            return t;
        
        }

        if(y.maala==2){
            if(delta>=0){
                double ans = (y.co2+sqrt(delta))/(2*y.co1);
                return x;
            }
            delta=-delta;
            double rr = sqrt(delta);

            double real= -y.co2/(2*y.co1);
            double img= rr/(2*y.co1);
            cout<<"real "<<real<<endl;
            cout<<"img "<<img<<endl;

             std::complex <double> v(real,img);


           return v;
        }
        cout<<"x "<<x<<endl;
        return x;
    }
}
int main(){
    using namespace std;
    using solver::solve; 
    using solver::RealVariable ;
    using solver::ComplexVariable;
    ComplexVariable x;
    solve((x^2)+(2*x)==-20);
    return 0;
}
