
#include <iostream>
#include <complex>
#include "solver.hpp"
#include <cmath>
#include <bits/exception.h>
namespace solver {

    RealVariable operator==( RealVariable a, const RealVariable b) {

        cout<<"=====aa=== r r"<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        cout<<a.maala<<endl;
        cout<<"=====b=== r r"<<endl;
        cout<<b.co1<<endl;
        cout<<b.co2<<endl;
        cout<<b.co3<<endl;
        cout<<b.maala<<endl;
        a.co1-=b.co1;
        a.co2-=b.co2;
        a.co3-=b.co3;

cout<<"=====after correction=== r r"<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        cout<<a.maala<<endl;
        return a;
    };
    RealVariable operator== (RealVariable a,  int b) {
        cout<<"============="<<endl;
        cout<<"before"<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        a.co3-=b;
         cout<<"after"<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        return a;
    }

    RealVariable operator+( RealVariable b,  RealVariable a) {
        cout<<"+++++++++++++++++a,second is "<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        cout<<endl;
        cout<<"b is"<<endl;
        cout<<b.co1<<endl;
        cout<<b.co2<<endl;
        cout<<b.co3<<endl;

        b.co1+=a.co1;
        b.co2+=a.co2;
        b.co3+=a.co3;
        
        b.maala = max(b.maala,a.maala);
        if(a.co1==0){a.maala=1;}
        if(a.co2==0){a.maala=0;}
        return b;
    }
    RealVariable operator+( RealVariable b,  int a) {
        cout<<"++++++++++++int!!"<<endl;
        cout<<"before"<<endl;
        cout<<b.co1<<endl;
        cout<<b.co2<<endl;
        cout<<b.co3<<endl;

        b.co3+=a;
        cout<<"after"<<endl;
        cout<<b.co1<<endl;
        cout<<b.co2<<endl;
        cout<<b.co3<<endl;

        return b;
    }
    RealVariable operator+( RealVariable b, double a) {
        b.co3+=a;
        return b;
    }
    
    
    RealVariable operator+( double a, RealVariable b) {
        b.co3+=a;
        return b;
    }
    RealVariable operator-( RealVariable a,const  RealVariable b) {
        cout<<"---------------- r r"<<endl;
        a.co1-=b.co1;
        a.co2-=b.co2;
        a.co3-=b.co3;
        if(a.co1==0){a.maala=1;}
        if(a.co2==0){a.maala=0;}

    return a;
    };
    RealVariable operator-( RealVariable a,  int b) {
        a.co3-=b;
        

    return a;
    };
        const RealVariable operator- (RealVariable const &a,  int const &b) {
        cout<<"--------------- r d"<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        cout<<a.maala<<endl;
        return a;
    }

    RealVariable operator*( RealVariable a,  RealVariable b) {
        
        return a;
    }
    const RealVariable operator/( RealVariable a, const RealVariable b) {
        cout<<"/////////bb?/////////////"<<endl;
        a.co2=a.co2/b.co1;
        return a;

    }
    const RealVariable operator/( RealVariable a,  int b) {
        cout<<"//////////////////////"<<endl;
        cout<<"before"<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        if(b==0){
            throw "wrong assignment";
        }
        if(a.co2!=0){
            a.co2=a.co2/b;
        }
         if(a.co1!=0){
            a.co1=a.co1/b;
        }
         cout<<"after"<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        return a;

    }

    RealVariable operator+(int b, RealVariable a) {
        cout<<"++++++++++++++++++"<<endl;
        cout<<"b is "<<b<<endl;
        cout<<"before"<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        a.co3+=b;
        cout<<"after"<<endl;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;

        return a;
    };

    RealVariable operator-(int b, RealVariable a) {
        
        return a;
    }


    RealVariable operator*(const int a, RealVariable b) {
        
        cout<<"***********"<<endl;
        if(b.maala==2){
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
       
        return a;

    }
    

    

    RealVariable operator-(double a, RealVariable b) {
        b.co3-=a;
        return b;
    }

    RealVariable operator*(double a, RealVariable b) {
        cout<<"****************"<<endl;
        cout<<"a is "<<a<<endl;
        cout<<"b.co1 is "<<b.co1<<endl;
        cout<<"b.co2 is "<<b.co2<<endl;
        cout<<"b.co3 is "<<b.co3<<endl;

        if(b.maala==1){
            b.co2=1;
            b.co2*=a;
        }
        if(b.maala==2){
            b.co2*=a;
        }
         cout<<"after"<<b.co1<<endl;

         cout<<"b.co1 is "<<b.co1<<endl;
        cout<<"b.co2 is "<<b.co2<<endl;
        cout<<"b.co3 is "<<b.co3<<endl;
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
        cout<<"a is  "<<a<<endl;
        if(a==0){
           throw "this is wrong";
        }
       
        if (a==2){
            b.co1=1;
        }
        cout<<b.co1<<endl;
        cout<<b.co2<<endl;
        cout<<b.co3<<endl;
        cout<<b.maala<<endl<<endl<<endl;

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

    double solve( RealVariable a) {
        if(a.maala==0){
            throw "the maala is 0";
        }
        RealVariable z(a.co1,a.co2,a.co3,a.maala);
       
        cout<<"########solve##########"<<endl;
        cout<<z.co1<<endl;
        cout<<z.co2<<endl;
        cout<<z.co3<<endl;
        cout<<z.maala<<endl;

        if (z.maala ==2){

            cout<<"maala ==2"<<endl;
            double delta =z.co2*z.co2-(4*z.co1*z.co3);
            cout<<"delta "<<delta<<endl;
            
                if(delta<0||z.co1==0){
                    
                    throw "there is no solution";
                     
                    }
            
            

            
        double x = (-z.co2+sqrt(delta))/(2*z.co1);
            cout<<"answer "<<x<<endl;
        a.co1=0;
        a.co2=1;
        a.co3=0;
        a.maala=1;
            return x;
        }
        if(a.co3==0||(a.co2==0&&a.co1==0)){
            throw "wrong assignment";
        }
        cout<<"#########aaa#########"<<endl;
        a.co1=0;
        a.co2=1;
        a.co3=0;
        a.maala=1;
        cout<<a.co1<<endl;
        cout<<a.co2<<endl;
        cout<<a.co3<<endl;
        cout<<a.maala<<endl;
        
        double answer = z.co3/(-z.co2);
        cout<<"answer "<<answer<<endl;
        
        return answer;
        
       
       
   
        
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
// int main(){
//     using namespace std;
//     using solver::solve; 
//     using solver::RealVariable ;
//     using solver::ComplexVariable;
//     ComplexVariable x;
//     solve((x^2)+(2*x)==-20);
//     return 0;
// }
