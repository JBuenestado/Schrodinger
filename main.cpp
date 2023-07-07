// Code by Jorge Buenestado, Schrödingers Problem
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>

using namespace std;

const int Num = 500;
const int ciclos = Num/25;
const int pasos = 500;
const double lambda = 0.3;
const double π = 3.14159265;
const complex<double> i (0,1);


int main(){

    int  j, l;
    double k, K0, Sm, norm, sig, x0;
    double V[Num];
    complex<double> alfa;
    complex<double> beta[Num], phi[Num+1], xhi[Num], A0[Num], b[Num];
    ofstream sch, norma;
    sch.open("DatosSchrodinger.txt");
    norma.open("norma.txt");

    cout<<"------------------Start of Program------------------"<<endl<<endl; 

    K0 = 2.0*π*ciclos/Num;
    Sm =  1.0/(4*K0*K0);
    alfa = 0;

    norm = 0;
    x0 = Num/4;
    sig = Num/16;

    //Potencial // A0
    for (j = 0; j < Num; j++){
        V[j] = 0.0;
        if((j>2*Num/5)&&(j<3*Num/5)){
            V[j] = lambda*pow(K0,2);
        }
        A0[j] = 2.0*i/Sm - V[j] - 2.0;
    }
    //Alfa
    for (j = 0; j < 10; j++){
        alfa = -1.0 / (A0[5] + alfa);
    }

    //Primer Phi
    phi[0] = 0.0;
    phi[Num] = 0.0;
    for (j = 1; j<Num-1; j++){
        k = j;
        phi[j] = exp(i*k*K0) * exp((1.0*pow( x0 - k ,2))/(2*sig*sig));
        norm += pow(abs(phi[j]),2);
    }

    for (j = 1; j<Num-1; j++){
        phi[j]=phi[j]/sqrt(norm);
    }

    //iteraciones
    for (j = 0; j < pasos ; j++){

        for (l = 0; l < Num ; l++){
            //sch << l << ", " << real(pow(abs(phi[l]),2)) << ", " << imag(pow(abs(phi[l]),2)) << ", " << pow(abs(phi[l]),2) << endl;
            sch << l << ", " << pow(abs(real(phi[l]/sqrt(norm))),2) << ", " << pow(abs(imag((phi[l]/sqrt(norm)))),2) << ", " << pow(abs(phi[l]/sqrt(norm)),2) << endl;
        }
        sch << endl;

        for (l = 0; l < Num ; l++){
            b[l] = 4.0*i*phi[l]/Sm;
        }
        
        beta[Num - 1] = 0.0;

        for (l = 0; l < Num - 1; l++){
            beta[Num - l - 2] = ( b[Num - l - 1]- beta[Num - l - 1])/(A0[Num - l - 1] + alfa);
        }

        xhi[0] = 0.0;
        xhi[Num - 1] = 0.0;

        for (l = 0; l < Num - 1; l++){
            xhi[l + 1] = alfa*xhi[l] + beta[l];
        }

        norm = 0;

        for (l = 0; l < Num; l++){
            phi[l] = xhi[l] - phi[l];
            norm += pow(abs(phi[l]),2);
        }
        norma << norm << endl;
    }
    sch.close();
    norma.close();

    cout<<"------------------End of Program------------------"<<endl<<endl;
    return 0;
}
