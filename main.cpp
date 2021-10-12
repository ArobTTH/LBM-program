#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <string>
#include <ctime>

using namespace std;

#pragma region constant
#define MX 600
#define MY 150
#define MX1 (MX+1)
#define MY1 (MY+1)
#define MX2 (MX+2)
#define MY2 (MY+2)

#define Q 9

#define Base_T 10
#define Base_B 10


#pragma endregion

#pragma region variable
int i,j,m,n,NMAX,Nplus;

double rho[MX1][MY1];
double f[MX1][MY1][Q],fcol[MX1][MY1][Q];
double u[MX1][MY1],v[MX1][MY1];
double du[MX1][MY1],dv[MX1][MY1];
double u_real[MX1][MY1],v_real[MX1][MY1];
double u_pre[MX1][MY1],v_pre[MX1][MY1];
double fai[MX1][MY1];
double p[MX1][MY1];


int SorL[MX1][MY1];



#pragma endregion

#pragma region function

void init();
void model();
void streaming();
void boundary();
void macro();
void collision();
void force();
void RK();
void output();


#pragma endregion

#pragma region main loop


int main(){
    init();
    model();
    init_psy();

    cout << "The max number of evolution: ";
    cin >> NMAX;

    for (n = 0; n < NMAX; n++){
        streaming();
        boundary();
        macro();
        force();
        collision();
        output();
    }

    return 0;
}


#pragma endregion

#pragma region function definition
void init(){


}
void model(){
    for(i = 0; i < MX2; i ++)
        for(j = 0; j < MY2; j ++)
            {
                if(j <Base_B || j> MY1 - Base_T)
                    SorL[i][j] = 1;
                else
                    SorL[i][j] = 0;
            }
}


void streaming(){

}
void boundary(){

}
void macro(){

}
void collision(){

}
void force(){

}
void RK(){

}
void output(){
    ostringstream name;
    const int a[] = {10,100,1000,10000,100000,1000000,10000000,100000000};
    int b[8];
    for(i = 0; i < 8; i ++)
    {
        b[i] = n % a[i] / (a[i] / 10);
    }
    name << "tecplot"
    << b[7] << b[6] << b[5] << b[4] << b[3] << b[2] << b[1] << b[0]
    << ".dat";
    ofstream out (name.str().c_str());
    out << "TITLE = \"Tecplot Data\"\n"
    << "VARIABLES = \"X\", \"Y\", \"rho\", \"u\", \"v\", \"SorL\", \"T\", \"p\" \n"
    << "ZONE T= \"Zone 1\", J=" << MX1 << ", I=" << MY1-2 << ", F=POINT\n";
    for(i = 0; i < MX1; i ++)
        for(j = 1; j < MY1-1; j ++)
        {
            out << i << "\t"
            << j << "\t"
            << rho[i][j] << "\t"
            << u_real[i][j] << "\t"
            << v_real[i][j] << "\t"
            << SorL[i][j] << "\t"
            << T[i][j] << "\t"
            << p[i][j] << "\t"
            << endl;
        }
}
#pragma endregion