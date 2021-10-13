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
#define NX 600
#define NY 150
#define NX1 (NX+1)
#define NY1 (NY+1)

#define Q 9

#define Base_T 10
#define Base_B 10
#define Liquid_H 90
#define Vapor_H (NY - Base_T - Base_B - Liquid_H)

#define dt 1.0
#define cc 1.0
#define c_squ (1.0/3.0)
#define rho_s 1.0
#define rho_v 0.3797
#define rho_l 6.4989
#define rho_a ((rho_v+rho_l)/2)

#define nu_l 0.1
#define nu_v (0.5/3.0)

#define tau_l (3 * nu_l + 0.5)
#define tau_v (3 * nu_v + 0.5)
#define tau_e (1.0/1.1)
#define tau_t (1.0/1.1)



#define a (3.0/49.0)
#define b (2.0/21.0)
#define R 1.0
#define ome 0.344

#define dT  0.0137
#define T_c (0.0778/0.45724*a/(b*R))
#define T_s (0.86*T_c)
#define T_b (T_s + dT)

#define g (3e-5)



#pragma endregion

#pragma region variable
int i,j,k,m,n,NMAX,Nplus;

double rho[NX1][NY1];
double P[NX1][NY1];

double f[NX1][NY1][Q],F[NX1][NY1][Q], mf[NX1][NY1][Q], MF[NX1][NY1][Q];

double force_b[NX1][NY1][2];
double force_m[NX1][NY1][2];
double force_s[NX1][NY1][2];
double force[NX1][NY1][2];

double u[NX1][NY1][2],u0[NX1][NY1][2];

double fai[NX1][NY1];





int area[NX1][NY1];

double T[NX1][NY1];

int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};

double s[Q];
double M0[Q][Q] = {{ 1.0,   1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0},
                   {-4.0,  -1.0, -1.0, -1.0, -1.0,  2.0,  2.0,  2.0,  2.0},
                   { 4.0,  -2.0, -2.0, -2.0, -2.0,  1.0,  1.0,  1.0,  1.0},
                   { 0.0,   1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,  1.0},
                   { 0.0,  -2.0,  0.0,  2.0,  0.0,  1.0, -1.0, -1.0,  1.0},
                   { 0.0,   0.0,  1.0,  0.0, -1.0,  1.0,  1.0, -1.0, -1.0},
                   { 0.0,   0.0, -2.0,  0.0,  2.0,  1.0,  1.0, -1.0, -1.0},
                   { 0.0,   1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0},
                   { 0.0,   0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0}};
double Mre[Q][Q] = { {1.0/9.0, -1.0/9.0,   1.0/9.0,     0,        0,          0,       0,           0,        0},
                     {1.0/9.0, -1.0/36.0, -1.0/18.0, 1.0/6.0,  -1.0/6.0,      0,       0,        1.0/4.0,     0},
                     {1.0/9.0, -1.0/36.0, -1.0/18.0,    0,        0,       1.0/6.0,  -1.0/6.0,  -1.0/4.0,     0},
                     {1.0/9.0, -1.0/36.0, -1.0/18.0, -1.0/6.0,  1.0/6.0,      0,       0,        1.0/4.0,     0},
                     {1.0/9.0, -1.0/36.0, -1.0/18.0,    0,        0,      -1.0/6.0,   1.0/6.0,  -1.0/4.0,     0},
                     {1.0/9.0,  1.0/18.0,  1.0/36.0,  1.0/6.0,  1.0/12.0,  1.0/6.0,   1.0/12.0,     0,     1.0/4.0},
                     {1.0/9.0,  1.0/18.0,  1.0/36.0, -1.0/6.0, -1.0/12.0,  1.0/6.0,   1.0/12.0,     0,    -1.0/4.0},
                     {1.0/9.0,  1.0/18.0,  1.0/36.0, -1.0/6.0, -1.0/12.0, -1.0/6.0,  -1.0/12.0,     0,     1.0/4.0},
                     {1.0/9.0,  1.0/18.0,  1.0/36.0,  1.0/6.0,  1.0/12.0, -1.0/6.0,  -1.0/12.0,     0,    -1.0/4.0}};
double w[Q]={ 4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};

double w2[Q]={0.0, 1.0/3, 1.0/3, 1.0/3, 1.0/3, 1.0/12, 1.0/12, 1.0/12, 1.0/12};



#pragma endregion

#pragma region function

void init();
void model();
void streaming();
void boundary();
void macro();
void collision();
double Peos(double rho, double T);


void RK();
void output();

double feq(int k, double rho, double u[2]);

#pragma endregion

#pragma region main loop


int main(){
    init();
    model();

    cout << "The max number of evolution: ";
    cin >> NMAX;

    for (n = 0; n < NMAX; n++){
        streaming();
        boundary();
        macro();

        collision();
        output();
    }

    return 0;
}


#pragma endregion

#pragma region function definition

void model(){
    for(i = 0; i <= NX; i ++)
        for(j = 0; j <= NY; j ++)
            {
                //solid
                if(j <= Base_B || j>= NY - Base_T) //TODO NEED TO CHEEK
                    area[i][j] = 1;
                //liquid
                if(j > Base_B  && j<= Liquid_H +Base_B  )
                    area[i][j] = 2;
                //vapor
                else
                    area[i][j] = 3;
            }
}

void init(){

//密度初始化
	for (i=0;i<=NX;i++)
	{
		for (j=0;j<=NY;j++)
		{
			if (area[i][j]==0)
                rho[i][j]= rho_s;
			else if (area[i][j]==1)
			{
				rho[i][j]=rho_l;
			}else if(area[i][j]==2){

                rho[i][j]=rho_v;
			}
		}
	}
//温度、压力初始化
for (i=0;i<=NX;i++)
	{
		for (j=0;j<=NY;j++)
		{
			if (area[i][j]==0)//固体
                T[i][j] = T_b;
			else if (area[i][j])//液体 & 气体
			{
				T[i][j] = T_s;
				P[i][j] = Peos(rho[i][j], T[i][j]);

			}
		}
	}
//分布函数、作用力初始化

    for (i=0;i<=NX;i++)
	{
		for (j=0;j<=NY;j++)
		{
			for (k=0;k<2;k++)
			{
				u[i][j][k]=0;
				u0[i][j][k]=0;
				force_b[i][j][k]=0;
				force_m[i][j][k]=0;
				force_s[i][j][k]=0;

			}
			for (k=0;k<Q;k++)
			{
				f[i][j][k]=feq(k,rho[i][j],u[i][j]);
				F[i][j][k]=f[i][j][k];
			}
		}
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

void RK(){

}

double Peos(double rho, double T){
    double p;
    double phi;

    phi = (1+(0.37464 + 1.54226*ome - 0.26992*ome*ome)*(1 - sqrt(T_s/T_c)));
    p = (rho * R * T)/(1-b*rho) - (a*phi)/(1+2*b*rho-b*b*rho*rho);

    return p;
}

double feq(int k, double rho, double u[2])
{
	double eu, uv, feq;
	eu= (e[k][0]*u[0]+e[k][1]*u[1]);
	uv= (u[0]*u[0]+u[1]*u[1]);
	feq=w[k]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
	return feq;
}

void output(){
    ostringstream name;
    const int aa[] = {10,100,1000,10000,100000,1000000,10000000,100000000};
    int bb[8];
    for(i = 0; i < 8; i ++)
    {
        bb[i] = n % aa[i] / (aa[i] / 10);
    }
    name << "tecplot"
    << bb[7] << bb[6] << bb[5] << bb[4] << bb[3] << bb[2] << bb[1] << bb[0]
    << ".dat";
    ofstream out (name.str().c_str());
    out << "TITLE = \"Tecplot Data\"\n"
    << "VARIABLES = \"X\", \"Y\", \"rho\", \"u\", \"v\", \"area\", \"T\", \"p\" \n"
    << "ZONE T= \"Zone 1\", J=" << NX1 << ", I=" << NY1-2 << ", F=POINT\n";
    for(i = 0; i < NX1; i ++)
        for(j = 1; j < NY1-1; j ++)
        {
            out << i << "\t"
            << j << "\t"
            << rho[i][j] << "\t"
            << u[i][j][0] << "\t"
            << u[i][j][1] << "\t"
            << area[i][j] << "\t"
            << T[i][j] << "\t"
            << P[i][j] << "\t"
            << endl;
        }
}
#pragma endregion