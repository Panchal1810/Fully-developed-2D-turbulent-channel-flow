#include<iostream>
#include<stdio.h>
#include<math.h>
#include <cmath>
#include<string.h>
#include<iomanip>
#include<fstream>
#include<string>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <vector>


//#define NI 97 //nodes in x direction
#define NJ 97 //nodes in y direction

using namespace std;


/////////////////////
/* geomeUric daUa */
////////////////////

float L, delta, H; 
int i,j;
//double hA, hC;

//vector<double> X1;
vector<double> Y1;

//double X[NJ][NI], Y[NJ][NI];
double Y[NJ];

//double  X_NODE[NJ+1][NI+1], Y_NODE[NJ+1][NI+1], DXF[4][NJ+1][NI+1];//DXF node distance
double  Y_NODE[NJ+1], DXF[4][NJ+1];//DXF node distance
double Y_plus[NJ+1];


double  S[4][NJ]; //surface area
double DELY[NJ];
double fy[NJ]; //non uniform mesh factors


//////////////////////////
/* turbulence data */
//////////////////////////


double AS[NJ+1], AN[NJ+1], AP[NJ+1];
double nu_t[NJ+1], nu_ts[NJ+1], nu_tn[NJ+1];
double SP[NJ+1]; 
double U_OLD[NJ+1],U[NJ+1];
double U_FACE[2][NJ+1];
double k[NJ+1], k_OLD[NJ+1];
double Pk[NJ+1];
double w[NJ+1], w_OLD[NJ+1];
double l_T[NJ+1];
double eps[NJ+1];
double uv[NJ+1];

double rho;  //properties
double Su; //source term zero here
double I, U_inf;
double nu, ustar, kappa, beta, c1_omega, c2_omega, sigma_k, sigma_omega;
double c_mu;
//double -dp/dx ;


///////////////////////////
/* convergence and error */
///////////////////////////

double RESIDUE[NJ+1];
double urf;
double MAX_ERROR, ABS_ERROR, RMS_ERROR, RRESIDUE;
int ITER = 0;

void SET_GEOMETRY_UNIFORM();
void SET_GEOMETRY_NONUNIFORM(); //here non uniform mesh is given
void APPLY_IC();
void APPLY_BC();
void UPDATE();
void CALC_U();
void CALC_EDDY_VIS();
void CALC_Pk();
void CALC_k();
void CALC_w();


double CALC_U_ERROR();



int FILE_WRITE1(); //writing post processing file in .dat format 
//int FILE_WRITE2();
void FILE_READ(); //reading geometric and u and v data
void INTERPOLATE(); // interpolating u and v data from node to face center (convective part)

double A;

int main ()
{
	delta = 1;
	L = 2; H = 1;//2*delta;
//	hA  = 0.068*H;
//	hC  = 0.068*H;
	rho  = 1;
    urf = 0.8 ;
    U_inf = 1;
    I = 0.05; //5 %

	nu=0.0025316456;
	
//	cout<<nu<<endl;
	ustar=1;
	kappa=0.41;
	c_mu = 0.08;
	
	// k-w model constants
	beta = 0.09;
	c1_omega=0.55;
	c2_omega=0.075;
	sigma_k=2;
	sigma_omega=2;

//    -dp/dx = 1;
	
//	ST_FACTOR_Y = 1.1;
//	ST_FACTOR_X = 1.08;
	MAX_ERROR = 1e-04;	
	
    FILE_READ();	
	
	SET_GEOMETRY_NONUNIFORM();
	
//	INTERPOLATE();

    APPLY_BC();
    APPLY_IC();


    UPDATE();
    ofstream file3("ITER_VS_RESIDUE.dat");

    
ITER:

   ++ITER;
   
   INTERPOLATE();
   
   CALC_EDDY_VIS();
   
       
   CALC_U();
   
   CALC_Pk();
   
   CALC_k();
   
    CALC_w();
    
	for (j=1;j<=NJ;j++)	
	{
//	cout<<Y_NODE[j]<<" "<<U[j]<<" "<<k[j]<<" "<<" "<<w[j]<<" "<<nu_t[j]<<" "<<Pk[j]<<endl;
	}

   
    APPLY_BC();
    
    RRESIDUE = CALC_U_ERROR();

    
    file3<<ITER<<"\t"<<RRESIDUE<<"\n";
//    cout<<ITER<<" "<<RRESIDUE<<" "<<endl;

	if(RRESIDUE > MAX_ERROR){
		UPDATE();
		goto ITER;
	}

    for(j=1;j<=NJ;j++){
//    	cout<<U[j]<<" "<<endl;
	}
    
    
    file3.close();

	FILE_WRITE1();
     


}

void FILE_READ(){
	
	double xdata, ydata, udata, vdata;
	int k;
	k = 0;
	

	while( file1 >> ydata)
	{
		Y1.push_back(ydata);
	}

	
	for(int i = 0; i <= Y1.size()-1;i++)
	{

	}

		
	file1.close();

		
}


void SET_GEOMETRY_NONUNIFORM()
{
	

	for(j=1;j<=NJ-1;j++)
	{

		Y_NODE[j]= Y1[j-1];

    } 
    
        Y_NODE[1] = 0;
	  	Y_NODE[NJ] = H;
    
    for(j=1;j<=NJ-1;j++){
    	if(j==1){
			Y[j] = Y_NODE[1];
		}
		
		else if(j==NJ-1){
			Y[j] = Y_NODE[NJ];
		}
		
		else{
			Y[j] = 0.5*(Y_NODE[j] + Y_NODE[j+1]);
		}
//            cout<<Y[j]<<" "<<" "<<Y_NODE[j]<<endl;

	}
       

	  	
    for(j=1;j<=NJ;j++){
    	Y_plus[j] = (Y_NODE[j])/nu;
 //   	cout<<Y_NODE[j]<<" "<<" "<<Y_plus[j]<<endl;
	}
    
 
	  	Y_NODE[1] = 0;

	  	Y_NODE[NJ] = H;
  
//NODE DISTANCE

for(j=2;j<=NJ-1;j++)
	{

		DXF[1][j] = (Y_NODE[j]-Y_NODE[j-1]); //SOUTH

		DXF[3][j] = (Y_NODE[j+1]-Y_NODE[j]); //NORTH

	}    
		
//SURFACE AREA

	for(j=2;j<=NJ-1;j++)
	{

		S[1][j] = 1;//(X[j-1][i] - X[j-1][i-1]);	//SOUTH FACE 

	    S[3][j] = 1;//(X[j][i] - X[j][i-1]);    //NORTH FACE
        
    }
	

	for(j=2;j<=NJ-1;j++)
	{
//     fx[j][i] = 0.5*DELX[j-1][i-1]/(X_NODE[j][i+1]-X_NODE[j][i]);
	   fy[j] = 0.5*DELY[j-1]/(Y_NODE[j+1]-Y_NODE[j]);
    }
 				

}  


void INTERPOLATE(){


//x-face velocity (in y direction)
	  for(j=1;j<=NJ-1;j++)
	  {	
	     if (j==1){
	     	U_FACE[1][1]= U[1];
		 }
		 else if (j==NJ-1){
		 	U_FACE[1][NJ-1] = U[NJ];
		 }
		 else{
         U_FACE[1][j] = fy[j]*U[j] + (1-fy[j])*U[j+1];	  		
	     }
	//     cout<<U_FACE[1][j][i]<<endl;
	  }
}




void APPLY_BC(){


//	for(i=1;i<=NI;i++) // top and bottom boundary condition
//	{
		U[NJ] = U[NJ-1]; //top //symmetric //zero gradient
		k[NJ] = k[NJ-1];
		w[NJ] = w[NJ-1];
		nu_t[NJ] = nu_t[NJ-1];
       
		
		U[1] =  0; //bottom
		k[1] = 0;
		nu_t[1] = 0;
		w[1] = 6*nu/(Y_NODE[2]*Y_NODE[2]*c2_omega);
        w[2] = 6*nu/(Y_NODE[2]*Y_NODE[2]*c2_omega);
//         cout<<w[1]<<endl;
//	}


}

void APPLY_IC(){
	
		for(j=2;j<=NJ-1;j++) // interior nodes
	   {	
		 U[j] = 1; 
		 k[j] = 1.5*(I*U_inf)*(I*U_inf);
		 l_T[j] = kappa*Y_NODE[j];//min((pow(c_mu,0.25)*pow(k[j][i],1.5)/eps),(kappa*S/U_secd))
        
        if(Y_plus[j] <=3){
		w[j] = 6*nu/(Y_NODE[j]*Y_NODE[j]*c2_omega);
        	
		}

        else {
		 w[j] = (pow(c_mu,0.25)*pow(k[j], 0.5))/ (beta*l_T[j]);
		}
		
//		cout<<U[j]<<" "<<k[j]<<" "<<w[j]<<endl;
}
}

void CALC_EDDY_VIS(){
	
	for(j=2;j<=NJ-1;j++)
	{

		
           if (ITER <= 100){
        
            l_T[j] = kappa*Y_NODE[j];
            eps[j] = c_mu*pow(k[j],1.5)/l_T[j];
           nu_t[j] = c_mu*k[j]*k[j]/eps[j]; //*******
       }
         else {
           nu_t[j] = k[j]/w[j];
         	
	   }

    }
        
}


void CALC_U(){

	for(j=2;j<=NJ-1;j++)
	{

          SP[j] = 0;
		  		    

		  AS[j] = (nu + nu_t[j-1]) * S[1][j] / (DXF[1][j]);

		  AN[j] = (nu+ nu_t[j]) * S[3][j] / (DXF[3][j]); 
		  AP[j] = AN[j] + AS[j] - SP[j];
		  

		}	
		
	for(j=2;j<=NJ-1;j++)
	{
		U[j] = (AN[j] * U[j+1] + AS[j] * U[j-1] + ((1/(rho))*1*DELY[j-1]))/ AP[j];			

	}

}

void CALC_Pk(){
	
	for(j=2;j<=NJ-1;j++)
	{

          Pk[j] = nu_t[j]*pow(((U_FACE[1][j] - U_FACE[1][j-1])/DELY[j-1]),2);
                    
	  }
   
}

void CALC_k(){

	for(j=2;j<=NJ-1;j++)
	{

       
          SP[j] = -DELY[j-1]*beta*w[j];

		  AS[j] = (nu + (nu_t[j-1])/sigma_k) * S[1][j] / (DXF[1][j]);

		  AN[j] = (nu+ (nu_t[j])/sigma_k) * S[3][j] / (DXF[3][j]);
		  AP[j] = AN[j] + AS[j] - SP[j] ;

	}
		
	for(j=2;j<=NJ-1;j++)
	{
		k[j] = (AN[j] * k[j+1] + AS[j] * k[j-1] + DELY[j-1]*Pk[j])/ AP[j];			
	}
		
}

void CALC_w(){


	for(j=3;j<=NJ-1;j++)
	{

          SP[j] =  DELY[j-1]*(((c1_omega*Pk[j])/k[j]) - (c2_omega*w_OLD[j]));

		  AS[j] = (nu + (nu_t[j-1])/sigma_omega) * S[1][j] / (DXF[1][j]);

		  AN[j] = (nu+ (nu_t[j])/sigma_omega) * S[3][j] / (DXF[3][j]);
		  AP[j] = AN[j] + AS[j] - SP[j] ;

	     }
	 
	for(j=3;j<=NJ-1;j++)
	{

		w[j] = (AN[j] * w[j+1] + AS[j] * w[j-1])/ AP[j];			
	   }


}



void UPDATE(){

	for (j=2;j<=NJ-1;j++)
	 {	
		 	 
	 U_OLD[j] =  urf*U_OLD[j] + (1-urf)*U[j];
	 k_OLD[j] =  urf*k_OLD[j] + (1-urf)*k[j];
	 w_OLD[j] =  urf*w_OLD[j] + (1-urf)*w[j];
   
	 }
}
  


double CALC_U_ERROR()
{	
	double RTOT = 0;
	for (j=2;j<=NJ-1;j++)
	{
		  RESIDUE[j] =(U[j] - U_OLD[j])*(U[j] - U_OLD[j]);
          RTOT=RTOT+RESIDUE[j];			
		}

	
	RTOT = sqrt(RTOT/(NJ-2));
	return RTOT;
}



int FILE_WRITE1(){  // how to write the data in the file ...

	//grid line
	ofstream file5("turbulence_k_w_data.dat");


	for(j=2;j<=NJ-1;j++)
	{	      
       uv[j] = (-1)*nu_t[j]*((U_FACE[1][j] - U_FACE[1][j-1])/DELY[j-1]); 
	   D_v[j] = nu*((k[j+1]-2*k[j]+k[j-1])/(DELY[j-1]*DELY[j-1]));
	   D_t[j] = (nu_t[j]/sigma_k)*((k[j+1]-2*k[j]+k[j-1])/(DELY[j-1]*DELY[j-1]));
	   Dis[j] =  beta*w[j]*k[j];               
	}
      
      uv[1] = 0;
      uv[NJ] = 0;
      
      
      

//	for (i=1;i<=NI-1;i++)
//	{
		for (j=1;j<=NJ;j++)	
		{
	 	file5<<Y_NODE[j]<<" "<<U[j]<<" "<<k[j]<<" "<<w[j]<<" "<<nu_t[j]<<" "<<uv[j]<<" "<<Pk[j]<<endl;
		}
//	}
	file5.close();


	ofstream file7("turbulence_k_w_data2.dat");
		for (j=1;j<=NJ;j++)	
		{
	 	file5<<Y_plus[j]<<" "<<Dis[j]<<" "<<Pk[j]<<" "<<D_t[j]<<" "<<D_v[j]<<endl;
		}
	
	file7.close();


}


