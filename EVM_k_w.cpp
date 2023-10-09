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

double  Y_NODE[NJ+1], DXF[4][NJ+1];//DXF node distance
double Y_plus[NJ+1];

double  S[4][NJ]; //surface area
double DELY[NJ];
double fy[NJ]; //non uniform mesh factors


//////////////////////////
/* turbulence data */
//////////////////////////

double AS[NJ+1], AN[NJ+1], AP[NJ+1];
double nu_t[NJ+1], nu_tf[NJ+1];
double SP[NJ+1]; 
double U_OLD[NJ+1],U[NJ+1];
double U_FACE[2][NJ+1];
double k[NJ+1], k_OLD[NJ+1];
double Pk[NJ+1];
double w[NJ+1], w_OLD[NJ+1];
double l_T[NJ+1];
double eps[NJ+1];
double uv[NJ+1];
double D_v[NJ+1];
double D_t[NJ+1];
double Dis[NJ+1];

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
    urf = 0. ;
    U_inf = 1;
    I = 0.05; //5 %

	nu=0.0025316456;
	
//	cout<<nu<<endl;
//	ustar=1;
	kappa=0.41;
	c_mu = 0.08;
	
	// k-w model constants
	beta = 0.09;
	c1_omega=0.55;
	c2_omega=0.075;
	sigma_k=2;
	sigma_omega=2;

	MAX_ERROR = 1e-04;	
	
    FILE_READ();	
	
	SET_GEOMETRY_NONUNIFORM();

    APPLY_BC();
    APPLY_IC();

//   INTERPOLATE();

    UPDATE();
    ofstream file3("ITER_VS_RESIDUE.dat");

    
ITER:

   ++ITER;
   
//   INTERPOLATE();
   
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

	for (j=1;j<=NJ;j++)	
	{
	cout<<Y_NODE[j]<<" "<<U[j]<<" "<<k[j]<<" "<<" "<<w[j]<<" "<<nu_t[j]<<" "<<Pk[j]<<endl;
	}

   cout<<"\n"<<"\n"<<endl;
      
   cout<<nu*((U[2] - U[1])/(Y_NODE[2] - Y_NODE[1]))<<endl;
  


}

void FILE_READ(){
	
	double xdata, ydata, udata, vdata;
	int k;
	k = 0;
	
	ifstream file1("y_dns.dat");

	while( file1 >> ydata)
	{
		Y1.push_back(ydata);
	}

	
	for(int i = 0; i <= Y1.size()-1;i++)
	{
//	cout <<Y1[i]<< endl;
	}

		
	file1.close();

		
}


void SET_GEOMETRY_NONUNIFORM()
{
	
	for(j=1;j<=NJ-1;j++)
	{
//	    X[j][i]= X1[i-1];
		Y_NODE[j]= Y1[j-1];
// 	  cout<<Y_NODE[j]<<" "<<" "<<Y[j]<<endl;
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
    
     cout<<"\n"<<endl;
       
	for(j=1;j<=NJ-2;j++)
	 {	
		
		DELY[j] = Y[j+1] - Y[j];
//		cout<<DELY[j]<<endl;
     }


	  	Y_NODE[1] = 0;
//	  	X_NODE[NJ][i] = X_NODE[NJ-1][i];
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

		U[NJ] = U[NJ-1]; //top //symmetric //zero gradient
		k[NJ] = k[NJ-1];
		w[NJ] = w[NJ-1];
		nu_t[NJ] = nu_t[NJ-1];
		nu_tf[NJ-1] = nu_t[NJ];
		
       
		
		U[1] =  0; //bottom
		k[1] = 0;
		nu_t[1] = 0;
		nu_tf[1] = 0;
		w[1] = 6*nu/(Y_NODE[2]*Y_NODE[2]*c2_omega);
        w[2] = 6*nu/(Y_NODE[2]*Y_NODE[2]*c2_omega);

//x-face velocity (in y direction)
	     	U_FACE[1][1]= U[1];
		 	U_FACE[1][NJ-1] = U[NJ];


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
		
		}


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

void CALC_EDDY_VIS(){
	
	for(j=2;j<=NJ-1;j++)
	{
/*
		
           if (ITER <= 100){
        
            l_T[j] = kappa*Y_NODE[j];
            eps[j] = c_mu*pow(k[j],1.5)/l_T[j];
           nu_t[j] = c_mu*k[j]*k[j]/eps[j]; //*******
       }
         else {
           nu_t[j] = k[j]/w[j];
         	
	   }
*/
           nu_t[j] = k[j]/w[j];

//            cout<<nu_t[j]<<endl;
    }
            nu_tf[1] = nu_t[1];
    		nu_tf[NJ-1] = nu_t[NJ];

    
	  for(j=2;j<=NJ-2;j++)
     {	    
       nu_tf[j] = fy[j]*nu_t[j] + (1-fy[j])*nu_t[j+1]; //*******
//       nu_tn[j] = fy[j]*nu_t[j] + (1-fy[j])*nu_t[j+1];//*******
      }
        
}


void CALC_U(){

	for(j=2;j<=NJ-1;j++)
	{

          SP[j] = 0;
		  		    
		  AS[j] = (nu + nu_tf[j-1]) * S[1][j] / (DXF[1][j]);
		  AN[j] = (nu+ nu_tf[j]) * S[3][j] / (DXF[3][j]); 
		  AP[j] = AN[j] + AS[j] - SP[j];
		  
//		  cout<<AN[j]<<" "<<AS[j]<<endl;
		}	
		
	for(j=2;j<=NJ-1;j++)
	{
		U[j] = (AN[j] * U[j+1] + AS[j] * U[j-1] + ((1/(rho))*1*DELY[j-1]))/ AP[j];			
//       cout<<AS[j]<<" "<<AN[j]<<" "<<AP[j]<<" "<<U[j]<<endl;
    }
    

	     	U_FACE[1][1]= U[1];
		 	U_FACE[1][NJ-1] = U[NJ];

	  for(j=2;j<=NJ-2;j++)
     {	
     U_FACE[1][j] = fy[j]*U[j] + (1-fy[j])*U[j+1];	  		
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
//		for(i=2;i<=NI-1;i++)
//		{

       
          SP[j] = -DELY[j-1]*beta*w[j];
		    
//		  AW[j][i] = kw[j][i] * S[0][j][i] / (DXF[0][j][i]);
		  AS[j] = (nu + (nu_tf[j-1])/sigma_k) * S[1][j] / (DXF[1][j]);
//		  AE[j][i] = ke[j][i] *S[2][j][i] / (DXF[2][j][i]);
		  AN[j] = (nu+ (nu_tf[j])/sigma_k) * S[3][j] / (DXF[3][j]);
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
		    
		  AS[j] = (nu + (nu_tf[j-1])/sigma_omega) * S[1][j] / (DXF[1][j]);

		  AN[j] = (nu+ (nu_tf[j])/sigma_omega) * S[3][j] / (DXF[3][j]);
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
//		  RESIDUE[j] =(U[j] - U_OLD[j])*(U[j] - U_OLD[j]);
          RESIDUE[j] =(U[j] - U_OLD[j]);
          RTOT=max(RESIDUE[j], RTOT);			
		}

    
	
//	RTOT = sqrt(RTOT/(NJ-2));
	return RTOT;
}



int FILE_WRITE1(){  // how to write the data in the file ...

	//grid line
	ofstream file5("turbulence_k_w_data.dat");
	ofstream file7("turbulence_k_w_data2.dat");


	for(j=2;j<=NJ-1;j++)
	{	      
       uv[j] = (-1)*nu_t[j]*((U_FACE[1][j] - U_FACE[1][j-1])/DELY[j-1]); 
	   D_v[j] = nu*((k[j+1]-2*k[j]+k[j-1])/(DELY[j-1]*DELY[j-1]));
	   D_t[j] = (nu_t[j]/sigma_k)*((k[j+1]-2*k[j]+k[j-1])/(DELY[j-1]*DELY[j-1]));
	   Dis[j] =  beta*w[j]*k[j];               
	}
      
      uv[1] = (-1)*nu_t[1]*((U[2] - U[1])/(Y_NODE[2] - Y_NODE[1])); //nu_t is zero at wall
      uv[NJ] = (-1)*nu_t[NJ]*((U[NJ]-U[NJ-1])/(Y_NODE[NJ] - Y_NODE[NJ-1]));// zero gredient
      Pk[NJ] =  nu_t[NJ]*pow(((U[NJ]-U[NJ-1])/(Y_NODE[NJ] - Y_NODE[NJ-1])),2); //zero gredient

      D_v[1] = nu*((k[1]- 2*k[2] + k[3])/pow((Y_NODE[2] - Y_NODE[1]),2));
      D_v[NJ] = nu*((k[NJ]-2*k[NJ-1] + k[NJ-2])/pow((Y_NODE[NJ] - Y_NODE[NJ-1]),2));
      D_t[1] = (nu_t[1]/sigma_k)*((k[1]- 2*k[2] + k[3])/pow((Y_NODE[2] - Y_NODE[1]),2));
      D_t[NJ] = (nu_t[NJ]/sigma_k)*((k[NJ]-2*k[NJ-1] + k[NJ-2])/pow((Y_NODE[NJ] - Y_NODE[NJ-1]),2));

      Dis[1] =  beta*w[1]*k[2];
      Dis[NJ] = beta*w[NJ]*k[NJ];

	for(j=1;j<=NJ;j++)	
	{
 	file5<<Y_NODE[j]<<" "<<U[j]<<" "<<k[j]<<" "<<w[j]<<" "<<nu_t[j]<<" "<<uv[j]<<" "<<Pk[j]<<endl;
	}

	file5.close();

	for(j=1;j<=NJ;j++)	
	{
 	file7<<Y_plus[j]<<" "<<Dis[j]<<" "<<Pk[j]<<" "<<D_t[j]<<" "<<D_v[j]<<endl;
	}
	
	file7.close();

}

