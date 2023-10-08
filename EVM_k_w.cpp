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


//float CELL_ZI[NI+1][NJ+1],CELL_ZEUA[NI+1][NJ+1],

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
/*
double  S[4][NJ][NI]; //surface area
double  VOL[NJ][NI];  //volume 
double DELX[NJ][NI], DELY[NJ][NI];
double fx[NJ][NI],  fy[NJ][NI]; //non uniform mesh factors
*/

double  S[4][NJ]; //surface area
double DELY[NJ];
double fy[NJ]; //non uniform mesh factors


//////////////////////////
/* turbulence data */
//////////////////////////


//double U[2][NJ+1][NI+1];//U-field given here U[0][j][i] means u and U[1][j][i] means v
//double U_FACE[2][NJ+1][NI+1];//U-field given face velocity

/*
double AW[NJ+1][NI+1], AS[NJ+1][NI+1], AE[NJ+1][NI+1], AN[NJ+1][NI+1], AP[NJ+1][NI+1];
double nu_t[NI+1][NJ+1], nu_ts[NI+1][NJ+1], nu_tn[NI+1][NJ+1];
double SP[NJ+1][NI+1]; 
double U_OLD[NJ+1][NJ+1],U[NJ+1][NI+1];
double U_FACE[2][NJ+1][NI+1];
double k[NJ+1][NI+1], k_OLD[NJ+1][NI+1];
double Pk[NJ+1][NI+1];
double w[NJ+1][NI+1], w_OLD[NJ+1][NI+1];
double l_T[NJ+1][NI+1];
double eps[NJ+1][NI+1];
*/

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

//void CALC_CONV(); //equation solved using gauss seidal
//double CALC_ABS_ERROR();
//double CALC_CONV_ERROR();
double CALC_U_ERROR();

//TDMA algorithm
//double a[NJ+1][NI+1], b[NJ+1][NI+1], c[NJ+1][NI+1], d[NJ+1][NI+1];
//double P[NJ+1][NI+1], Q[NJ+1][NI+1];
//void CALC_CONV_VTDMA(); //equation solved using vertical TDMA
//void CALC_CONV_HTDMA(); //equation solved using horizontal TDMA


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

/*

    for(j=1;j<=NJ;j++){
    	cout<<U[j]<<" "<<k[j]<<" "<<eps[j]<<endl;
	}
*/   

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
	cout<<Y_NODE[j]<<" "<<U[j]<<" "<<k[j]<<" "<<" "<<w[j]<<" "<<nu_t[j]<<" "<<Pk[j]<<endl;
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
	
//	ifstream file1("xc.dat");
	ifstream file1("y_dns.dat");
//	ifstream file3("u.dat");
//	ifstream file4("v.dat");
	
/*
	while ( file1 >> xdata)
	{
		X1.push_back(xdata);
	}
*/	
	while( file1 >> ydata)
	{
		Y1.push_back(ydata);
	}

/*	
	while ( file3 >> udata)
	{
		u.push_back(udata);
	}
	
	while ( file4 >> vdata)
	{
		v.push_back(vdata);
	}

*/
	
	for(int i = 0; i <= Y1.size()-1;i++)
	{
//	cout <<Y1[i]<< endl;
	}

		
	file1.close();

		
}


void SET_GEOMETRY_NONUNIFORM()
{
	
/*	                     // in compuational domain generates 1x1 squre domain.
	delx = (L/(NI-2));    // X is horizontal direction in physical domain
	dely = (H/(NJ-2));  // Y is vertical direction in physical domain
	delv = delx * dely ;
//	cout<<delv<<endl ;
//	cout<<delx<<" "<<dely<<endl;
*/
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


    for(j=1;j<=NJ-1;j++){
//    	cout<<Y[j]<<" "<<DELY[j]<<endl;
	}

/*	 	
	  for(j=2;j<=NJ-1;j++)
	  {
//		X_NODE[j][i] = (X[j][i]+X[j][i-1])/2;
		Y_NODE[j] = (Y[j]+Y[j-1])/2;       
//        cout<<Y_NODE[j]<<endl;
	  } 
*/

/*
for(j=2;j<=NJ-1;j++) // left and right boundary
	  {
//	  	X_NODE[j][1] = 0;
	  	Y_NODE[j] = Y_NODE[j];
//	  	X_NODE[j][NI] = L;
	  	Y_NODE[j] =Y_NODE[j];
	  }
*/

//for(i=2;i<=NI-1;i++) //top and bottom boundary
//	  {
//	  	X_NODE[1][i] = X_NODE[2][i];
	  	Y_NODE[1] = 0;
//	  	X_NODE[NJ][i] = X_NODE[NJ-1][i];
	  	Y_NODE[NJ] = H;
	  	
//	  }
	  
/*	  
	  Y_NODE[1][1] = 0;  //four corner of mesh
	  Y_NODE[1][NI]=0;
	  Y_NODE[NJ][1] = H;
	  Y_NODE[NJ][NI] = H;
*/
	  
	  
//	  X_NODE[1][1] = 0; 
//	  X_NODE[1][NI] = L; 
//	  X_NODE[NJ][1] = 0;  
//	  X_NODE[NJ][NI] = L; 

//face center
/*
for(i=1;i<=NI-1;i++)
	{
	for(j=2;j<=NJ-1;j++)
	{
	//	Y_FACE_CENTER_X[j][i]= X[j][i];
	//	cout<<Y_FACE_CENTER_X[j][i]<<endl;
	    } 
	}	
*/	
	  
//NODE DISTANCE

for(j=2;j<=NJ-1;j++)
	{
//		DXF[0][j][i] = (X_NODE[j][i]-X_NODE[j][i-1]); //WEST
		DXF[1][j] = (Y_NODE[j]-Y_NODE[j-1]); //SOUTH
//		DXF[2][j][i] = (X_NODE[j][i+1]-X_NODE[j][i]); //EAST
		DXF[3][j] = (Y_NODE[j+1]-Y_NODE[j]); //NORTH
	  //  cout<<DXF[0][j][i]<<" "<<DXF[1][j][i]<<" "<<DXF[2][j][i]<<" "<<DXF[3][j][i]<<endl;
	  //  cout<<"\n"<<endl;
	}    
		
//SURFACE AREA

	for(j=2;j<=NJ-1;j++)
	{
//		S[0][j][i] = (Y[j][i-1] - Y[j-1][i-1]);//WEST FACE
		S[1][j] = 1;//(X[j-1][i] - X[j-1][i-1]);	//SOUTH FACE 
//	    S[2][j][i] = (Y[j][i] - Y[j-1][i]); //EAST FACE
	    S[3][j] = 1;//(X[j][i] - X[j][i-1]);    //NORTH FACE
	    
//	    VOL[j][i] = S[2][j][i] * S[3][j][i] ;
	//    cout<<S[2][j][i]<<" "<<S[3][j][i]<<endl;         
    }
	

	for(j=2;j<=NJ-1;j++)
	{
//     fx[j][i] = 0.5*DELX[j-1][i-1]/(X_NODE[j][i+1]-X_NODE[j][i]);
	   fy[j] = 0.5*DELY[j-1]/(Y_NODE[j+1]-Y_NODE[j]);
    }
 				

}  


void INTERPOLATE(){
/*	
 //y-face velocity (in x direction)
    for(i=1;i<=NI-1;i++)
	{
	  for(j=2;j<=NJ-1;j++)
	  {	
	     if (i==1)
	     {
	     	U_FACE[0][j][1] = U[0][j][1] ;
		 }
		 else if (i==NI-1){
		 	U_FACE[0][j][NI-1] = U[0][j][NI] ;
		 }
		 else{
         U_FACE[0][j][i] = fx[j][i]*U[0][j][i] + (1-fx[j][i])*U[0][j][i+1];	  		
	     }
	  //   cout<<U_FACE[0][j][i]<<endl;
	  }
    }
*/

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
/*	
	 for(j=1;j<=NJ;j++) // left and Right boundary condition
	{
		 U[j] = 1; //left
		 k[j] = 1.5*(I*U_inf)*(I*U_inf);
		 l_T[j] = kappa*Y_NODE[j];//min((pow(c_mu,0.25)*pow(k[j][i],1.5)/eps),(kappa*S/U_secd))
		 w[j] = (pow(c_mu,0.25)*pow(k[j], 0.5))/ (beta*l_T[j]);
		 
		 
		U[j][NI] = U[j][NI-1]; //Right //zero gradient
		k[j][NI] = k[j][NI-1]; //Right		
		w[j][NI] = w[j][NI-1];

	}

*/

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

//           nu_t[j] = k[j]/w[j];

//           nu_ts[j][i] = c_mu*k[j-1][i]*k[j-1][i]/eps[j][i]; //*******
//           nu_tn[j][i] = c_mu*k[j][i]*k[j-1][i]/eps[j][i]; //*******
//            cout<<nu_t[j]<<endl;
    }
        
}


void CALC_U(){

	for(j=2;j<=NJ-1;j++)
	{

          SP[j] = 0;
		  		    
//		  AW[j][i] = kw[j][i] * S[0][j][i] / (DXF[0][j][i]);
		  AS[j] = (nu + nu_t[j-1]) * S[1][j] / (DXF[1][j]);
//		  AE[j][i] = ke[j][i] *S[2][j][i] / (DXF[2][j][i]);
		  AN[j] = (nu+ nu_t[j]) * S[3][j] / (DXF[3][j]); 
		  AP[j] = AN[j] + AS[j] - SP[j];
		  
//		  cout<<AN[j]<<" "<<AS[j]<<endl;
		}	
		
	for(j=2;j<=NJ-1;j++)
	{
		U[j] = (AN[j] * U[j+1] + AS[j] * U[j-1] + ((1/(rho))*1*DELY[j-1]))/ AP[j];			
//       cout<<AS[j]<<" "<<AN[j]<<" "<<AP[j]<<" "<<U[j]<<endl;
	}
	cout<<"\n"<<endl;	
}

void CALC_Pk(){
	
	for(j=2;j<=NJ-1;j++)
	{
	
//		  kw[j][i] = 5 * (1 + (100/L) *Y_FACE_CENTER_X[j][i-1]);
//        Pk[j][i] = (nu_t[j+1][i]*pow(((U[j+1][i] - U[j][i])/DXF[3][j][i]),2)) - (nu_t[j-1][i]*pow(((U[j][i] - U[j-1][i])/DXF[1][j][i]),2));
      
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
		  AS[j] = (nu + (nu_t[j-1])/sigma_k) * S[1][j] / (DXF[1][j]);
//		  AE[j][i] = ke[j][i] *S[2][j][i] / (DXF[2][j][i]);
		  AN[j] = (nu+ (nu_t[j])/sigma_k) * S[3][j] / (DXF[3][j]);
		  AP[j] = AN[j] + AS[j] - SP[j] ;

         // cout<<kw[j][i]<<" "<<ks[j][i]<<" "<<ke[j][i]<<" "<<kn[j][i]<<endl;
       // cout<<AW[j][i]<<" "<<AS[j][i]<<" "<<AE[j][i]<<" "<<AN[j][i]<<" "<<AP[j][i]<<endl;
//		}	
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
		    
//		  AW[j][i] = kw[j][i] * S[0][j][i] / (DXF[0][j][i]);
		  AS[j] = (nu + (nu_t[j-1])/sigma_omega) * S[1][j] / (DXF[1][j]);
//		  AE[j][i] = ke[j][i] *S[2][j][i] / (DXF[2][j][i]);
		  AN[j] = (nu+ (nu_t[j])/sigma_omega) * S[3][j] / (DXF[3][j]);
		  AP[j] = AN[j] + AS[j] - SP[j] ;

         // cout<<kw[j][i]<<" "<<ks[j][i]<<" "<<ke[j][i]<<" "<<kn[j][i]<<endl;
       // cout<<AW[j][i]<<" "<<AS[j][i]<<" "<<AE[j][i]<<" "<<AN[j][i]<<" "<<AP[j][i]<<endl;

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
//    file5<<"VARIABLES = "<<'"'<<"Y"<<'"'<<endl;
//    file5<<"ZONE J="<<NJ-1<<", F=POINT"<<endl;

//	for (i=2;i<=NI-1;i++)
//	{
/*		
		for (j=2;j<=NJ-1;j++)	
		{
//       cout<<AW[j][i]<<" "<<AS[j][i]<<" "<<AE[j][i]<<" "<<AN[j][i]<<" "<<AP[j][i]<<endl;
		}
//	}
*/

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

/*
    //temperature values at CV and heat flux (face values)
	ofstream writer("structural_grid_04_HTDMA.dat");
	if (!writer){
		cout<<"error"<<endl;
	}
	else{
		writer<<"VARIABLES = "<<'"'<<"Y"<<'"'<<", "<<'"'<<"U"<<'"'<<", "<<'"'<<"k"<<'"'<<", "<<'"'<<"w"<<'"'<<", "<<'"'<<"nu_t"<<'"'<<", "<<'"'<<"Pk"<<'"'<<endl;
		writer<<"ZONE J="<<NJ<<", F=POINT"<<endl;
//     for (i=1;i<=NI;i++)	
//		{
			for (j=1;j<=NJ;j++)
			{
	 			writer<<Y_NODE[j]<<" "<<U[j]<<" "<<k[j]<<" "<<" "<<w[j]<<" "<<nu_t[j]<<" "<<Pk[j]<<endl;
			}
//		}
			writer.close();
	}
*/
}

/*
int FILE_WRITE2(int n){  
	int aa;
	char ch[80];
 
	char str_interface[80];

    aa = sprintf (ch, "%d",n);
			
		strcpy (str_interface,"2d_heat_conduction_with_heat_gen");
		strcat (str_interface,ch);
		strcat (str_interface,".dat");
		
  
         ofstream file(str_interface);
		file<<"VARIABLES = "<<'"'<<"X"<<'"'<<", "<<'"'<<"Y"<<'"'<<endl;
		file<<"ZONE I="<<NI<<", J="<<NI<<", F=POINT"<<endl;
     for(j=1;j<=NI;j++) 	
		{
			for(i=1;i<=NI;i++)
			{
				file<<X_NODE[j][i]<<" "<<Y_NODE[j][i]<<" "<<endl;
			}
		}
		//	file.close();
	}
*/
