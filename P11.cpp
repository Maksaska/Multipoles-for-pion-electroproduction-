#include <iostream> 
#include "header.h"
#include <complex>

using namespace std;

const double Mp(0.93827), Mpip(0.13957);

double phi_P11(const double& W, const double& Q2)
{
	double Value, arg;
		
	arg = (W - Mpip - Mp)/0.1;
	
	Value = -4.661*arg + 0.349*pow(arg,2) - 0.24*pow(arg,3);
	Value *= (1 + 9.977*Q2); 
	Value /= (1 + 2.384*Q2); 

	return Value;
}

complex<double> mult_P11(const double& W, const double& Q2)
{
	// P[:] = WR, Gamma_R, beta_pi; I[:] = 2*j, l, S(1535), n
	complex<double> buff;  
	vector<double> P; vector<int> I; 
	
	P.push_back(1.44); P.push_back(0.35); P.push_back(0.7); 
	I.push_back(1); I.push_back(1); I.push_back(0); I.push_back(1);
	
	buff = 0.078*f_gN(W,P[0],1,Q2)*f_PiN(W,P,I,Q2)*Gamma_tot(W,P,I)*P[0]*exp(1.0_i*phi_P11(W,Q2))*(k_mod(W, Q2)/k(W))*exp(-Q2/(6*0.229*0.229))*(P[0]*P[0] - W*W + 1.0_i*P[0]*Gamma_tot(W,P,I))/((pow(P[0]*P[0] - W*W,2) + pow(P[0]*Gamma_tot(W,P,I),2))*sqrt(3)); 	
	
	P.clear(); I.clear();
	 
	return buff;
}

void P11_table(vector<vector<double>>& Result, vector<double>& WQ2, int l_max) 
{
	complex<double> M;  
	vector<double> buff;
	double c; c = sqrt(2);
	
	for(double Q2 = WQ2[3]; Q2 <= WQ2[4]; Q2 = Q2 + WQ2[5])
	{
		for(double W = WQ2[0]; W <= WQ2[1]; W = W + WQ2[2])
		{
			M = mult_P11(W, Q2);		
			
			buff.push_back(W);
			buff.push_back(Q2);
			
			for(int i = 0; i < 6*(l_max+1); i++)
			{
				buff.push_back(0);
			}
			
			buff.push_back(0);
			buff.push_back(0);
			
			if(l_max > 0)
			{
				buff.push_back(c*M.real()); 
				buff.push_back(c*M.imag());
				for(int i = 4; i < 2*(l_max + 1); i++){buff.push_back(0);}			
			}
			
			for(int i = 0; i < 4*(l_max + 1); i++)
			{
				buff.push_back(0); 
			}
			
			for(int i = 0; i < 6*(l_max + 1); i++)
			{
				buff.push_back(0);
			}
			
			buff.push_back(0);
			buff.push_back(0);
			
			if(l_max > 0)
			{
				buff.push_back(M.real()); 
				buff.push_back(M.imag());
				for(int i = 4; i < 2*(l_max + 1); i++){buff.push_back(0);}		
			}
			
			for(int i = 0; i < 4*(l_max+1); i++)
			{
				buff.push_back(0); 
			}
			
			Result.push_back(buff);
			
			buff.clear();	
		}
	}	
}
