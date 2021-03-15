#include <iostream> 
#include "header.h"

using namespace std;

const double Mp(0.93827), Mpip(0.13957);

double FD(const double& Q2)
{
	return 1/((1 + Q2/0.71)*(1 + Q2/0.71));
}

double AB(const double& W, vector<string>& MES)
{
	string path("final.csv");
	vector<vector<double>> Values;
	vector<double> x, y, s; 
	Reading(path,Values);
	double buff, xx; xx = W*1000; int id_y;
	
	if(MES[0] == "M" and MES[1] == "A"){id_y = 1; buff = 0.6471;}
	if(MES[0] == "M" and MES[1] == "B"){id_y = 3; buff = 0.7963;}
	if(MES[0] == "E" and MES[1] == "A"){id_y = 5; buff = 1.604;}
	if(MES[0] == "E" and MES[1] == "B"){id_y = 7; buff = 1.5;}
	if(MES[0] == "L" and MES[1] == "A"){id_y = 9; buff = 0.2979;}
	if(MES[0] == "L" and MES[1] == "B"){id_y = 11; buff = 0.21;}
	
	for(long unsigned int i = 0; i < Values[0].size(); i++)
	{
		x.push_back(Values[0][i]);
		y.push_back(Values[id_y][i]);
		s.push_back(Values[id_y+1][i]);
	}
	
	for(long unsigned int i = 1; i < x.size(); i++)
	{
		if(W >= x[i-1] and W <= x[i])
		{
			buff = s[i-1]*pow(x[i] - xx,3)/(6*(x[i] - x[i-1])) + s[i]*pow(xx - x[i-1],3)/(6*(x[i] - x[i-1])) + (y[i-1] - s[i-1]*pow(x[i] - x[i-1], 2)/6)*(x[i] - xx)/(x[i] - x[i-1]) + (y[i] - s[i]*pow(x[i] - x[i-1], 2)/6)*(xx - x[i-1])/(x[i] - x[i-1]);
			Values.clear(); x.clear(); y.clear(); s.clear();
			return buff;
		}
	}
	
	Values.clear(); x.clear(); y.clear(); s.clear();
	return buff;
}

double phi_delta(const double& W, const double& Q2, const string& mult)
{
	double Value, arg; string buff; 
	vector<string> V;
	
	arg = (W - Mpip - Mp)/0.1;
	
	if(mult == "M")
	{
		buff = "M"; V.push_back(buff);
		Value = 22.13*arg - 3.769*pow(arg,2) + 0.184*pow(arg,3);
		buff = "A"; V.push_back(buff); 
		Value *= (1 + AB(W, V)*Q2); V.clear(); 
		buff = "M"; V.push_back(buff); buff = "B"; V.push_back(buff);
		Value /= (1 + AB(W, V)*Q2); 
	}
	
	if(mult == "E")
	{
		buff = "E"; V.push_back(buff);
		Value = 83.336*arg - 28.457*pow(arg,2) + 3.356*pow(arg,3) - 0.122*pow(arg,4);
		buff = "A"; V.push_back(buff);
		Value *= (1 + AB(W, V)*Q2); V.clear();
		buff = "E"; V.push_back(buff); buff = "B"; V.push_back(buff);
		Value /= (1 + AB(W, V)*Q2);	
	}
	
	if(mult == "L")
	{
		buff = "L"; V.push_back(buff);
		Value = 43.668*arg - 2.872*pow(arg,2) - 1.91*pow(arg,3) + 0.23*pow(arg,4);
		buff = "A"; V.push_back(buff);
		Value *= (1 + AB(W, V)*Q2); V.clear();
		buff = "L"; V.push_back(buff); buff = "B"; V.push_back(buff);
		Value /= (1 + AB(W, V)*Q2);	
	}
	
	V.clear();
	return Value;
}

vector<complex<double>> mult_delta(const double& W, const double& Q2) // M first
{
	vector<complex<double>> res; // P[:] = WR, Gamma_R, beta_pi; I[:] = 2*j, l, S(1535), n
	complex<double> buff; string mult("M"); 
	vector<double> P; vector<int> I; 
	
	P.push_back(1.235); P.push_back(0.13); P.push_back(1); 
	I.push_back(3); I.push_back(1); I.push_back(0); I.push_back(2);
	
	buff = 323*f_gN(W,P[0],2)*f_PiN(W,P,I)*sqrt(1.5)*Gamma_tot(W,P,I)*P[0]*exp(1i*phi_delta(W,Q2,mult))*(k_mod(W, Q2)/k(W))*exp(-0.24*Q2)*FD(Q2)/(P[0]*P[0] - W*W - 1i*P[0]*Gamma_tot(W,P,I)); 	
	res.push_back(buff); I.clear();
	
	I.push_back(3); I.push_back(1); I.push_back(0); I.push_back(1); mult = "E";
	buff = -17*f_gN(W,P[0],1)*f_PiN(W,P,I)*sqrt(1.5)*Gamma_tot(W,P,I)*P[0]*exp(1i*phi_delta(W,Q2,mult))*(k_mod(W, Q2)/k(W))*exp(-0.24*Q2)*FD(Q2)/(P[0]*P[0] - W*W - 1i*P[0]*Gamma_tot(W,P,I));
	res.push_back(buff);
	
	mult = "L";
	buff = -17*f_gN(W,P[0],1)*f_PiN(W,P,I)*sqrt(1.5)*Gamma_tot(W,P,I)*P[0]*exp(1i*phi_delta(W,Q2,mult))*(k_mod(W, Q2)/k(W))*exp(-0.24*Q2)*FD(Q2)/(P[0]*P[0] - W*W - 1i*P[0]*Gamma_tot(W,P,I)); // This is S1+
	res.push_back(buff);
	
	P.clear(); I.clear();
	 
	return res;
}

void P33_table(vector<vector<double>>& Result, vector<double>& WQ2, int l_max)
{
	vector<complex<double>> M;  
	vector<double> buff;
	double c; c = sqrt(2);
	
	for(double W = WQ2[0]; W <= WQ2[1]; W = W + WQ2[2])
	{
		for(double Q2 = WQ2[3]; Q2 <= WQ2[4]; Q2 = Q2 + WQ2[5])
		{
			M = mult_delta(W, Q2);		
			
			buff.push_back(W);
			buff.push_back(Q2);
			buff.push_back(0);
			buff.push_back(0);
			
			if(l_max > 0)
			{
				buff.push_back(-c*M[1].real()/3); // E
				buff.push_back(-c*M[1].imag()/3);
				for(int i = 4; i < 2*(l_max + 1); i++){buff.push_back(0);}			
			}
			
			for(int i = 0; i < 2*(l_max + 1); i++)
			{
				buff.push_back(0); 
			}
			
			buff.push_back(0);
			buff.push_back(0);
			
			if(l_max > 0)
			{
				buff.push_back(-c*M[0].real()/3); // M
				buff.push_back(-c*M[0].imag()/3);
				for(int i = 4; i < 2*(l_max+1); i++){buff.push_back(0);}			
			}
			
			for(int i = 0; i < 2*(l_max+1); i++)
			{
				buff.push_back(0); 
			}
			
			buff.push_back(0);
			buff.push_back(0);
			
			if(l_max > 0)
			{
				buff.push_back(-c*M[2].real()/3); // S
				buff.push_back(-c*M[2].imag()/3);
				for(int i = 4; i < 2*(l_max+1); i++){buff.push_back(0);}			
			}
			
			for(int i = 0; i < 2*(l_max+1); i++)
			{
				buff.push_back(0); 
			}			
			
			buff.push_back(0);
			buff.push_back(0);
			
			if(l_max > 0)
			{
				buff.push_back(2*M[1].real()/3); // E
				buff.push_back(2*M[1].imag()/3);
				for(int i = 4; i < 2*(l_max+1); i++){buff.push_back(0);}			
			}
			
			for(int i = 0; i < 2*(l_max+1); i++)
			{
				buff.push_back(0); 
			}
			
			buff.push_back(0);
			buff.push_back(0);
			
			if(l_max > 0)
			{
				buff.push_back(2*M[0].real()/3); // M
				buff.push_back(2*M[0].imag()/3);
				for(int i = 4; i < 2*(l_max+1); i++){buff.push_back(0);}			
			}
			
			for(int i = 0; i < 2*(l_max+1); i++)
			{
				buff.push_back(0); 
			}
			
			buff.push_back(0);
			buff.push_back(0);
			
			if(l_max > 0)
			{
				buff.push_back(2*M[2].real()/3); // S
				buff.push_back(2*M[2].imag()/3);
				for(int i = 4; i < 2*(l_max+1); i++){buff.push_back(0);}			
			}
			
			for(int i = 0; i < 2*(l_max+1); i++)
			{
				buff.push_back(0); 
			}
			
			Result.push_back(buff);
			
			buff.clear();
			M.clear(); 	
		}
	}	
}
