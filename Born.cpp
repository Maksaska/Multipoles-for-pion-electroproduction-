#include <iostream> 
#include "header.h"

using namespace std;

const double Mp(0.93827), Mn(0.93957), Mpip(0.13957), Mpiz(0.13498);

double F1p_f(const double& Q2)
{
	double F1p, Ge, Gm, tau, e; e = sqrt(4*M_PI/137.0);
	
	tau = Q2/(4*Mp*Mp);
	Ge = 1/((1 + Q2/0.71)*(1 + Q2/0.71));	
	Gm = 2.79/((1 + Q2/0.71)*(1 + Q2/0.71));
	F1p = e*(Ge + tau*Gm)/(1 + tau); 
	
	return F1p;
}

double F2p_f(const double& Q2)
{
	double F2p, Ge, Gm, tau, e; e = sqrt(4*M_PI/137.0);
	
	tau = Q2/(4*Mp*Mp);
	Ge = 1/((1 + Q2/0.71)*(1 + Q2/0.71));	
	Gm = 2.79/((1 + Q2/0.71)*(1 + Q2/0.71));
	F2p = e*1.79*(Gm - Ge)/(2*Mp*1.79*(1 + tau)); 
	
	return F2p;
}

double F1n_f(const double& Q2)
{
	double F1n, Ge, Gm, tau, e; e = sqrt(4*M_PI/137.0);
	
	tau = Q2/(4*Mn*Mn);	
	Gm = -1.91/((1 + Q2/0.71)*(1 + Q2/0.71));
	Ge = -Gm*tau/(1 + 5.6*tau);
	F1n = e*(Ge + tau*Gm)/(1 + tau);		
	
	return F1n;
}

double F2n_f(const double& Q2)
{
	double F2n, Ge, Gm, tau, e; e = sqrt(4*M_PI/137.0);
	
	tau = Q2/(4*Mn*Mn);	
	Gm = -1.91/((1 + Q2/0.71)*(1 + Q2/0.71));
	Ge = -Gm*tau/(1 + 5.6*tau);
	F2n = -1.91*e*(Gm - Ge)/(-1.91*2*Mp*(1 + tau));	
	
	return F2n;
}

double Fpi_f(const double& Q2)
{
	double Fpi;
	
	Fpi = 0.85/(1 + 0.431*Q2/(6*0.04)) + 0.15/((1 + 0.411*Q2/(12*0.04))*(1 + 0.411*Q2/(12*0.04)));	
	
	return Fpi;
}

double F1v_f(const double& Q2){ return F1p_f(Q2) - F1n_f(Q2);}
double F2v_f(const double& Q2){ return F2p_f(Q2) - F2n_f(Q2);}
double F1s_f(const double& Q2){ return F1p_f(Q2) + F1n_f(Q2);}
double F2s_f(const double& Q2){ return F2p_f(Q2) + F2n_f(Q2);}

vector<double> Gamma(const double& Q2, const double& t, const int& index)
{
	vector<double> G;
	double buff, g;
	
	g = sqrt(4*M_PI*14.4);	
	
	if(index == 1 or index == -1)
	{
		buff = g*F1v_f(Q2)/2;
		G.push_back(buff);
		
		buff = -g*F1v_f(Q2)/(t - Mpip*Mpip);
		G.push_back(buff);
		
		buff = -g*F2v_f(Q2)/2; 
		G.push_back(buff);
		G.push_back(buff);
		
		buff = -0.5*g*F1v_f(Q2)/(t - Mpip*Mpip);
		G.push_back(buff);
		G.push_back(0);
		
		buff = 2*g*(Fpi_f(Q2) - F1v_f(Q2))/Q2;
		G.push_back(buff);			
	} else 
	{
		buff = g*F1s_f(Q2)/2;
		G.push_back(buff);
		
		buff = -g*F1s_f(Q2)/(t - Mpip*Mpip);
		G.push_back(buff);
		
		buff = -g*F2s_f(Q2)/2; 
		G.push_back(buff);
		G.push_back(buff);
		
		buff = -0.5*g*F1s_f(Q2)/(t - Mpip*Mpip);
		G.push_back(buff);
		G.push_back(0);
		
		buff = 2*g*(Fpi_f(Q2) - F1v_f(Q2))/Q2;
		G.push_back(buff);
	}
	
	return G;
}

int ksi(const int& index)
{
	if(index == 1 or index == 0)
	{
		return 1;
	} 
	return -1;
}

double Q(const int& n, const double& x)
{
	double Value;
	
	if (n == 0)
	{
		return log((x + 1)/(x - 1))/2;
	} else if (n == 1)
	{
		return x*log((x + 1)/(x - 1))/2 - 1;
	} else if (n < 0)
	{
		return 0;
	} else 
	{
		Value = (2*n - 1)*x*Q(n-1, x)/n - (n - 1)*Q(n-2, x)/n;
	}
	
	return Value;
}

int Parity(const int& m)
{
	if(m % 2 == 0){ return 1;}
	return -1;
}

double RN(const int& l, const double& E2_barred)
{
	return Parity(l)*(  Q(l+1, E2_barred) - Q(l-1, E2_barred) )/(2*l+1);
}

double RPi(const int& l, const double& q0_barred)
{
	return (  Q(l+1, q0_barred) - Q(l-1, q0_barred) )/(2*l+1);
}

double T(const int& l, const double& Q2, const double& W, const int& pm_index)
{
	double k, k0, q, q0, E1, E2, E2_barred;
	
	k0 = (W*W - Q2 - Mp*Mp)/(2*W); q0 = (W*W + Mpip*Mpip - Mp*Mp)/(2*W);
	
	k = sqrt(k0*k0 + Q2); q = sqrt(q0*q0 - Mpip*Mpip);
	
	E1 = (W*W + Mp*Mp + Q2)/(2*W); E2 = (W*W - Mpip*Mpip + Mp*Mp)/(2*W);
	
	E2_barred = (2*k0*E2 + Q2)/(2*k*q);
	
	if(pm_index == 1)
	{
		return Parity(l)*(Q(l, E2_barred)/(k*q) - (W + Mp)*Q(l+1, E2_barred)/((E1 + Mp)*(E2 + Mp)*(W - Mp)));	
	} else 
	{
		return Parity(l)*(Q(l, E2_barred)/(k*q) - (W + Mp)*Q(l-1, E2_barred)/((E1 + Mp)*(E2 + Mp)*(W - Mp)));
	}
}

int delta_f(const int& q1,const int& q2)
{
	if(q1 == q2){ return 1;}
	return 0;
}

vector<vector<double>> EMS(const double& W, const double& Q2, const double& theta, const int& limiter, const int& id)
{
	vector<vector<double>> V;
	vector<double> buff;
	vector<double> G_f; 
	double val, Const(1); 
	
	double k, k0, q, q0, E1, E2, E2_barred, q0_barred, t;
	
	k0 = (W*W - Q2 - Mp*Mp)/(2*W); q0 = (W*W + Mpip*Mpip - Mp*Mp)/(2*W);
	
	k = sqrt(k0*k0 + Q2); q = sqrt(q0*q0 - Mpip*Mpip);
	
	E1 = (W*W + Mp*Mp + Q2)/(2*W); E2 = (W*W - Mpip*Mpip + Mp*Mp)/(2*W);
	
	E2_barred = (2*k0*E2 + Q2)/(2*k*q); 
	
	q0_barred = (2*k0*q0 + Q2)/(2*k*q); 
	
	t = -2*k*q*(q0_barred - cos(theta)) + Mpip*Mpip;
	
	G_f = Gamma(Q2, t, id);
	
	for(int l = 0; l < limiter; l++) // El+
	{
		val = Const*(W - Mp)*(sqrt((E1 + Mp)*(E2 + Mp))/(4*Mp*(l+1)))*(  2*delta_f(l, 0)*(G_f[0] + (W - Mp)*G_f[2] - ksi(id)*(W + Mp)*G_f[2])/(W*W - Mp*Mp) - ksi(id)*(G_f[0] - 2*Mp*G_f[2])*T(l, Q2, W, 1) - 0.5*(1 - ksi(id))*(4*G_f[0] + Q2*G_f[6])*(l*RPi(l, q0_barred)/((E1 + Mp)*(W - Mp)) - q*(l+1)*RPi(l+1, q0_barred)/(k*(E2 + Mp)*(W - Mp))) - 2*ksi(id)*(G_f[0] - (W + Mp)*G_f[2])*l*RN(l, E2_barred)/((E1+Mp)*(W - Mp)) + 2*ksi(id)*(G_f[0] + (W - Mp)*G_f[2])*q*(l+1)*RN(l+1, E2_barred)/(k*(E2 + Mp)*(W-Mp))  );
		buff.push_back(val);
	}
	V.push_back(buff);
	buff.clear();
	
	for(int l = 0; l < limiter; l++) // El-
	{
		if(l != 0)
		{
		val = Const*(W - Mp)*(sqrt((E1 + Mp)*(E2 + Mp))/(4*Mp*l))*( - ksi(id)*(G_f[0] - 2*Mp*G_f[2])*T(l, Q2, W, -1) + 0.5*(1 - ksi(id))*(4*G_f[0] + Q2*G_f[6])*((l+1)*RPi(l, q0_barred)/((E1 + Mp)*(W - Mp)) - q*l*RPi(l-1, q0_barred)/(k*(E2 + Mp)*(W - Mp))) + 2*ksi(id)*(G_f[0] - (W + Mp)*G_f[2])*(l+1)*RN(l, E2_barred)/((E1+Mp)*(W - Mp)) - 2*ksi(id)*(G_f[0] + (W - Mp)*G_f[2])*q*l*RN(l-1, E2_barred)/(k*(E2 + Mp)*(W-Mp))  );
		} else {val = 0;}
		buff.push_back(val);
	}
	V.push_back(buff);
	buff.clear();
	
	for(int l = 0; l < limiter; l++) // Ml+
	{
val = Const*(W - Mp)*(sqrt((E1 + Mp)*(E2 + Mp))/(4*Mp*(l+1)))*( - ksi(id)*(G_f[0] - 2*Mp*G_f[2])*T(l, Q2, W, 1) + 0.5*(1 - ksi(id))*(4*G_f[0] + Q2*G_f[6])*RPi(l, q0_barred)/((E1 + Mp)*(W - Mp)) + 2*ksi(id)*(G_f[0] - (W + Mp)*G_f[2])*RN(l, E2_barred)/((E1+Mp)*(W - Mp)) );
		buff.push_back(val);
	}
	V.push_back(buff);
	buff.clear();
	
	for(int l = 0; l < limiter; l++) // Ml-
	{
		if(l != 0)
		{
val = Const*(W - Mp)*(sqrt((E1 + Mp)*(E2 + Mp))/(4*Mp*l))*(   2*q*k*delta_f(l, 1)*(-G_f[0] + (W + Mp)*G_f[2] - ksi(id)*(W - Mp)*G_f[2])/((E1 + Mp)*(E2 + Mp)*(W - Mp)*(W - Mp)) + ksi(id)*(G_f[0] - 2*Mp*G_f[2])*T(l, Q2, W, -1) - 0.5*(1 - ksi(id))*(4*G_f[0] + Q2*G_f[6])*RPi(l, q0_barred)/((E1 + Mp)*(W - Mp)) - 2*ksi(id)*(G_f[0] - (W + Mp)*G_f[2])*RN(l, E2_barred)/((E1+Mp)*(W - Mp))  );
		} else {val = 0;}
		buff.push_back(val);
	}
	V.push_back(buff);
	buff.clear();
	
	for(int l = 0; l < limiter; l++) // Sl+
	{
		val = Const*(W - Mp)*(sqrt((E1 + Mp)*(E2 + Mp))/(4*Mp*(l+1)))*(2*k*delta_f(l, 0)*(-G_f[0] + (E1 + Mp)*G_f[2] + ksi(id)*(W + Mp)*G_f[2] + 0.25*(1 - ksi(id))*k0*(W + Mp)*G_f[6])/((E1 + Mp)*(W*W - Mp*Mp)) + 0.5*(1 - ksi(id))*(2*q0 - k0)*(0.5*Q2*G_f[6] + 2*G_f[0])*(Q(l, q0_barred)/(q*(E1 + Mp)) - Q(l+1, q0_barred)/(k*(E2 + Mp)))/(W - Mp) + Parity(l+1)*ksi(id)*(-Mp*(2*q0 - k0)*G_f[2] + (2*q0 - W)*G_f[0])*(Q(l, E2_barred)/(q*(E1 + Mp)) + Q(l+1, E2_barred)/(k*(E2 + Mp)))/(W - Mp) + Parity(l+1)*ksi(id)*(Mp*G_f[0] + (k0*W + Q2 - Mpip*Mpip)*G_f[2])*(Q(l, E2_barred)/(q*(E1 + Mp)) - Q(l+1, E2_barred)/(k*(E2 + Mp)))/(W - Mp));
		buff.push_back(val);
	}
	V.push_back(buff);
	buff.clear();
	
	for(int l = 0; l < limiter; l++) // Sl-
	{
		if(l != 0)
		{
		val = Const*(W - Mp)*(sqrt((E1 + Mp)*(E2 + Mp))/(4*Mp*l))*(2*q*delta_f(l, 1)*(G_f[0] + (E1 - Mp)*G_f[2] + ksi(id)*(W - Mp)*G_f[2] - 0.25*(1 - ksi(id))*k0*(W - Mp)*G_f[6])/((E2 + Mp)*(W - Mp)*(W - Mp)) + 0.5*(1 - ksi(id))*(2*q0 - k0)*(0.5*Q2*G_f[6] + 2*G_f[0])*(Q(l, q0_barred)/(q*(E1 + Mp)) - Q(l-1, q0_barred)/(k*(E2 + Mp)))/(W - Mp) + Parity(l+1)*ksi(id)*(-Mp*(2*q0 - k0)*G_f[2] + (2*q0 - W)*G_f[0])*(Q(l, E2_barred)/(q*(E1 + Mp)) + Q(l-1, E2_barred)/(k*(E2 + Mp)))/(W - Mp) + Parity(l+1)*ksi(id)*(Mp*G_f[0] + (k0*W + Q2 - Mpip*Mpip)*G_f[2])*(Q(l, E2_barred)/(q*(E1 + Mp)) - Q(l-1, E2_barred)/(k*(E2 + Mp)))/(W - Mp));
		} else {val = 0;}
		buff.push_back(val);
	}
	V.push_back(buff);
	buff.clear(); G_f.clear();

	return V;
}

void Born_table(vector<vector<double>>& Result, vector<double>& WQ2, int l_max)
{
	vector<vector<double>> M1, M2, M3;  
	vector<double> buff;
	double theta(0), c; c = sqrt(2);
	
	for(double Q2 = WQ2[3]; Q2 <= WQ2[4]; Q2 = Q2 + WQ2[5])	
	{
		for(double W = WQ2[0]; W <= WQ2[1]; W = W + WQ2[2])
		{
			M1 = EMS(W, Q2, theta, l_max+1, 0);			
			M2 = EMS(W, Q2, theta, l_max+1, 1);			
			M3 = EMS(W, Q2, theta, l_max+1, -1);
			
			buff.push_back(W);
			buff.push_back(Q2);
			
			for (long unsigned int i = 0; i < M1.size(); i++)
			{
				for(long unsigned int j = 0; j < M1[i].size(); j++)
				{
					buff.push_back(c*(M1[i][j] + M3[i][j]));
					buff.push_back(0);	
				}
			}
		
			for (long unsigned int i = 0; i < M1.size(); i++)
			{
				for(long unsigned int j = 0; j < M1[i].size(); j++)
				{
					buff.push_back(M1[i][j] + M2[i][j]);
					buff.push_back(0);	
				}
			}
			
			Result.push_back(buff);
			
			buff.clear();
			M1.clear(); M2.clear(); M3.clear();	
		}
	}	
}
