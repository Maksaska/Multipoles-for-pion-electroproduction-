#include <iostream> 
#include "header.h"

using namespace std;

const double Mp(0.93827), Mn(0.93957), Mpip(0.13957), Mpiz(0.13498), X(0.5), Meta(0.547);

double q(const double& W)
{
	return sqrt((W*W + Mpip*Mpip - Mp*Mp)*(W*W + Mpip*Mpip - Mp*Mp)/(4*W*W) - Mpip*Mpip);
}

double q2p(const double& W)
{
	return sqrt((W*W + 4*Mpip*Mpip - Mp*Mp)*(W*W + 4*Mpip*Mpip - Mp*Mp)/(4*W*W) - 4*Mpip*Mpip);
}

double qeta(const double& W)
{
	return sqrt((W*W + Meta*Meta - Mp*Mp)*(W*W + Meta*Meta - Mp*Mp)/(4*W*W) - Meta*Meta);
}

double k(const double& W)
{
	return (W*W - Mp*Mp)/(2*W);
}

double k_mod(const double& W, const double& Q2)
{
	return sqrt(Q2 + pow(W*W - Mp*Mp - Q2, 2)/(4*W*W));
}

double Gamma_PiN(const double& W, vector<double>& P, const int& l)
{
	return P[1]*P[2]*pow((q(W)/q(P[0])), 2*l+1)*pow(((X*X + q(P[0])*q(P[0]))/(X*X + q(W)*q(W))), l)*P[0]/W;
}

double Gamma_in(const double& W, vector<double>& P, const int& l)
{
	if(W >= 1.22){return (1 - P[2])*P[1]*pow((q2p(W)/q2p(P[0])),2*l+4)*pow(((X*X + q2p(P[0])*q2p(P[0]))/(X*X + q2p(W)*q2p(W))),l+2);}
	return 0;
}

double Gamma_eta(const double& W, vector<double>& P, const int& l)
{
	return P[1]*P[2]*pow((qeta(W)/qeta(P[0])), 2*l+1)*pow(((X*X + qeta(P[0])*qeta(P[0]))/(X*X + qeta(W)*qeta(W))), l)*P[0]/W;
}

double Gamma_tot(const double& W, vector<double>& P, vector<int>& I)
{
	if(I[2] == 1)
	{
		return Gamma_PiN(W, P, I[1]) + Gamma_in(W, P, I[1]) + Gamma_eta(W, P, I[1]);
	} 
	
	return Gamma_PiN(W, P, I[1]) + Gamma_in(W, P, I[1]);
}

double f_PiN(const double& W, vector<double>& P, vector<int>& I, const double& Q2) 
{ 
	return sqrt(Gamma_PiN(W, P, I[1])*(k_mod(W,Q2)/q(W))*(Mp/W)/((I[0]+1)*M_PI*Gamma_tot(W, P, I)*Gamma_tot(W, P, I)));
}

double f_gN(const double& W, const double& WR, const int& n, const double& Q2) 
{
	return pow(k_mod(W,Q2)/k_mod(WR,Q2), n)*((X*X + k_mod(WR,Q2)*k_mod(WR,Q2))/(X*X + k_mod(W,Q2)*k_mod(W,Q2)));
}
