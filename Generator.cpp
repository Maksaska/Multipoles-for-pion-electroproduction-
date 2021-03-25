#include <iostream> 
#include <fstream>
#include "header.h"

using namespace std;

const double Mp(0.93827), Mn(0.93957), Mpip(0.13957), Mpiz(0.13498);
const int polarization(0); // Beam polarizataion 
int ch(0); 
vector<bool> P1;

double P(const int& der,const int& n, const double& theta) 
{
	double Value(0);

	switch(n)
	{
	
		case 0: if(der == 0)
			{
				Value = 1;
			}
			else
			{
				Value = 0;
			}
			break;

		case 1: if(der == 0){Value = cos(theta);} 
			if(der == 1){Value = 1;}
			if(der == 2){Value = 0;}
			break;

		case 2: if(der == 0){Value = 0.5*(3*pow(cos(theta),2) - 1);} 
			if(der == 1){Value = 3*cos(theta);}
			if(der == 2){Value = 3;}
			break;

		case 3: if(der == 0){Value = 0.5*(5*pow(cos(theta),3) - 3*cos(theta));} 
			if(der == 1){Value = 7.5*pow(cos(theta),2)-1.5;}
			if(der == 2){Value = 15*cos(theta);}
			break;

		case 4: if(der == 0){Value = 0.125*(35*pow(cos(theta),4) - 30*pow(cos(theta),2) + 3);} 
			if(der == 1){Value = 17.5*pow(cos(theta),3) - 7.5*cos(theta);}
			if(der == 2){Value = 52.5*pow(cos(theta),2) - 7.5;}
			break;

		case 5: if(der == 0){Value = 0.125*(63*pow(cos(theta),5) - 70*pow(cos(theta),3) + 15*cos(theta));} 
			if(der == 1){Value = 315*pow(cos(theta),4)/8 - 210*pow(cos(theta),2)/8 + 15/8;}
			if(der == 2){Value = 315*pow(cos(theta),3)/2 - 210*cos(theta)/4;}
			break;

		case 6: if(der == 0){Value = 231*pow(cos(theta),6)/16 - 315*pow(cos(theta),4)/16 + 105*pow(cos(theta),2)/16 - 5/16;} 
			if(der == 1){Value = 3*231*pow(cos(theta),5)/8 - 315*pow(cos(theta),3)/4 + 105*cos(theta)/8  ;}
			if(der == 2){Value = 15*231*pow(cos(theta),4)/8 - 945*pow(cos(theta),2)/4 + 105/8;}
			break;

		default: Value = 0; break;
	}
	
	return Value;
}

vector<complex<double>> Merger(const double& W,  const double& Q2)
{
	vector<complex<double>> AMP;
	complex<double> buff;
	
	for(int i = 0; i < 42; i++)
	{
		buff = 0;
		AMP.push_back(buff);
	}
	
	if(ch == 1 and P1[0]) // Pi+n born
	{
		vector<vector<double>> M1, M2, M3; 
		
		M1 = EMS(W, Q2, 0, 6, 0);			
		M2 = EMS(W, Q2, 0, 6, 1);			
		M3 = EMS(W, Q2, 0, 6, -1);
	
		for (long unsigned int i = 0; i < M1.size(); i++)
		{
			for(long unsigned int j = 0; j < M1[i].size(); j++)
			{
				AMP[i*7 + j] += sqrt(2)*(M1[i][j] + M3[i][j]);	
			}
		}
		
		M1.clear(); M2.clear(); M3.clear();	
			
	} else if(ch == 0 and P1[0]) // Pi0p born
	{
		vector<vector<double>> M1, M2, M3;  
		
		M1 = EMS(W, Q2, 0, 6, 0);			
		M2 = EMS(W, Q2, 0, 6, 1);			
		M3 = EMS(W, Q2, 0, 6, -1);
		
		for (long unsigned int i = 0; i < M1.size(); i++)
		{
			for(long unsigned int j = 0; j < M1[i].size(); j++)
			{
				AMP[i*7 + j] += M1[i][j] + M2[i][j];	
			}
		}
		
		M1.clear(); M2.clear(); M3.clear();		
	} 
	
	if(ch == 1 and P1[1]) // Pi+n P33(1232)
	{
		vector<complex<double>> M; 
		M = mult_delta(W, Q2);
		AMP[1] += (-sqrt(2)*M[1]/3.0);
		AMP[15] += (-sqrt(2)*M[0]/3.0);
		AMP[29] += (-sqrt(2)*M[2]/3.0);		
		M.clear();
	} else if(ch == 0 and P1[1]) // Pi0p P33(1232)
	{
		vector<complex<double>> M; 
		M = mult_delta(W, Q2);
		AMP[1] += (2.0*M[1]/3.0);
		AMP[15] += (2.0*M[0]/3.0);
		AMP[29] += (2.0*M[2]/3.0);		
		M.clear();	
	}
	
	if(ch == 1 and P1[2]) // Pi+n P11(1440)
	{ 
		buff = mult_P11(W, Q2); 
		AMP[22] += sqrt(2)*buff;
	} else if(ch == 0 and P1[2]) // Pi0p P11(1440)
	{
		buff = mult_P11(W, Q2); 
		AMP[22] += buff;	
	} 
		
	return AMP;
}

vector<complex<double>> Helicity_amplitudes(const double& W,  const double& Q2, const double& theta)
{
	vector<complex<double>> AMP, H; double l;
	complex<double> buff; buff = 0;
	
	AMP = Merger(W, Q2);
	
	for(int i = 0; i < 6; i++) // H1
	{
		l = double(i);
		buff += (AMP[l] - AMP[l+14] - AMP[l+8] - AMP[l+22])*(P(2,i,theta) - P(2,i+1,theta));
	} 
	buff = buff*sin(theta)*cos(theta/2)/sqrt(2);
	H.push_back(buff); buff = 0;
	
	for(int i = 0; i < 6; i++) // H2
	{
		l = double(i);
		buff += ((l+2)*AMP[l] + l*AMP[l+14] + l*AMP[l+8] - (l+2)*AMP[l+22])*(P(1,i,theta) - P(1,i+1,theta));
	} 
	buff = buff*cos(theta/2)/sqrt(2);
	H.push_back(buff); buff = 0;
	
	for(int i = 0; i < 6; i++) // H3
	{
		l = double(i);
		buff += (AMP[l] - AMP[l+14] + AMP[l+8] + AMP[l+22])*(P(2,i,theta) + P(2,i+1,theta));
	} 
	buff = buff*sin(theta)*sin(theta/2)/sqrt(2);
	H.push_back(buff); buff = 0;
	
	for(int i = 0; i < 6; i++) // H4
	{
		l = double(i);
		buff += ((l+2)*AMP[l] + l*AMP[l+14] - l*AMP[l+8] + (l+2)*AMP[l+22])*(P(1,i,theta) + P(1,i+1,theta));
	} 
	buff = buff*sin(theta/2)/sqrt(2);
	H.push_back(buff);buff = 0;
	
	for(int i = 0; i < 6; i++) // H5
	{
		l = double(i);
		buff += (l+1)*(AMP[l+28] + AMP[l+36])*(P(1,i,theta) - P(1,i+1,theta));
	} 
	buff = buff*cos(theta/2)*sqrt(Q2)/k_mod(W,Q2);
	H.push_back(buff);buff = 0;
	
	for(int i = 0; i < 6; i++) // H6
	{
		l = double(i);
		buff += (l+1)*(AMP[l+28] - AMP[l+36])*(P(1,i,theta) + P(1,i+1,theta));
	} 
	buff = buff*sin(theta/2)*sqrt(Q2)/k_mod(W,Q2);
	H.push_back(buff); 
	
	AMP.clear();
	return H;
}

double Section(const double& W,  const double& Q2, const double& theta,  const double& phi, const double& E0)
{
	double S_t, S_l, S_tt, S_lt, S_lt_pr, S, Gamma_flux, eps, nu;
	vector<complex<double>> H;
	
	nu =  (W*W + Q2 - Mp*Mp)/(2*Mp);
	eps = 1/(1 + 2*(nu*nu + Q2)/(4*(E0 - nu)*E0 - Q2));		
	Gamma_flux = (E0 - nu)*(W*W - Mp*Mp)/(137*4*M_PI*M_PI*Mp*E0*(1 - eps)*Q2);
	
	H = Helicity_amplitudes(W, Q2, theta);
	
	S_t = q(W)*(abs(H[0])*abs(H[0]) + abs(H[1])*abs(H[1]) + abs(H[2])*abs(H[2]) + abs(H[3])*abs(H[3]))/(2*k(W));
	S_l = q(W)*(abs(H[4])*abs(H[4]) + abs(H[5])*abs(H[5]))/k(W);
	S_tt = q(W)*real(H[2]*conj(H[1]) - H[3]*conj(H[0]))/k(W);
	S_lt = -q(W)*real((H[0] - H[3])*conj(H[4]) + (H[1] + H[2])*conj(H[5]))/(sqrt(2)*k(W));
	S_lt_pr = -q(W)*imag((H[0] - H[3])*conj(H[4]) + (H[1] + H[2])*conj(H[5]))/(sqrt(2)*k(W));

	S = S_t + eps*S_l + eps*S_tt*cos(2*phi) + sqrt(2*eps*(1+eps))*S_lt*cos(phi) + polarization*S_lt_pr*sin(phi);

	Gamma_flux = 1;

	S = Gamma_flux*S; 
	
	return S;
}

void Generate_lund(vector<double>& V, int channel)
{
	double theta, W, Q2, x, y, z, S, phi, ang1, ang2, mpi, mp, Epi, Ep, p, nu; int N(-1), NN; bool value;
	TLorentzVector e, adron, meson, gamma1, gamma2; ch = channel;
	TVector3 beta;
	char FileName[100];
	
	double E0(6.5), Radius_c(10); //Target's radius 10 cm; Beam energy 6.5 GeV
	
	srand(time(NULL));
	
	if(channel == 1)
	{	
		mpi = Mpip; mp = Mn;
		sprintf(FileName,"pin_W_%g_%g_Q2_%g_%g_.lund", V[0], V[1], V[2], V[3]);
	}
	else
	{	
		mpi = Mpiz; mp = Mp;
		sprintf(FileName,"pi0p_W_%g_%g_Q2_%g_%g_.lund", V[0], V[1], V[2], V[3]);
	}
	
	cout << "\n";
	
	cout << "What are we going to include?\nBorn terms:";
	value = Response(); P1.push_back(value); if(value){N++;}
	
	cout << "P33(1232):"; 
	value = Response(); P1.push_back(value); if(value){N++;}
	
	cout << "P11(1440):"; 
	value = Response(); P1.push_back(value); if(value){N++;}	
		
	if(N == -1){return;} N = -1;
		
	cout << "How many events?\nN = ";
	
	while(N <= 0)
	{
		cin >> N; if(N <= 0){cout << "Try again. Value should be >0\nN = ";}
	} 	NN = N;
	
	cout << "\n\tThe name of output file will be " << FileName << endl;
	
	ofstream File;
	File.open(FileName);
	
	while(N > 0)
	{
		W = fRand(V[0], V[1]); 
		Q2 = fRand(V[2], V[3]); 
		theta = fRand(0, M_PI); 
		phi = fRand(0, 2*M_PI);		
		ang1 = fRand(0, M_PI); 
		ang2 = fRand(0, 2*M_PI); S = -1.0;
		
		S = Section(W, Q2, theta, phi, E0);	
		
		z = 0.5*fRand(-100, 100); //Target's length 100 cm		
		x = Radius_c*2; y = Radius_c*2;
		while(x*x + y*y > Radius_c*Radius_c)
		{
			x = fRand(-Radius_c, Radius_c);
			y = fRand(-Radius_c, Radius_c);
		}
		
		if((100*N)%NN == 0)
		{
			cout << 100*(1-double(N)/double(NN)) << "%" << endl;
		}
		
		Epi = (W*W + mpi*mpi - mp*mp)/(2*W);
		Ep = (W*W + mp*mp - mpi*mpi)/(2*W);
		p = sqrt(Epi*Epi - mpi*mpi);
		nu =  (W*W + Q2 - Mp*Mp)/(2*Mp);

		meson.SetPxPyPzE(p*cos(phi)*sin(theta), p*sin(phi)*sin(theta), p*cos(theta), Epi);
		adron.SetPxPyPzE(-p*cos(phi)*sin(theta), -p*sin(phi)*sin(theta), -p*cos(theta), Ep);
		gamma1.SetPxPyPzE(mpi*cos(ang2)*sin(ang1)/2, mpi*sin(ang2)*sin(ang1)/2, mpi*cos(ang1)/2, mpi/2);
		gamma2.SetPxPyPzE(-mpi*cos(ang2)*sin(ang1)/2 ,-mpi*sin(ang2)*sin(ang1)/2 ,-mpi*cos(ang1)/2 ,mpi/2);

		beta.SetXYZ(0., 0., p/Epi);
		gamma1.Boost(beta);
		gamma2.Boost(beta);
		
		gamma1.RotateY(theta);
		gamma2.RotateY(theta);
		gamma1.RotateZ(phi); 
		gamma2.RotateZ(phi);
		
		beta.SetXYZ(0., 0., sqrt(nu*nu + Q2)/(nu + mp));
		
		adron.Boost(beta);
		meson.Boost(beta);
		gamma1.Boost(beta);
		gamma2.Boost(beta);
		
      		adron.RotateY(acos((Q2 + 2*E0*nu)/(2*E0*sqrt(nu*nu + Q2))));
		meson.RotateY(acos((Q2 + 2*E0*nu)/(2*E0*sqrt(nu*nu + Q2))));
		gamma1.RotateY(acos((Q2 + 2*E0*nu)/(2*E0*sqrt(nu*nu + Q2))));
		gamma2.RotateY(acos((Q2 + 2*E0*nu)/(2*E0*sqrt(nu*nu + Q2))));
		
		e.SetPxPyPzE((E0 - nu)*sqrt(1 - pow(1 - Q2/(2*E0*(E0 - nu)),2)), 0, (E0 - nu)*(1 - Q2/(2*E0*(E0 - nu))), E0 - nu);
		
		ang2 = fRand(0, 2*M_PI);
		
		adron.RotateZ(ang2);
		meson.RotateZ(ang2);
		gamma1.RotateZ(ang2);
		gamma2.RotateZ(ang2);
		e.RotateZ(ang2);
		
		if(channel == 0)
		{
	File << "4 0 0 0 " << polarization << " 11 " << E0 << " 2212 0 " << S << endl;
	File << "0 0 0 11 0 0 " << e.Px() << " " << e.Py() << " " << e.Pz() << " " << e.E() << " 0.0005 " << x << " " << y << " " << z << endl;
	File << "0 0 0 2212 0 0 " <<  adron.Px() << " " << adron.Py() << " " << adron.Pz() << " " << adron.E() << " " << mp << " " << x << " " << y << " " << z << endl;
	File << "0 0 0 22 0 0 " << gamma1.Px() << " " << gamma1.Py() << " " << gamma1.Pz() << " " << gamma1.E() << " 0 " << x << " " << y << " " << z << endl;
	File << "0 0 0 22 0 0 " << gamma2.Px() << " " << gamma2.Py() << " " << gamma2.Pz() << " " << gamma2.E() << " 0 " << x << " " << y << " " << z << endl;
		} else 
		{
	File << "3 0 0 0 " << polarization << " 11 " << E0 << " 2212 0 " << S << endl;
	File << "0 0 0 11 0 0 " << e.Px() << " " << e.Py() << " " << e.Pz() << " " << e.E() << " 0.0005 " << x << " " << y << " " << z << endl;
	File << "0 0 0 2112 0 0 " <<  adron.Px() << " " << adron.Py() << " " << adron.Pz() << " " << adron.E() << " " << mp << " " << x << " " << y << " " << z << endl;
	File << "0 0 0 211 0 0 " << meson.Px() << " " << meson.Py() << " " << meson.Pz() << " " << meson.E() << " " << mpi << " " << x << " " << y << " " << z << endl;
		}

		N--;
	}
	
	cout << "\tCompleted 100%" << endl;
	File.close();	
}
