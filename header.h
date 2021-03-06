#include <vector>
#include <iostream> 
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <ctime>
#include "TLorentzVector.h"

using namespace std;

constexpr std::complex<double> operator""_i(unsigned long long d)
{
    return std::complex<double>{0.0, static_cast<double>(d)};
}
constexpr std::complex<double> operator""_i(long double d)
{
    return std::complex<double>{0.0, static_cast<double>(d)};
}

/*Output file*/

void Reading(string Path,vector<vector<double>>& V2); // Reading function for csv double grid
void PrintBiggy(vector<vector<double>>& V); //Output function
bool Response(); // "yes" or "no" response recognition algorithm
void Greetings(vector<double>& V, vector<bool>& P1, vector<int>& P2, int &B);	// Greetings and input parameters
void Initialization(vector<vector<double>>& V, vector<double>& WQ2, int l); // Initialization of the final table
void Merge(vector<vector<double>>& V, vector<vector<double>>& A); // Combining the final table and the selected contribution 
void PrintCSV(vector<vector<double>>& V, vector<int>& P); // Output procedure

/*general.cpp*/

double fRand(const double& fMin, const double& fMax); // Random function with max resolution
double q(const double& W); // Pion momentum
double q2p(const double& W); // Double pion momentum
double qeta(const double& W); // Eta momentum
double Gamma_PiN(const double& W, vector<double>& P, const int& l); // partial width res -> PiN decay
double Gamma_in(const double& W, vector<double>& P, const int& l); // partial width res -> 2PiN decay
double Gamma_eta(const double& W, vector<double>& P, const int& l); // partial width res -> etaN decay
double Gamma_tot(const double& W, vector<double>& P, vector<int>& I); // total partial width
double k(const double& W); // Photon momentum; Q2 = 0
double k_mod(const double& W, const double& Q2); // Photon momentum; Q2 != 0
double f_PiN(const double& W, vector<double>& P, vector<int>& I, const double& Q2); // Breit-Wigner factor describing the decay of the N*resonance
double f_gN(const double& W, const double& WR, const int& n, const double& Q2); //Parametrization of the yNN* vertex beyond the resonance peak

/*P33 resonance*/

double AB(const double& W,vector<string>& MES); // Q2 parametrization for a phase dependence
double phi_delta(const double& W, const double& Q2, const string& mult); // phase function 
vector<complex<double>> mult_delta(const double& W, const double& Q2); // Multipoles evaluation 
void P33_table(vector<vector<double>>& Result, vector<double>& WQ2, int l_max); // This function provide with the table of multipoles for given (W, Q2) grid (P33 - resonance)

/*P11 resonance*/

double phi_P11(const double& W, const double& Q2); // phase function 
complex<double> mult_P11(const double& W, const double& Q2); // Multipoles evaluation pM1-(1/2)
void P11_table(vector<vector<double>>& Result, vector<double>& WQ2, int l_max); // This function provide with the table of multipoles for given (W, Q2) grid (P11 - resonance)

/*Born terms*/

double F1p_f(const double& Q2); //Form-factors
double F2p_f(const double& Q2);
double F1n_f(const double& Q2);
double F2n_f(const double& Q2);
double Fpi_f(const double& Q2);
double F1v_f(const double& Q2);
double F2v_f(const double& Q2);
double F1s_f(const double& Q2);
double F2s_f(const double& Q2);

double Q(const int& n, const double& x); 				       	// Legendre polinomials of second kind (Berends 1967)
vector<double> Gamma(const double& Q2, const double& t, const int& index);		// Gamma functions for multipoles (Born termes)
int ksi(const int& index);						        	// ksi function for multipoles (Born terms)
int Parity(const int& m); 						        	// Parity function (-1)^l
double RN(const int& l, const double& E2_barred);    					// RN expression for multipoles (Born terms)
double RPi(const int& l, const double& q0_barred);                     		// RPi expression for multipoles (Born terms)
double T(const int& l, const double& Q2, const double& W, const int& pm_index);	// T expression for multipoles (Born terms)
int delta_f(const int& q1,const int& q2);						// Kronecker symbol

vector<vector<double>> EMS(const double& W, const double& Q2, const double& theta, const int& limiter, const int& id);
// This function calculates all 6 multipoles for Born terms; limiter provide a number of members in a row; id - isospin of multipole (+, -, 0)

void Born_table(vector<vector<double>>& Result, vector<double>& WQ2, int l_max);
// This function provide with the table of multipoles for given (W, Q2) grid (Born terms)

/* Generator.cpp*/

double P(const int& der,const int& n, double& theta); //legendre polynomials P(cos)n
vector<complex<double>> Merger(const double& W,  const double& Q2); // Multipoles compound
vector<complex<double>> Helicity_amplitudes(const double& W,  const double& Q2, const double& theta); // Helicity amplitudes evaluation for chosen contributions
double Section(const double& W,  const double& Q2, const double& theta,  const double& phi, const double& E0); // diff section evaluation
void Generate_lund(vector<double>& V, int chanel); // Main sequence for generator
