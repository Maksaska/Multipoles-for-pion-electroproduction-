#include <iostream> 
#include "header.h"

using namespace std;

const double Mp(0.93827), Mn(0.93957), Mpip(0.13957), Mpiz(0.13498);

int main(int argc, char **argv)
{	
	vector<double> Area; 
	vector<bool> Contributions;
	vector<int> Parameters;
	int Purpose(0);
	
	Greetings(Area, Contributions, Parameters, Purpose); 
	
	if(Contributions.size() == 0){return 0;}
	if(Parameters.size() == 0){return 0;}
	
	if(Purpose == 0)
	{
		vector<vector<double>> Multipoles;
		
		Initialization(Multipoles, Area, Parameters[0]); 
		
		if(Contributions[0])
		{
			vector<vector<double>> Born;
			Born_table(Born, Area, Parameters[0]);
			Merge(Multipoles, Born);
			Born.clear();		
		}

		if(Contributions[1])
		{
			vector<vector<double>> P33; 
			P33_table(P33, Area, Parameters[0]);
			Merge(Multipoles, P33);
			P33.clear();		
		}
		
		if(Contributions[2])
		{
			vector<vector<double>> P11; 
			P11_table(P11, Area, Parameters[0]);
			Merge(Multipoles, P11);
			P11.clear();		
		}
		
		PrintCSV(Multipoles, Parameters);		
		Multipoles.clear();	
	}
	
	if(Purpose == 1)
	{
		Generate_lund(Area, Parameters[0]);
	}
	
	Area.clear(); Contributions.clear(); Parameters.clear(); 
   	
	return 0;
}
