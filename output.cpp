#include <iostream> 
#include <fstream>
#include "header.h"

using namespace std; 

const double Mpip(0.13957), Mp(0.93827);

void Reading(string Path,vector<vector<double>>& V2)
{
	string line; stringstream ss;
	ifstream File;
	double dub; vector<double> Numbers;

	File.open(Path,fstream::in | fstream::out | fstream::app);

	if (!File.is_open())
	{
		cout << "Can't open " << Path << " !" << endl;
	}
	else
	{	
		while(!File.eof())
		{		
			getline(File,line); 	
			ss << line;

			if(File.eof())
			{
				break;
			}

			while(ss>>dub)
			{							
				Numbers.push_back(dub);
				while(ss.peek() == ',')
            			ss.ignore();
			}
	
			V2.push_back(Numbers);			
			Numbers.clear();
			ss.clear();
		}
	}
	File.close();
}

void Initialization(vector<vector<double>>& V, vector<double>& WQ2, int l)
{
	vector<double> buff; 
	
	for(double Q2 = WQ2[3]; Q2 <= WQ2[4]; Q2 = Q2 + WQ2[5])
	{
		for(double W = WQ2[0]; W <= WQ2[1]; W = W + WQ2[2])
		{
			buff.push_back(W);
			buff.push_back(Q2);
			for(int i = 0; i < 24*(l+1); i++)
			{
				buff.push_back(0);
			}
			
			V.push_back(buff);			
			buff.clear();
		}
	} 
}

void Merge(vector<vector<double>>& V, vector<vector<double>>& A) 
{
	for (long unsigned int i = 0; i < V.size(); i++)
	{
		for(long unsigned int j = 2; j < V[i].size(); j++)
		{
			V[i][j] += A[i][j];  
		}
	}
}

void PrintCSV(vector<vector<double>>& V, vector<int>& P) 
{
	double Const(1), C_L(1); long unsigned int j_1, j_2;
	
	string Sh = "All multipoles are in 1/GeV";
	
	if(P[2] == 1)
	{
		Const = Mpip*1000;
		Sh = "All multipoles are in 10^-3/MPi+";
	}
	
	ofstream File;
	File.open("Table.csv");
	
	File << Sh << endl;
	File << endl; 
	
	if(P[1] == 0){ j_1 = 1 + V[0].size()/2; j_2 = V[0].size(); File << "Pi0p" << endl;}
	else if(P[1] == 1){ j_1 = 2; j_2 = 1 + V[0].size()/2; File << "Pi+n" << endl;}
	else if(P[1] == 2){ j_1 = 1 + V[0].size()/2; j_2 = V[0].size(); File << "Pi0p" << endl;}
	
	File << "W_GeV,Q2_GeV^2,E0+_Re,0_Im,";
	
	for(int i = 1; i < (P[0]+1); i++)
	{
		File << i << "_Re," << i << "_Im,";
	}
	
	File << "E0-_Re,0_Im" << ","; 
	
	for(int i = 1; i < (P[0]+1); i++)
	{
		File << i << "_Re," << i << "_Im,";
	}
	
	File << "M0+_Re,0_Im" << ","; 
	
	for(int i = 1; i < (P[0]+1); i++)
	{
		File << i << "_Re," << i << "_Im,";
	}
	
	File << "M0-_Re,0_Im" << ","; 
	
	for(int i = 1; i < (P[0]+1); i++)
	{
		File << i << "_Re," << i << "_Im,";
	}
	
	if(P[3] == 0){File << "S0+_Re,0_Im" << ",";}
	else{File << "L0+_Re,0_Im" << ",";}									
	
	for(int i = 1; i < (P[0]+1); i++)
	{
		File << i << "_Re," << i << "_Im,";
	}
	
	if(P[3] == 0){File << "S0-_Re,0_Im" << ",";}
	else{File << "L0-_Re,0_Im" << ",";}	
	
	for(int i = 1; i < (P[0]+1); i++)
	{
		File << i << "_Re," << i << "_Im,";
	}	
	File << endl;
	
	for(long unsigned int i = 0; i < V.size(); i++)
	{
		File << V[i][0] << "," << V[i][1] << ","; 
		for(long unsigned int j = j_1; j < j_2; j++)
		{
		if(j >= j_1 + 8*(P[0]+1) and P[3] == 1){C_L = (V[i][0]*V[i][0] - Mp*Mp - V[i][1])/(2*V[i][0]*k_mod(V[i][0],V[i][1]));}
			File << C_L*Const*V[i][j] << ","; C_L = 1;
		}
		File << endl;	
	} 
	
	if(P[1] == 2)
	{
		j_1 = 2; j_2 = 1 + V[0].size()/2; File << "\nPi+n" << endl;
		
		for(long unsigned int i = 0; i < V.size(); i++)
		{
			File << V[i][0] << "," << V[i][1] << ","; 
			for(long unsigned int j = j_1; j < j_2; j++)
			{
			if(j >= j_1 + 8*(P[0]+1) and P[3] == 1){C_L = (V[i][0]*V[i][0] - Mp*Mp - V[i][1])/(2*V[i][0]*k_mod(V[i][0],V[i][1]));}
				File << C_L*Const*V[i][j] << ","; C_L = 1;
			}
			File << endl;	
		}
	}
	
	File.close();
}

void PrintBiggy(vector<vector<double>>& V)
{
	for (long unsigned int i = 0; i < V.size(); i++)
	{
		for(long unsigned int j = 0; j < V[i].size(); j++)
		{
			cout << V[i][j] << "\t";
		}
		cout << endl;
	}
}

bool Response()
{
	string answer;
	while(true)
	{
		cin >> answer;
		if(answer == "1"){ return true;}
		else if(answer == "positive"){ return true;}
		else if(answer == "Positive"){ return true;}
		else if(answer == "y"){ return true;}
		else if(answer == "Y"){ return true;}
		else if(answer == "yes"){ return true;}
		else if(answer == "Yes"){ return true;}
		else if(answer == "YES"){ return true;}
		else if(answer == "0"){ return false;}
		else if(answer == "no"){ return false;}
		else if(answer == "No"){ return false;}
		else if(answer == "n"){ return false;}
		else if(answer == "N"){ return false;}
		else if(answer == "NO"){ return false;}
		else if(answer == "negative"){ return false;}
		else if(answer == "Negative"){ return false;}
		else
		{
			cout << "\n\nCan't recognize your response. Please, use Yes/No words.\n Your answer:";
			
		}	
	}
}

void Greetings(vector<double>& V, vector<bool>& P1, vector<int>& P2, int &B) 	
{
	double buff(0), buff_2; bool value; int l(-1);  
	cout << " ------------------------------------------------------------------- " << endl;
	cout << "| Welcome to the multipole evaluation program for the pion       | \n| electroproduction!                                       |       \n|                                                                   |\n|     Authors: Davydov M. - MSU, Physics dep.                       |\n|              Isupov E.  - MSU, SINP                               |\n|                                                   Version 1.0    |\n| https://github.com/Maksaska/none |\n ------------------------------------------------------------------- " << endl;
	
	cout << endl;
	
	cout << "What do you want?\nCSV table with multipoles:";
	value = Response();
	
	if(!value)
	{
		cout << "PiN event generator:";
		value = Response(); B = 1;
	}
	
	if(!value){return;}
	
	cout << endl;
	
	cout << "Enter the kinematic area of (W, Q2) in GeV and GeV^2 respectively.\n\nW_min = ";
	
	while(buff < 1.1)
	{
		cin >> buff; if(buff < 1.1){cout << "\nW should be >= 1.1 GeV. Try again.\nW_min = ";} 	
	} V.push_back(buff); buff_2 = buff;
	
	cout << "W_max = "; buff = 0;
	
	while(buff <= buff_2)
	{
		cin >> buff; if(buff < buff_2){cout << "\nYour W_max < W_min. Try again.\nW_max = ";}	
	} V.push_back(buff); buff = 0;
	
	if(B == 0)
	{
		cout << "Increment of W: "; 
	
		while(buff <= 0)
		{
			cin >> buff; if(buff <= 0){cout << "\nIncorrect value. Try again.\nIncrement of W: ";} 	
		} V.push_back(buff); cout << endl; buff = 0;
	}	

	cout << "Q2_min = ";	
	
	while(buff <= 0)
	{
		cin >> buff; if(buff <= 0){cout << "\nNon-physical value. Try again.\nQ2_min = ";} 	
	} V.push_back(buff); buff_2 = buff; buff = 0;
	
	cout << "Q2_max = ";
	
	while(buff <= buff_2)
	{
		cin >> buff; if(buff < buff_2){cout << "\nYour Q2_max < Q2_min. Try again.\nQ2_max = ";}	
	} V.push_back(buff); buff = 0;
	
	if(B == 0)
	{
		cout << "Increment of Q2: ";
		
		while(buff <= 0)
		{
			cin >> buff; if(buff <= 0){cout << "\nIncorrect value. Try again.\nIncrement of Q2: ";} 	
		} V.push_back(buff); cout << endl; buff = 0;
	}
	
	if(B == 0)
	{
		cout << "\n";
		
		cout << "What are we going to include?\nBorn terms:";
		value = Response(); P1.push_back(value); if(value){l++;}
		
		cout << "P33(1232):"; 
		value = Response(); P1.push_back(value); if(value){l++;}
		
		cout << "P11(1440):"; 
		value = Response(); P1.push_back(value); if(value){l++;}	
		
		if(l == -1){return;} 
		
		l = -1;	
	}	
	
	if(B == 0)
	{
		cout << "\n\nEnter the max value of orbital momentum.\nl = ";
	
		while(l < 0)
		{
			cin >> l; if(l < 0){ cout << "\nYou need at least one member of a row (l >= 0). Try again.\nl = ";} 	
		} P2.push_back(l); cout << endl; l = -1;
		
		cout << "Choose the channel:\t\t(Pi0p - 0, Pi+n - 1, Both - 2)\nAnswer:";
		
		while(l != 0 and l != 1 and l != 2)
		{
	cin >>l; if(l != 0 and l != 1 and l != 2){ cout << "\nCan't recognize your response. Try again.\t\t(Pi0p - 0, Pi+n - 1, Both - 2)\n Answer:";}
		} P2.push_back(l); cout << endl; l = -1;
		
		cout << "Select multipole units: 1/GeV - 0 or 10^-3/MPi+ - 1\nAnswer:";

		while(l != 0 and l != 1)
		{
	cin >>l; if(l != 0 and l != 1){ cout << "\nCan't recognize your response. Try again.\t\t(1/GeV - 0 or 10^-3/MPi+ - 1)\nAnswer:";}
		} P2.push_back(l); cout << endl; l = -1;
		
		cout << "Select multipole set: EMS - 0 or EML - 1\nAnswer:";
		
		while(l != 0 and l != 1)
		{
	cin >>l; if(l != 0 and l != 1){ cout << "\nCan't recognize your response. Try again.\t\t(EMS - 0 or EML - 1)\nAnswer:";}
		} P2.push_back(l); cout << endl;
	} 

	if(B == 1)
	{	
		cout << "Choose the channel:\t\t(Pi0p - 0, Pi+n - 1)\nAnswer:";
		
		while(l != 0 and l != 1)
		{
	cin >>l; if(l != 0 and l != 1){ cout << "\nCan't recognize your response. Try again.\t\t(Pi0p - 0, Pi+n - 1)\n Answer:";}
		} P2.push_back(l); cout << endl; l = -1; P1.push_back(true);
	} 	
	
}
