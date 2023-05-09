#include <iostream>
#include <cmath>
#include <random>
#include <tuple>
#include <fstream>

//SYSTEM CONSTANTS
const int N = 75000000; //number of monte carlo steps
// const int nskip = 15; //
const int sitesx = 100; //number of sites along x
const int sitesy = 100; //number of sites along y
const int q = 3; //possible values of spins
const double kT = 0.4; //boltzmann constant times temperature
const double J = 1.0; //coupling constant

//random number generator
std::random_device rd;
std::uniform_int_distribution<int> dist(0,q-1);
std::uniform_int_distribution<int> sitex(0,sitesx-1);
std::uniform_int_distribution<int> sitey(0,sitesy-1);
std::uniform_real_distribution<float> w(0.0,1.0);

//global array to store our 2D system
int sites[sitesx][sitesy] = {0};

//initialize our system with random spins given by q
void init(){
	for(int i=0;i<sitesx;i++){
		for(int j=0;j<sitesy;j++){
			sites[i][j] = dist(rd);
		}		
	}
}

//the weight to be used for the Metropolis Monte-Carlo algorithm
double weight(double dH){
    return exp(-dH/kT);
}

//returns 1 if the are the same and 0 if they are not
bool kroneckerDelta(int s1, int s2){
	return (s1==s2);
}

std::tuple<bool,int> metropolis(int x, int y){
	//to ensure PBCs
	int s0 = sites[x][y];
	int sn = dist(rd);
	int l = (x-1)*(x!=0)+(sitesx-1)*(x==0);//left
	int r = (x+1)*(x!=sitesx-1);//right
	int d = (y-1)*(y!=0)+(sitesy-1)*(y==0);//down
	int u = (y+1)*(y!=sitesy-1);//up
	double oldH = -J*(kroneckerDelta(s0,sites[l][y])+kroneckerDelta(s0,sites[r][y])+kroneckerDelta(s0,sites[x][u])+kroneckerDelta(s0,sites[x][d]));
	//std::cout << oldH <<std::endl;
	double newH = -J*(kroneckerDelta(sn,sites[l][y])+kroneckerDelta(sn,sites[r][y])+kroneckerDelta(sn,sites[x][u])+kroneckerDelta(sn,sites[x][d]));
	//std::cout << newH <<std::endl;
	float we = w(rd);
	//std::cout << weight(newH-oldH) << " " << we <<std::endl;
	double delH = newH-oldH;
	if(delH <= 0){
		return std::make_pair(true,sn);
	}else if(weight(delH)>=we){
		return std::make_pair(true,sn);
	}
	return std::make_pair(false,s0);
}

int main(){
	init();
	int count = 0;
	for(int i=0;i<N;i++){
		int x = sitex(rd);
		int y = sitey(rd);
		//std::cout << x << " " << y << std::endl;
		auto[change,sn] = metropolis(x,y);
		//std::cout << change << " " << sn << std::endl;
		sites[x][y] = sn;
		//return 0;
		if(change){
			count++;
		}
		if(i==0||i==12500000-1||i==N-1){
			std::cout << std::endl << i << std::endl;
			for (int i=0; i<sitesx; i++){
		        for (int j=0; j<sitesy; j++){
		        	std::cout << sites[i][j] << " ";
		    	}
		    	std::cout << std::endl;
    		}
		}
	}
	
	
}