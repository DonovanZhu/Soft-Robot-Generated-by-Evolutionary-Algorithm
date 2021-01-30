#include <iostream>
#include <fstream>  
#include <vector>
#include <chrono>
#include <random>
#include <ctime>
#include <iomanip>
using namespace std;

#define POINT_NUM 1000          //1000 random points
#define GENERATION 1000000      //the number of sampling

//declear sub function:
void read();
vector<int> Random_list();
double Fitness(vector<int> Ind, double x_coord[], double y_coord[]);
vector<int> Operator();

// global variable:
double xlist[POINT_NUM] = { 0 };
double ylist[POINT_NUM] = { 0 };
double fitnesslist[GENERATION] = { 0 };

int main()
{
	vector<int> Index;
	vector<int> BestOrder;
	

	Index.clear();
	BestOrder.clear();

	read();                   //read data from tsp.txt
	BestOrder = Operator();   //output the order of points which has best fitness
	
	/*

	ofstream outputoder("order_in_1000000_big.txt");     //import order into a txt file

	for (int i = 0; i < POINT_NUM; i++)
		outputoder << setprecision(9) << BestOrder[i] << endl;

	outputoder.close();
	*/

	ofstream outputfitness("fitness_1000000_5_big.txt");  //import fitness into a txt file

	for (int i = 0; i < GENERATION; i++)
		outputfitness << setprecision(9) << fitnesslist[i] << endl;

	outputfitness.close();
	/*
	ofstream outputpoint("point_200000_1.txt");      //import points into a txt file

	for (int i = 0; i < POINT_NUM; i++)
	{
		//make sure the precision:
		outputpoint << setprecision(9) << xlist[BestOrder[i]] << "	"; 
		outputpoint << setprecision(9) << ylist[BestOrder[i]] << endl;
	}

	outputpoint.close();
	*/
	return 0;
}

//sub function :
//read data from tsp.txt
void read()
{
	int i = 0;
	ifstream input("tsp.txt");

	double a, b;
	while (input >> a >> b) {
		xlist[i] = a;
		ylist[i] = b;
		i++;
	}
	input.close();
}

//shuffle a list from 0 to 999 randomly
vector<int> Random_list()
{
	vector<int> randomlist;
	for (int i = 0;i < POINT_NUM; i++)
		randomlist.push_back(i);            //a list from 0 to 999
	
	//shuffle the list:
	unsigned seed = (unsigned)chrono::system_clock::now().time_since_epoch().count();
	shuffle(randomlist.begin(), randomlist.end(), default_random_engine(seed));
	
	return randomlist;
}

//calculate the fitness of different order
double Fitness(vector<int> Ind, double x_coord[], double y_coord[])
{
	double fitness = 0;
	double X_distance;
	double Y_distance;
	
	//the fitness is the distance between each point in order
	for (int i = 1; i < POINT_NUM; i++)
	{
		X_distance = double(x_coord[Ind[i]]) - double(x_coord[Ind[i - 1]]);
		Y_distance = double(y_coord[Ind[i]]) - double(y_coord[Ind[i - 1]]);
		fitness += sqrt(X_distance * X_distance + Y_distance * Y_distance);
	}
	X_distance = double(x_coord[Ind[0]]) - double(x_coord[Ind[Ind.size() - 1]]);
	Y_distance = double(y_coord[Ind[0]]) - double(y_coord[Ind[Ind.size() - 1]]);
	fitness += sqrt(X_distance * X_distance + Y_distance * Y_distance);
	return fitness;
}

//choose the best order and its fitness
vector<int> Operator()
{
	double fitness;
	double BestFitness;
	vector<int> Ind, BestInd;

	Ind = Random_list();
	BestFitness = Fitness(Ind, xlist, ylist);
	BestInd = Ind;

	for (int i = 0; i < GENERATION; i++)
	{
		Ind = Random_list();
		fitness = Fitness(Ind, xlist, ylist);

		if (fitness >= BestFitness)
		{
			BestFitness = fitness;
			BestInd.assign(Ind.begin(), Ind.end());
		}
		fitnesslist[i] = BestFitness;
		
	}
	return BestInd;
}