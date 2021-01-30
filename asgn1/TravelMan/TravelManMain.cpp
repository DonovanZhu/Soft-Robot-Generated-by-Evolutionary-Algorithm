#include "TravelMan.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <algorithm>
using namespace std;


int main()
{
	Population p;
	int cha[GENERATION];
	static double pathlength[GENERATION] = { 0.0 };
	int i= 0; 
	srand(int(time(0)));

	ifstream input("tsp.txt");

	double a , b;
	while (input >> a >> b) {
		p.xlist[i] = a;
		p.ylist[i] = b;
		i++;
	}
	input.close();
	
	p.Initial_Population(p.xlist, p.ylist);

	p.BestIndividual.push_back(p.Population[0]);
	p.pastfitness = p.BestIndividual[0].Fitness;

	srand(int(time(0)));
	for ( p.loop = 0; p.loop < GENERATION; p.loop++)
	{

		double fit = p.BestIndividual[0].Fitness;
		if (p.loop % 200 == 0 && p.loop != 0)
		{
			if (p.pastfitness == p.BestIndividual[0].Fitness)
			{
				p.Initial_Population(p.xlist, p.ylist);
				cout << "Catastrophic" << endl;
			}
			p.pastfitness = p.BestIndividual[0].Fitness;
		}

		p.TotalFit();
		p.Selection();
		p.Crossover(p.xlist, p.ylist);
		p.Mutation(p.xlist, p.ylist);
		p.Population.clear();
		sort(p.NewPopulation.begin(), p.NewPopulation.end());

		for (int s = 0; s < POPULATION_SIZE; s++)
			p.Population.push_back(p.NewPopulation[s]);		

		p.CompareFit();
		if (p.BestIndividual[0].Fitness < fit)
		{
			cout << "Generation" << p.loop + 1 << ":" << p.BestIndividual[0].Fitness << endl;
		}
		pathlength[p.loop] = p.BestIndividual[0].Fitness;

		for (int i = 0; i < POPULATION_SIZE; i++)
		{
			int x = 0;
			if (p.NewPopulation[i].Fitness < 90.0)
			{
				x++;
			}
			cha[p.loop] = x;
		}
	}
	/*
	ofstream outputorder("Order_R2A.txt"); 

	for (int i = 0; i < POINT_NUM; i++)
		outputorder << p.xlist[p.BestIndividual[0].Genome[i]]<< "	" << p.ylist[p.BestIndividual[0].Genome[i]] << endl;

	outputorder.close();
	*/
	ofstream outputorder("chara.txt");

	for (int i = 0; i < GENERATION; i++)
		outputorder << cha[i] << endl;

	outputorder.close();
	/*
	ofstream outputfitness("fitness_R2A.txt");

	for (int i = 0; i < GENERATION; i++)
	{
		outputfitness << setprecision(9) << pathlength[i] << endl;
	}

	outputfitness.close();
	*/
	p.NewPopulation.clear();
	p.Population.clear();
	return 0;
}