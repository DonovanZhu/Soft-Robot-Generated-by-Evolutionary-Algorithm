#include "GP.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <algorithm>
#include <string>
using namespace std;


int main()
{
	Population p;
	static double pathlength[GENERATION] = { 0.0 };
	int i = 0;
	srand(int(time(0)));

	ifstream input("data.txt");

	double a, b;
	static double x[1000] = { 0 }, y[1000] = { 0 };
	while (input >> a >> b) {
		x[i] = a;
		y[i] = b;
		i++;
	}
	input.close();
	int g = 0;
	for (int i = 0; i < 1000; i += 10)
	{
		p.xlist[g] = x[i];
		p.ylist[g] = y[i];
		g++;
	}


	p.Initial_Population(POPULATION_SIZE);

	p.BestIndividual.push_back(p.Population[0]);

	srand(int(time(0)));
	ofstream outputdot("conv.txt");
	for (p.loop = 0; p.loop < GENERATION; p.loop++)
	{

		double fit = p.BestIndividual[0].Fitness;
		p.TotalFit();
		p.Selection();
		outputdot << a << endl;
		p.Crossover();
		p.Mutation();
		p.Population.clear();
		sort(p.NewPopulation.begin(), p.NewPopulation.end());
		for (int s = 0; s < POPULATION_SIZE; s++)
			p.Population.push_back(p.NewPopulation[s]);
		
		p.CompareFit();
		p.Simplify();

		if (p.BestIndividual[0].Fitness < fit)
		{
			cout << "Generation" << p.loop + 1 << ":" << p.BestIndividual[0].Fitness << endl;
			cout << p.To_string(p.BestIndividual[0].Genome, 1) << endl;
		}
		pathlength[p.loop] = p.BestIndividual[0].Fitness;
		if (p.loop % 10 == 0)
		{
			cout << "Generation" << p.loop + 1 << endl;
		}
		for (int i = 0; i < POPULATION_SIZE; i++)
		{
			for (int j = 0; j < 256; j++)
			{
				if ((p.NewPopulation[i].Genome[j] == "+" || p.NewPopulation[i].Genome[j] == "-" || p.NewPopulation[i].Genome[j] == "/" 
					|| p.NewPopulation[i].Genome[j] == "*" ) && p.NewPopulation[i].Genome[2 * j] == "0" && p.NewPopulation[i].Genome[2*j+1] == "0")
				{
					if (i == 0)
					{
						p.NewPopulation[i] = p.NewPopulation[i + 1];
					}
					else
					{
						p.NewPopulation[i] = p.NewPopulation[0];

					}
				}
			}
		}

		int a = 0;
		for (int i = 0; i < POPULATION_SIZE; i++)
		{
			if (p.NewPopulation[i].Fitness < 0.3)
			{
				a++;
			}
		}

	}

	outputdot.close();
	ofstream outputfitness("fitness.txt");

	for (int i = 0; i < GENERATION; i++)
	{
		outputfitness << setprecision(9) << pathlength[i] << endl;
	}

	outputfitness.close();

	p.NewPopulation.clear();
	p.Population.clear();
	return 0;
}