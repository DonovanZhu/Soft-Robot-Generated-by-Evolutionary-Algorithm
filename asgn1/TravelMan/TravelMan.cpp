#include<iostream>
#include<fstream>
#include<ctime>
#include<random>
#include<vector>
#include<map>
#include<algorithm>
#include<typeinfo>
#include <chrono>
#include "TravelMan.h"


using namespace std;


void Population::Initial_Population(double x_coord[], double y_coord[])
{
	vector<Individual> Genometemp;

	Genometemp.clear();
	for (int i = 0; i < POPULATION_SIZE; i++)
	{

		for (int i = 0; i < POPULATION_SIZE; i++)
		{
			Genometemp.push_back(Individual());
			for (int j = 0; j < POINT_NUM; j++)
				Genometemp[i].Genome.push_back(j);            //a list from 0 to 999

			//shuffle the list:
			unsigned seed = (unsigned)chrono::system_clock::now().time_since_epoch().count();
			shuffle(Genometemp[i].Genome.begin(), Genometemp[i].Genome.end(), default_random_engine(seed));
		}

		for (int i = 0; i < POPULATION_SIZE; ++i)
		{
			Genometemp[i].Fitness = Fitness(Genometemp[i].Genome, xlist, ylist);
		}
		sort(Genometemp.begin(), Genometemp.end());
		Population.push_back(Genometemp[POPULATION_SIZE - 1]);
		Genometemp.clear();
	}
}


double Population::Fitness(vector<int> weight,double x_coord[], double y_coord[])
{

	int I1;
	int I2;
	double fitness = 0.0;

	for (int j = 1; j < POINT_NUM; j++)
	{
		I1 = weight[j];
		I2 = weight[j - 1];
		fitness += sqrt((x_coord[I1] - x_coord[I2]) * (x_coord[I1] - x_coord[I2]) + \
			(y_coord[I1] - y_coord[I2]) * (y_coord[I1] - y_coord[I2]));
	}

	return fitness;
}

#if 0
//Truncation Selection
void Population::Selection()
{
	sort(Population.begin(), Population.end());
	NewPopulation.clear();
	NewPopulation.push_back(BestIndividual[0]);
	NewPopulation.push_back(Population[POPULATION_SIZE - 1]);
	for (int i = 0; i < (POPULATION_SIZE / 2 - 2); ++i)
	{
		NewPopulation.push_back(Population[POPULATION_SIZE - 1 - i]);
	}
}
#endif

#if 1
void Population::Selection()
{
	sort(Population.begin(), Population.end());
	NewPopulation.clear();
	double sum;
	double point;
	NewPopulation.push_back(BestIndividual[0]);
	NewPopulation.push_back(Population[POPULATION_SIZE - 1]);

	for (int i = 0; i < (POPULATION_SIZE / 2 - 2); i++)
	{
		sum = 0.0;
		point = ((double)rand() / (double)RAND_MAX) * totalfitnessInv;

		vector<Individual>::iterator it;
		for (it = Population.begin(); it != Population.end();it++)
		{
			sum += 1.0 / (*(it)).Fitness;
			if (sum >= point)
			{
				NewPopulation.push_back(*it);
				break;
			}
		}

	}
}
#endif

void Population::TotalFit()
{
	totalfitnessInv = 0.0;
	for (int i = 0; i < POPULATION_SIZE; i++)
	{
		totalfitnessInv += 1.0 / Population[i].Fitness;
	}
}

void Population::TotalFit2()
{
	totalfitness = 0.0;
	for (unsigned int i = 0; i < NewPopulation.size(); i++)
	{
		totalfitness += NewPopulation[i].Fitness;
	}
	averagefitness = totalfitness / (double)NewPopulation.size();
}

#if 1
//PMX
void Population::Crossover(double x_coord[], double y_coord[])
{
	int a = 0, b = 0, x, y, c;
	bool same = true;
	double Pc_;
	int Tank[POINT_NUM] = { 0 };
	int i = 0;



	TotalFit2();
	sort(NewPopulation.begin(), NewPopulation.end());
	Adaptive(x_coord, y_coord);

	do
	{
		NewPopulation.push_back(NewPopulation[i]);
		i++;
	} while (NewPopulation.size() != POPULATION_SIZE);

	for (int i = 0; i < POPULATION_SIZE / 4; i++)
	{
		do
		{
			while (true)
			{
				x = rand() % (POPULATION_SIZE / 2);
				Pc_ = (double)rand() / (double)RAND_MAX;
				if (Pc_ <= PC[x])
					break;
			}
			while (true)
			{
				y = rand() % (POPULATION_SIZE / 2);
				Pc_ = (double)rand() / (double)RAND_MAX;
				if (Pc_ <= PC[y])
					break;
			}
		} while (x == y);
		do
		{
			a = rand() % POINT_NUM;
			b = rand() % POINT_NUM;
		} while (a == b);

		if (a > b)
		{
			c = a;
			a = b;
			b = c;
		}
		
		for (int k = a; k <= b; k++)
			Tank[k - a] = NewPopulation[x].Genome[k];

		for (int j = a; j <= b; j++)
		{
			NewPopulation[x].Genome[j] = NewPopulation[y].Genome[j];
			NewPopulation[y].Genome[j] = Tank[j - a];
			Tank[j - a] = 0;
		}

		int ind = 0;
		while (ind < a)
		{
			for (int k = a; k <= b; k++)
			{
				if (NewPopulation[x].Genome[ind] == NewPopulation[x].Genome[k])
				{
					NewPopulation[x].Genome[ind] = NewPopulation[y].Genome[k];
					same = true;
					break;
				}
				else
					same = false;
			}
			if (!same)
				ind++;
		}

		if (b < POINT_NUM - 1)
		{
			ind = b + 1;
			while (ind < POINT_NUM)
			{
				for (int k = a; k <= b; k++)
				{
					if (NewPopulation[x].Genome[ind] == NewPopulation[x].Genome[k])
					{
						NewPopulation[x].Genome[ind] = NewPopulation[y].Genome[k];
						same = true;
						break;
					}
					else
					same = false;
				}
				if (!same)
					ind++;
			}
		}
		ind = 0;
		while (ind < a)
		{
			for (int k = a; k <= b; k++)
			{
				if (NewPopulation[y].Genome[ind] == NewPopulation[y].Genome[k])
				{
					NewPopulation[y].Genome[ind] = NewPopulation[x].Genome[k];
					same = true;
					break;
				}
				else
					same = false;
			}
			if (!same)
				ind++;
		}
		if (b < POINT_NUM - 1)
		{
			ind = b + 1;
			while (ind < POINT_NUM)
			{
				for (int k = a; k <= b; k++)
				{
					if (NewPopulation[y].Genome[ind] == NewPopulation[y].Genome[k])
					{
						NewPopulation[y].Genome[ind] = NewPopulation[x].Genome[k];
						same = true;
						break;
					}
					else
						same = false;
				}
				if (!same)
					ind++;
			}
		}
		NewPopulation[x].Fitness = Fitness(NewPopulation[x].Genome, x_coord, y_coord);
		NewPopulation[y].Fitness = Fitness(NewPopulation[y].Genome, x_coord, y_coord);
	}
	
	for (int i = 0; i < POPULATION_SIZE; i++)
	{
		Pc_ = (double)rand() / (double)RAND_MAX;
		if (Pc_ < 0.05)
		{
			reverse(NewPopulation[i].Genome.begin(), NewPopulation[i].Genome.end());
			NewPopulation[i].Fitness = Fitness(NewPopulation[i].Genome, x_coord, y_coord);
		}
	}
}
#endif

#if 0
void Population::Crossover(double x_coord[], double y_coord[])
{
	int a = 0, b = 0, x, y;
	bool same = true;
	double Pc_;
	int Tank[POINT_NUM] = { 0 };
	int i = 0;

	TotalFit2();
	sort(NewPopulation.begin(), NewPopulation.end());
	Adaptive(x_coord, y_coord);

	do
	{
		NewPopulation.push_back(NewPopulation[i]);
		i++;
	} while (NewPopulation.size() != POPULATION_SIZE);

	for (int i = 0; i < POPULATION_SIZE / 4; i++)
	{
		do
		{
			while (true)
			{
				x = rand() % (POPULATION_SIZE / 2);
				Pc_ = (double)rand() / (double)RAND_MAX;
				if (Pc_ <= PC[x])
					break;
			}
			while (true)
			{
				y = rand() % (POPULATION_SIZE / 2);
				Pc_ = (double)rand() / (double)RAND_MAX;
				if (Pc_ <= PC[y])
					break;
			}
			
		} while (x == y);

		a = rand() % POINT_NUM;

		for (int k = a; k < POINT_NUM; k++)
			Tank[k - a] = NewPopulation[x].Genome[k];

		for (int j = a; j < POINT_NUM; j++)
		{
			NewPopulation[x].Genome[j] = NewPopulation[y].Genome[j];
			NewPopulation[y].Genome[j] = Tank[j - a];
			Tank[j - a] = 0;
		}

		int ind = 0;
		while (ind < a)
		{
			for (int k = a; k < POINT_NUM; k++)
			{
				if (NewPopulation[x].Genome[ind] == NewPopulation[x].Genome[k])
				{
					NewPopulation[x].Genome[ind] = NewPopulation[y].Genome[k];
					same = true;
					break;
				}
				else
					same = false;
			}
			if (!same)
				ind++;
		}

		
		ind = 0;
		while (ind < a)
		{
			for (int k = a; k < POINT_NUM; k++)
			{
				if (NewPopulation[y].Genome[ind] == NewPopulation[y].Genome[k])
				{
					NewPopulation[y].Genome[ind] = NewPopulation[x].Genome[k];
					same = true;
					break;
				}
				else
					same = false;
			}
			if (!same)
				ind++;
		}
		
		NewPopulation[x].Fitness = Fitness(NewPopulation[x].Genome, x_coord, y_coord);
		NewPopulation[y].Fitness = Fitness(NewPopulation[y].Genome, x_coord, y_coord);
	}

	for (int i = 0; i < POPULATION_SIZE; i++)
	{
		Pc_ = (double)rand() / (double)RAND_MAX;
		if (Pc_ < 0.05)
		{
			reverse(NewPopulation[i].Genome.begin(), NewPopulation[i].Genome.end());
			NewPopulation[i].Fitness = Fitness(NewPopulation[i].Genome, x_coord, y_coord);
		}
	}
}
#endif

#if 1
void Population::Mutation(double x_coord[], double y_coord[])
{
	double Pm_;
	int c;
	int ran;

	sort(NewPopulation.begin(), NewPopulation.end());
	TotalFit2();
	Adaptive(x_coord, y_coord);

	for (int i = 0; i < POPULATION_SIZE; i++) 
	{
		for (int j = 0; j < POINT_NUM; j++)
		{
			Pm_ = (double)rand() / (double)RAND_MAX;

			ran = rand() % POINT_NUM;
			if (Pm_ <= PM[i])
			{
					c = NewPopulation[i].Genome[j];
					NewPopulation[i].Genome[j] = NewPopulation[i].Genome[ran];
					NewPopulation[i].Genome[ran] = c;
			}
		}
		NewPopulation[i].Fitness = Fitness(NewPopulation[i].Genome, x_coord, y_coord);
	}
}
#endif

#if 0
void Population::Mutation(double x_coord[], double y_coord[])
{
	double Pm_;
	int c;
	int ran;

	sort(NewPopulation.begin(), NewPopulation.end());
	TotalFit2();
	Adaptive(x_coord, y_coord);

	for (int i = 0; i < POPULATION_SIZE; i++)
	{
		for (int j = 0; j < POINT_NUM; j++)
		{
			Pm_ = (double)rand() / (double)RAND_MAX;

			ran = rand() % POINT_NUM;
			if (Pm_ <= PM[i])
			{
				if (ran < POINT_NUM - 1)
				{
					c = NewPopulation[i].Genome[ran + 1];
					NewPopulation[i].Genome[ran + 1] = NewPopulation[i].Genome[ran];
					NewPopulation[i].Genome[ran] = c;
				}
				else
				{
					c = NewPopulation[i].Genome[0];
					NewPopulation[i].Genome[0] = NewPopulation[i].Genome[ran];
					NewPopulation[i].Genome[ran] = c;
				}
			}
		}
		NewPopulation[i].Fitness = Fitness(NewPopulation[i].Genome, x_coord, y_coord);
	}
}
#endif


void Population::CompareFit()
{
	sort(NewPopulation.begin(), NewPopulation.end());
	if (BestIndividual[0].Fitness >= NewPopulation[POPULATION_SIZE - 1].Fitness)
	{
		BestIndividual[0] = NewPopulation[POPULATION_SIZE - 1];
	}

}

void Population::Adaptive(double x_coord[], double y_coord[])
{
	PC.clear();
	PM.clear();
	int x = NewPopulation.size();
	for (int i = 0; i < x; i++)
	{
		if (NewPopulation[i].Fitness > averagefitness)
		{
			PC.push_back(0.9);
			PM.push_back(0.1);
		}
		else
		{
			PC.push_back(0.9 - (0.9 - 0.6) * (averagefitness - NewPopulation[i].Fitness) / (averagefitness - NewPopulation[x - 1].Fitness));
			PM.push_back(0.1 - (0.1 - 0.001) * (averagefitness - NewPopulation[i].Fitness) / (averagefitness - NewPopulation[x - 1].Fitness));
		}
	}
}