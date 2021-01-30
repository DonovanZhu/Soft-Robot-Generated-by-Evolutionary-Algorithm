#define _TRAVELMAN_H_
#include<vector>
#include<cmath>

#define POPULATION_SIZE 50
#define POINT_NUM 1000
#define GENERATION 1000000

using namespace std;

class Individual
{
public:
	friend class Population;

	Individual() :Fitness(0) {}
	Individual(vector<int> genome, double fitness) : Genome(genome), Fitness(fitness) {}
	~Individual() {}

	vector<int> Genome;

	double Fitness;

	bool operator < (const Individual& I)const
	{
		return Fitness < I.Fitness;
	}
};

class Population
{
public:
	void Initial_Population(double x_coord[], double y_coord[]);
	void Crossover(double x_coord[], double y_coord[]);
	void Selection();
	void Mutation(double x_coord[], double y_coord[]);
	void Adaptive(double x_coord[], double y_coord[]);
	void TotalFit();
	void TotalFit2();
	void CompareFit();

	double Fitness(vector<int> weight, double x_coord[], double y_coord[]);
	double totalfitness = 0.0;
	double totalfitnessInv = 0.0;
	double averagefitness;
	double pastfitness;
	double P = 0.02;
	double xlist[POINT_NUM] = { 0.0 }, ylist[POINT_NUM] = { 0.0 };

	int Generation;
	int loop;
	int pressure = 2;
	vector<double> PC;
	vector<double> PM;
	vector<Individual> Population;
	vector<Individual> NewPopulation;
	vector<Individual> BestIndividual;
};

#pragma once
