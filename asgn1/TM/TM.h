#ifndef _TRAVELMAN_H_
#define _TRAVELMAN_H_
#include<vector>
#include<cmath>

#define POPULATION_SIZE 100
#define POINT_NUM 1000
#define GENERATION 10000
#define Pc 0.8
#define Pm 0.2

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
		return Fitness > I.Fitness;
	}
};

class Population
{
public:
	void Initial_Population(double x_coord[], double y_coord[]);
	vector<int> Order(vector<double> list);
	//void Seqencing(std::vector<Individual>& vecpop);
	double Fitness(vector<int> weight, double x_coord[], double y_coord[]);

	void CrossoverLow(double x_coord[], double y_coord[]);
	void CrossoverHigh(double x_coord[], double y_coord[]);
	void Selection();

	void Mutation(double x_coord[], double y_coord[]);
	void TotalFit();
	void TotalFit2();
	void CompareFit();

	int Generation;
	int loop = 0;
	double totalfitness = 0.0;
	double totalfitnessInv = 0.0;
	double averagefitness;
	double pastfitness;
	double P = 0.02;
	int pressure = 2;

	double xlist[POINT_NUM] = { 0 }, ylist[POINT_NUM] = { 0 };

	vector<Individual> Population;
	vector<Individual> NewPopulation;
	vector<Individual> BestIndividual;
};
#endif
#pragma once
