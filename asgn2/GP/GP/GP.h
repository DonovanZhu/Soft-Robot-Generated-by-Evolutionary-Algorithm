#define _GP_H_
#include<vector>
#include<cmath>
#include<string>

#define POPULATION_SIZE 100
#define POINT_NUM 100
#define GENERATION 1000
#define MUTATION 0.2
#define CROSSOVER 0.7

using namespace std;

class Individual
{
public:
	friend class Population;

	Individual() :Fitness(0) {}
	Individual(vector<string> genome, double fitness) : Genome(genome), Fitness(fitness) {}
	~Individual() {}

	vector<string> Genome;

	double Fitness;

	bool operator < (const Individual& I)const
	{
		return Fitness < I.Fitness;
	}
};

class Population
{
public:
	void Initial_Population(int P);
	void Crossover();
	void Selection();
	void Mutation();
	void Adaptive();
	void TotalFit();
	void TotalFit2();
	void CompareFit();
	void Clear(int i, int x);
	void Save(int i, int x);
	void Copy(int x, int y, int node1, int node2);
	void Copy2(int x, int node1, int node2);
	int IsLong(int Ind, int node);
	string To_string(vector<string> gene, int x);
	double Distance(vector<string> gene1, vector<string> gene2);
	void Simplify();
	//double Similar(vector<string> gene1, vector<string> gene2);

	double Fitness(vector<string> gene);
	double totalfitness = 0.0;
	double totalfitnessInv = 0.0;
	double averagefitness;
	double pastfitness;
	double P = 0.02;
	double xlist[POINT_NUM] = { 0.0 }, ylist[POINT_NUM] = { 0.0 };

	int Generation;
	int loop;
	vector<double> PC;
	vector<double> PM;
	vector<string> Cage;
	vector<Individual> Population;
	vector<Individual> NewPopulation;
	vector<Individual> BestIndividual;
	string operator_dic[8] = { "+","-","*","/","sin","cos","x","a" };
	string calculator[4] = { "+","-","*","/" };
	string tri[2] = { "sin","cos" };
	string cons[2] = { "x","a" };
};



