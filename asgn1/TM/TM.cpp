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

#if 0
vector<int> Population::Order(vector<double> list)
{
	vector<int> Order;
	map<int, double> iMap;

	Order.clear();
	iMap.clear();

	for (int i = 0; i < POINT_NUM; i++)
	{
		iMap[i] = list[i];
	}

	vector<pair<int, double>> vtMap;
	vtMap.clear();

	for (auto it = iMap.begin(); it != iMap.end(); it++)
		vtMap.push_back(make_pair(it->first, it->second));

	sort(vtMap.begin(), vtMap.end(),
		[](const pair<int, double>& x, const pair<int, double>& y) -> int {
			return x.second < y.second;
		});

	for (auto it = vtMap.begin(); it != vtMap.end(); it++)
		Order.push_back(it->first);
	return Order;
}
#endif

#if 0
void Population::Initial_Population(double x_coord[], double y_coord[])
{
	srand(int(time(0)));
	vector<Individual> Genometemp;
	vector<int> WeightOrder;

	Genometemp.clear();
	WeightOrder.clear();
	for (int i = 0; i < POPULATION_SIZE; i++)
	{

		for (int i = 0; i < POPULATION_SIZE; i++)
		{
			vector<double> weight;

			for (int k = 0; k < POINT_NUM; k++)
			{
				weight.push_back(((double)rand() / (double)RAND_MAX) * 5.0);
			}

			WeightOrder = Population::Order(weight);

			for (int j = 0; j < POINT_NUM; j++)
			{
				Genometemp.push_back(Individual());
				Genometemp[i].Genome.push_back(weight[j]);
			}
			weight.clear();
		}

		for (int i = 0; i < POPULATION_SIZE; ++i)
		{
			Genometemp[i].Fitness = Fitness(Genometemp[i].Genome, xlist, ylist);
		}
		sort(Genometemp.begin(), Genometemp.end());
		Population.push_back(Genometemp[POPULATION_SIZE - 1]);
		Genometemp.clear();
		WeightOrder.clear();
	}
	cout << "Initial Done" << endl;
}

#endif

#if 1
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
				Genometemp[i].Genome.push_back(j);

			unsigned seed = (unsigned)chrono::system_clock::now().time_since_epoch().count();
			shuffle(Genometemp[i].Genome.begin(), Genometemp[i].Genome.end(), default_random_engine(seed));
		}

		for (int i = 0; i < POPULATION_SIZE; ++i)
		{
			Genometemp[i].Fitness = Fitness(Genometemp[i].Genome, x_coord, y_coord);
		}
		sort(Genometemp.begin(), Genometemp.end());
		Population.push_back(Genometemp[POPULATION_SIZE - 1]);
		Genometemp.clear();
	}
	std::cout << "Initial Done" << endl;
}

#endif

double Population::Fitness(vector<int> weight, double x_coord[], double y_coord[])
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

#if 1
void Population::Selection()
{
	if (loop % 10 == 0 && loop > 1)
	{
		if (pastfitness == BestIndividual[0].Fitness)
		{
			pressure -= 2;
			if (pressure < 0)
				pressure = 2;
		}
		else
		{
			pressure += 2;
			if (pressure > 4)
				pressure = 2;
		}
	}
	std::cout << pressure << endl;
	sort(Population.begin(), Population.end());
	if (pressure == 2)
	{
		NewPopulation.clear();
		for (int i = 0; i < (POPULATION_SIZE / 2); ++i)
		{
			NewPopulation.push_back(Population[POPULATION_SIZE - 1 - i]);
		}
		std::cout << "Selection Done" << endl;
	}

	else if (pressure == 4)
	{
		NewPopulation.clear();

		for (int i = 0; i < (POPULATION_SIZE / 2); ++i)
		{
			NewPopulation.push_back(Population[POPULATION_SIZE - 1 - i]);
		}
		std::cout << "Selection Done" << endl;
	}
	else
	{
		NewPopulation.clear();

		for (int i = 0; i < (POPULATION_SIZE); ++i)
		{
			NewPopulation.push_back(Population[POPULATION_SIZE - 1 - i]);
		}
		std::cout << "Selection Done" << endl;
	}
#endif

#if 0
	void Population::Selection()
	{
		NewPopulation.clear();
		double sum;
		double point;
		sort(Population.begin(), Population.end());
		NewPopulation.push_back(Population[POPULATION_SIZE - 1]);
		NewPopulation.push_back(Population[POPULATION_SIZE - 2]);
		int z = NewPopulation.size();

		while (z != POPULATION_SIZE / 2)
		{
			sum = 0.0;
			point = ((double)rand() / (double)RAND_MAX) * totalfitnessInv;

			vector<Individual>::iterator it;
			for (it = Population.begin(); it != Population.end(); it++)
			{
				sum += (1.0 / (*(it)).Fitness);
				if (sum >= point)
				{
					NewPopulation.push_back(*it);
					totalfitnessInv -= 1.0 / (*(it)).Fitness;
					Population.erase(it);
					break;
				}
			}
			z = NewPopulation.size();
		}
		cout << "Selection Done" << endl;
	}
#endif
#if 0
	void Population::SelectionLowPresure()
	{

		NewPopulation.clear();
		sort(Population.begin(), Population.end());
		for (int i = 0; i < (POPULATION_SIZE / 2); ++i)
		{
			NewPopulation.push_back(Population[POPULATION_SIZE - 1 - i]);
		}
		cout << "Selection Done" << endl;
	}
#endif
#if 0
	void Population::SelectionHighPresure()
	{
		NewPopulation.clear();
		sort(Population.begin(), Population.end());
		for (int i = 0; i < (POPULATION_SIZE / 4); ++i)
		{
			NewPopulation.push_back(Population[POPULATION_SIZE - 1 - i]);
		}
		cout << "Selection Done" << endl;
	}
#endif
#if 0
	void Population::Selection()
	{
		NewPopulation.clear();
		sort(Population.begin(), Population.end());
		for (int i = 0; i < (POPULATION_SIZE / 8); ++i)
		{
			NewPopulation.push_back(Population[POPULATION_SIZE - 1 - i]);
		}
		for (int i = 0; i < (POPULATION_SIZE / 8); ++i)
		{
			NewPopulation.push_back(NewPopulation[i]);
		}
		cout << "Selection Done" << endl;
	}
#endif
#if 0
	void Population::TotalFit2()
	{
		totalfitnessInv = 0.0;
		totalfitness = 0.0;
		for (int i = 0; i < POPULATION_SIZE; i++)
		{
			totalfitnessInv += 1.0 / Population[i].Fitness;
			totalfitness += Population[i].Fitness;
		}
		averagefitness = totalfitness / (double)POPULATION_SIZE;
	}
#endif
	void Population::TotalFit()
	{
		totalfitness = 0.0;
		for (int i = 0; i < POPULATION_SIZE; i++)
		{
			totalfitness += Population[i].Fitness;
		}
		averagefitness = totalfitness / (double)POPULATION_SIZE;
	}

	void Population::CrossoverLow(double x_coord[], double y_coord[])
	{
		int a = 0, b = 0, x, y, c;
		bool same = true;
		double Pc_, Pc_1;
		vector<int> Tank;
		int len = NewPopulation.size();
		int z = len;

		while (NewPopulation.size() != POPULATION_SIZE)
		{
			Tank.clear();
			do
			{
				x = rand() % len;
				y = rand() % len;
			} while (x == y);
			NewPopulation.push_back(NewPopulation[x]); //NewPopulation[z]
			NewPopulation.push_back(NewPopulation[y]); //NewPopulation[z + 1]

			Pc_ = (double)rand() / (double)RAND_MAX;
			Pc_1 = (double)rand() / (double)RAND_MAX;
			if (Pc_1 < 0.2)
				reverse(NewPopulation[z].Genome.begin(), NewPopulation[z].Genome.end());

			if (Pc_ < 0.2)
				reverse(NewPopulation[z + 1].Genome.begin(), NewPopulation[z + 1].Genome.end());

			if (Pc_ < Pc)
			{
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

				int ind = 0;
				while (ind < POINT_NUM)
				{
					for (int k = a; k <= b; k++)
					{
						if (NewPopulation[y].Genome[ind] == NewPopulation[z].Genome[k])
						{
							same = true;
							break;
						}
						else
							same = false;
					}
					if (!same)
					{
						Tank.push_back(NewPopulation[y].Genome[ind]);
					}
					ind++;
				}


				ind = 0;
				if (a > 0)
					for (int i = 0; i < a; i++)
					{
						NewPopulation[z].Genome[i] = Tank[ind];
						ind++;
					}
				if (b < POINT_NUM - 1)
					for (int i = b + 1; i < POINT_NUM; i++)
					{
						NewPopulation[z].Genome[i] = Tank[ind];
						ind++;
					}

				Tank.clear();

				ind = 0;
				while (ind < POINT_NUM)
				{
					for (int k = a; k <= b; k++)
					{
						if (NewPopulation[x].Genome[ind] == NewPopulation[z + 1].Genome[k])
						{
							same = true;
							break;
						}
						else
							same = false;
					}
					if (!same)
					{
						Tank.push_back(NewPopulation[x].Genome[ind]);
					}
					ind++;
				}

				ind = 0;
				if (a > 0)
					for (int i = 0; i < a; i++)
					{
						NewPopulation[z + 1].Genome[i] = Tank[ind];
						ind++;
					}
				if (b < POINT_NUM - 1)
					for (int i = b + 1; i < POINT_NUM; i++)
					{
						NewPopulation[z + 1].Genome[i] = Tank[ind];
						ind++;
					}
				NewPopulation[z].Fitness = Fitness(NewPopulation[z].Genome, x_coord, y_coord);
				NewPopulation[z + 1].Fitness = Fitness(NewPopulation[z + 1].Genome, x_coord, y_coord);
				z += 2;
			}
		}
		cout << "Crossover Done" << endl;
	}
#if 0
	//PMX
	void Population::CrossoverHigh(double x_coord[], double y_coord[])
	{
		int a = 0, b = 0, x, y, c;
		bool same = true;
		double Pc_, Pc_1;
		int Tank[POINT_NUM] = { 0 };
		int len = NewPopulation.size();
		int z = len;

		while (NewPopulation.size() != POPULATION_SIZE)
		{
			do
			{
				x = rand() % len;
				y = rand() % len;
			} while (x == y);
			NewPopulation.push_back(NewPopulation[x]); //NewPopulation[z]
			NewPopulation.push_back(NewPopulation[y]); //NewPopulation[z + 1]


			Pc_ = (double)rand() / (double)RAND_MAX;
			Pc_1 = (double)rand() / (double)RAND_MAX;
			if (Pc_1 < 0.3)
				reverse(NewPopulation[z].Genome.begin(), NewPopulation[z].Genome.end());

			if (Pc_ < 0.3)
				reverse(NewPopulation[z + 1].Genome.begin(), NewPopulation[z + 1].Genome.end());

			if (Pc_ < Pc)
			{
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
					Tank[k - a] = NewPopulation[z].Genome[k];


				for (int j = a; j <= b; j++)
				{
					NewPopulation[z].Genome[j] = NewPopulation[z + 1].Genome[j];
					NewPopulation[z + 1].Genome[j] = Tank[j - a];
					Tank[j - a] = 0;
				}

				int ind = 0;
				while (ind < a)
				{
					for (int k = a; k <= b; k++)
					{
						if (NewPopulation[z].Genome[ind] == NewPopulation[z].Genome[k])
						{
							NewPopulation[z].Genome[ind] = NewPopulation[z + 1].Genome[k];
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
							if (NewPopulation[z].Genome[ind] == NewPopulation[z].Genome[k])
							{
								NewPopulation[z].Genome[ind] = NewPopulation[z + 1].Genome[k];
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
						if (NewPopulation[z + 1].Genome[ind] == NewPopulation[z + 1].Genome[k])
						{
							NewPopulation[z + 1].Genome[ind] = NewPopulation[z].Genome[k];
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
							if (NewPopulation[z + 1].Genome[ind] == NewPopulation[z + 1].Genome[k])
							{
								NewPopulation[z + 1].Genome[ind] = NewPopulation[z].Genome[k];
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
				NewPopulation[z].Fitness = Fitness(NewPopulation[z].Genome, x_coord, y_coord);
				NewPopulation[z + 1].Fitness = Fitness(NewPopulation[z + 1].Genome, x_coord, y_coord);
				z += 2;
			}
		}
		cout << "Crossover Done" << endl;
	}
#endif

#if 0
	void Population::CrossoverHigh(double x_coord[], double y_coord[])
	{
		int a = 0, b = 0, x, y, c;
		bool same = true;
		double Pc_, Pc_1;
		int Tank[POINT_NUM] = { 0 };
		int len = NewPopulation.size();
		int z = len;

		while (NewPopulation.size() != POPULATION_SIZE)
		{
			do
			{
				x = rand() % len;
				y = rand() % len;
			} while (x == y);
			NewPopulation.push_back(NewPopulation[x]); //NewPopulation[z]
			NewPopulation.push_back(NewPopulation[y]); //NewPopulation[z + 1]


			Pc_ = (double)rand() / (double)RAND_MAX;
			Pc_1 = (double)rand() / (double)RAND_MAX;
			if (Pc_1 < 0.1)
				reverse(NewPopulation[z].Genome.begin(), NewPopulation[z].Genome.end());

			if (Pc_ < 0.1)
				reverse(NewPopulation[z + 1].Genome.begin(), NewPopulation[z + 1].Genome.end());

			if (Pc_ < Pc)
			{
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
					Tank[k - a] = NewPopulation[z].Genome[k];


				for (int j = a; j <= b; j++)
				{
					NewPopulation[z].Genome[j] = NewPopulation[z + 1].Genome[j];
					NewPopulation[z + 1].Genome[j] = Tank[j - a];
					Tank[j - a] = 0;
				}

				int ind = 0;
				while (ind < a)
				{
					for (int k = a; k <= b; k++)
					{
						if (NewPopulation[z].Genome[ind] == NewPopulation[z].Genome[k])
						{
							NewPopulation[z].Genome[ind] = NewPopulation[z + 1].Genome[k];
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
							if (NewPopulation[z].Genome[ind] == NewPopulation[z].Genome[k])
							{
								NewPopulation[z].Genome[ind] = NewPopulation[z + 1].Genome[k];
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
						if (NewPopulation[z + 1].Genome[ind] == NewPopulation[z + 1].Genome[k])
						{
							NewPopulation[z + 1].Genome[ind] = NewPopulation[z].Genome[k];
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
							if (NewPopulation[z + 1].Genome[ind] == NewPopulation[z + 1].Genome[k])
							{
								NewPopulation[z + 1].Genome[ind] = NewPopulation[z].Genome[k];
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
				NewPopulation[z].Fitness = Fitness(NewPopulation[z].Genome, x_coord, y_coord);
				NewPopulation[z + 1].Fitness = Fitness(NewPopulation[z + 1].Genome, x_coord, y_coord);
				z += 2;
			}
		}
		cout << "Crossover Done" << endl;
	}
#endif

#if 0
	void Population::Crossover(double x_coord[], double y_coord[])
	{
		int a = 0, b = 0, c, parent;
		bool flag = true;
		srand(int(time(0)));
		double Pc_;
		double Tank[1000] = { 0 };


		for (int i = 0; i < POPULATION_SIZE / 2; i++)
		{
			//cout << NewPopulation.size() << endl;
			NewPopulation.push_back(NewPopulation[i]);
			Pc_ = (double)rand() / (double)RAND_MAX;
			if (Pc_ < Pc)
			{
				if (flag)
				{
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
					flag = !flag;

					for (int k = a; k < b; k++)
					{

						Tank[k - a] = NewPopulation[i].Genome[k];
					}
					parent = i;
				}
				else
				{
					for (int j = a; j < b; j++)
					{
						NewPopulation[parent].Genome[j] = NewPopulation[i].Genome[j];
						NewPopulation[i].Genome[j] = Tank[j - a];
						Tank[j - a] = 0;
					}

					flag = true;

					//Change the Fitness
					//cout << NewPopulation[parent].Fitness << "1" << endl;
					NewPopulation[parent].Fitness = Fitness(NewPopulation[parent].Genome, x_coord, y_coord);
					//cout << NewPopulation[parent].Fitness << "2" << endl;
					//cout << NewPopulation[i].Fitness << "1" << endl;
					NewPopulation[i].Fitness = Fitness(NewPopulation[i].Genome, x_coord, y_coord);
					//cout << NewPopulation[i].Fitness << "2" << endl;
				}
			}
		}
		cout << "Crossover Done" << endl;
	}

#endif

#if 1
	void Population::Mutation(double x_coord[], double y_coord[])
	{
		double Pm_, P;
		int c, x;
		int ran;
		sort(NewPopulation.begin(), NewPopulation.end());
		x = 0;
		for (int i = 1; i < POPULATION_SIZE; i++)
		{
			if (NewPopulation[i].Fitness == NewPopulation[i - 1].Fitness)
				x++;
			if (pressure < 2)
			{
				P = 0.1;
			}
			else if (pressure == 2)
			{
				P = 0.02;
			}
			else
			{
				P = 0.01;
			}
		}
		cout << P << endl;

		for (int i = 0; i < POPULATION_SIZE; i++)
		{
			int j = 0;
			for (; j < POINT_NUM; j++)
			{
				Pm_ = (double)rand() / (double)RAND_MAX;


				if (Pm_ <= P)
				{
					do
					{
						ran = rand() % POINT_NUM;
					} while (sqrt((x_coord[NewPopulation[i].Genome[ran]] - x_coord[NewPopulation[i].Genome[j]]) * \
						(x_coord[NewPopulation[i].Genome[ran]] - x_coord[NewPopulation[i].Genome[j]]) + \
						(y_coord[NewPopulation[i].Genome[ran]] - y_coord[NewPopulation[i].Genome[j]]) * \
						(y_coord[NewPopulation[i].Genome[ran]] - y_coord[NewPopulation[i].Genome[j]])) > 0.15);
					c = NewPopulation[i].Genome[ran];
					NewPopulation[i].Genome[ran] = NewPopulation[i].Genome[j];
					NewPopulation[i].Genome[j] = c;
					//				if (ran == POINT_NUM - 1)
					//				{
					//					c = NewPopulation[i].Genome[ran];
					//					NewPopulation[i].Genome[ran] = NewPopulation[i].Genome[0];
					//					NewPopulation[i].Genome[0] = c;
					//				}
					//				else
					//				{
					//					c = NewPopulation[i].Genome[ran];
					//					NewPopulation[i].Genome[ran] = NewPopulation[i].Genome[ran + 1];
					//					NewPopulation[i].Genome[ran + 1] = c;
					//				}
				}
			}
			NewPopulation[i].Fitness = Fitness(NewPopulation[i].Genome, x_coord, y_coord);
		}
		cout << "Mutation Done" << endl;
	}
#endif

#if 0
	void Population::Mutation()
	{
		double Pm_;
		int c;
		int i = 0;

		if (loop % 50 == 0 && loop > 1)
		{
			if (pastfitness == BestIndividual[0].Fitness)
			{
				P += 0.01;
				if (P >= 0.15)
				{
					P = 0.02;
				}
			}
			else
			{
				P -= 0.04;
				if (P <= 0.02)
				{
					P = 0.02;
				}
			}
		}
		for (i = 0; i < POPULATION_SIZE; i++)
		{
			for (int j = 0; j < POINT_NUM; j++)
			{
				Pm_ = (double)rand() / (double)RAND_MAX;
				if (Pm_ <= P)
				{
					int ran = rand() % POINT_NUM;
					c = NewPopulation[i].Genome[j];
					NewPopulation[i].Genome[j] = NewPopulation[i].Genome[ran];
					NewPopulation[i].Genome[ran] = c;
				}
			}
			NewPopulation[i].Fitness = Fitness(NewPopulation[i].Genome, xlist, ylist);
		}
		cout << P << endl;
		cout << "Mutation Done" << endl;
	}
#endif

#if 0
	void Population::Mutation()
	{
		double Pm_;
		int c;
		int i = 0;
		for (i = 0; i < POPULATION_SIZE; i++)
		{

			for (int j = 0; j < POINT_NUM; j++)
			{
				Pm_ = (double)rand() / (double)RAND_MAX;

				if (BestIndividual[0].Fitness >= 460.0)
				{
					if (Pm_ <= 0.02)
					{
						int ran = rand() % POINT_NUM;
						c = NewPopulation[i].Genome[j];
						NewPopulation[i].Genome[j] = NewPopulation[i].Genome[ran];
						NewPopulation[i].Genome[ran] = c;
					}
				}

				else if (BestIndividual[0].Fitness < 460.0 && BestIndividual[0].Fitness >= 360.0)
				{
					if (Pm_ <= 0.02)
					{
						int ran = rand() % POINT_NUM;
						c = NewPopulation[i].Genome[j];
						NewPopulation[i].Genome[j] = NewPopulation[i].Genome[ran];
						NewPopulation[i].Genome[ran] = c;
					}
				}
#if 1
				else if (BestIndividual[0].Fitness < 360.0 && BestIndividual[0].Fitness >= 260.0)
				{

					if (Pm_ <= 0.07)
					{
						int ran = rand() % POINT_NUM;
						c = NewPopulation[i].Genome[j];
						NewPopulation[i].Genome[j] = NewPopulation[i].Genome[ran];
						NewPopulation[i].Genome[ran] = c;
					}
				}
#endif
				else
				{

					if (Pm_ <= 0.1)
					{
						int ran = rand() % POINT_NUM;
						c = NewPopulation[i].Genome[j];
						NewPopulation[i].Genome[j] = NewPopulation[i].Genome[ran];
						NewPopulation[i].Genome[ran] = c;
					}
				}
			}
			NewPopulation[i].Fitness = Fitness(NewPopulation[i].Genome, xlist, ylist);
		}
		cout << "Mutation Done" << endl;
	}
#endif


	void Population::CompareFit()
	{
		sort(NewPopulation.begin(), NewPopulation.end());
		if (BestIndividual[0].Fitness >= NewPopulation[POPULATION_SIZE - 1].Fitness)
		{
			BestIndividual[0] = NewPopulation[POPULATION_SIZE - 1];
		}

		Population[0] = BestIndividual[0];

	}