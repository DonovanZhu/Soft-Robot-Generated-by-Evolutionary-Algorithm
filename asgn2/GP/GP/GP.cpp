#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <random>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>
#include <typeinfo>
#include <chrono>
#include <string>
#include "GP.h"


using namespace std;


void Population::Initial_Population(int P)
{
	vector<Individual> Genometemp;
	Genometemp.clear();
	for (int i = 0; i < P; i++)
	{
		int num;
		for (num = 0; num < POPULATION_SIZE; num++)
		{
			Genometemp.push_back(Individual());
			for (int j = 0; j < 256; j++)
				Genometemp[num].Genome.push_back("0");

			Genometemp[num].Genome[1] = operator_dic[rand() % 6];
			for (int i = 1; i < 7; i++)
			{
				for (int j = (int)pow(2, i); j < (int)pow(2, i + 1); j++)
				{
					for (int k = 0; k < 4; k++)
					{
						if (Genometemp[num].Genome[int(j / 2.0)] == calculator[k])
						{
							Genometemp[num].Genome[j] = operator_dic[rand() % 8];
							break;
						}
					}

					if (Genometemp[num].Genome[int(j / 2)] == tri[0] || Genometemp[num].Genome[int(j / 2)] == tri[1])
					{
						if (j % 2 == 0)
							Genometemp[num].Genome[j] = operator_dic[rand() % 8];
					}
				}
			}

			for (int i = (int)pow(2, 7); i < (int)pow(2, 8); i++)
			{
				if (Genometemp[num].Genome[int(i / 2)] != (string)"0" && Genometemp[num].Genome[int(i / 2)] != (string)"x" && Genometemp[num].Genome[int(i / 2)] != (string)"a")
				{
					if (Genometemp[num].Genome[int(i / 2)] == (string)"cos" || Genometemp[num].Genome[int(i / 2)] == (string)"sin")
					{
						if (i % 2 == 0)
						{
							Genometemp[num].Genome[i] = operator_dic[6 + rand() % 2];
						}
					}
					else
					{
						Genometemp[num].Genome[i] = operator_dic[6 + rand() % 2];
					}
				}
			}

			for (int i = 0; i < 256; i++)
			{
				if (Genometemp[num].Genome[i] == "a")
				{
					char buffer[20];
					sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
					string str = buffer;
					Genometemp[num].Genome[i] = str;
				}
			}

		}

		for (int i = 0; i < POPULATION_SIZE; ++i)
		{
			Genometemp[i].Fitness = Fitness(Genometemp[i].Genome);
		}
		sort(Genometemp.begin(), Genometemp.end());
		Population.push_back(Genometemp[0]);
		Genometemp.clear();
	}
}


double Population::Fitness(vector<string> gene)
{
	double error = 0.0;
	double fitness[256];
	for (int k = 0; k < POINT_NUM; k++)
	{
		for (int i = 0; i < 256; i++)
			fitness[i] = 0.0;


		for (int i = 255; i > 0; i--)
		{

			if (gene[i] == "x")
				fitness[i] = xlist[k];

			else if (gene[i] == "+" && i < 158)
				fitness[i] = fitness[2 * i] + fitness[2 * i + 1];

			else if (gene[i] == "-" && i < 158)
				fitness[i] = fitness[2 * i] - fitness[2 * i + 1];

			else if (gene[i] == "*" && i < 158)
				fitness[i] = fitness[2 * i] * fitness[2 * i + 1];

			else if (gene[i] == "/" && i < 158)
				fitness[i] = fitness[2 * i] / fitness[2 * i + 1];

			else if (gene[i] == "sin" && i < 158)
				fitness[i] = sin(fitness[2 * i]);

			else if (gene[i] == "cos" && i < 158)
				fitness[i] = cos(fitness[2 * i]);

			else if (gene[i] == "0")
				continue;
			else
				fitness[i] = atof(gene[i].c_str());
		}
		error += fabs(fitness[1] - ylist[k])/POINT_NUM;
	}
	return error;
}

#if 1
//Truncation Selection
void Population::Selection()
{
	sort(Population.begin(), Population.end());
	NewPopulation.clear();
    NewPopulation.push_back(BestIndividual[0]);
	vector<Individual>::iterator it;
	for (it = Population.begin(); it != Population.end();)
	{
		if (it - Population.begin() >= 80)
		{
			it = Population.erase(it);
		}
		else
		{
			it++;
		}
	}
	int a = Population.size();
	Initial_Population(POPULATION_SIZE - a);

	
	for (int i = 0; i < POPULATION_SIZE; i++)
	{
		if ((int)(i / 10) % 2 == 0)
		{
			NewPopulation.push_back(Population[i]);
		}
	}
	
}
#endif

#if 0
void Population::Selection()
{
	double sum;
	double point;
	sort(Population.begin(), Population.end());
	NewPopulation.clear();
//	NewPopulation.push_back(BestIndividual[0]);
	NewPopulation.push_back(Population[0]);

	while (NewPopulation.size() != POPULATION_SIZE / 2)
	{
		sum = 0.0;
		point = ((double)rand() / (double)RAND_MAX) * totalfitnessInv;

		vector<Individual>::iterator it;
		for (it = Population.begin(); it != Population.end(); it++)
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

#if 0
//PMX
void Population::Crossover()
{
	int a = 0, b = 0, x, y;
	bool same = true;
	double Pc_;
	int i = 0;
	sort(NewPopulation.begin(), NewPopulation.end());
	TotalFit2();
	Adaptive();

	do
	{
		NewPopulation.push_back(NewPopulation[i]);
		i++;
	} while (NewPopulation.size() != POPULATION_SIZE);

	for (int i = 0; i < POPULATION_SIZE / 5; i++)
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

		int node1,node2;
		do
		{
			do
			{
				node1 = rand() % 255 + 1;
			} while (NewPopulation[x].Genome[node1] == "0");
			
			do
			{
				node2 = rand() % 255 + 1;
			} while (NewPopulation[y].Genome[node2] == "0");
			

		} while ( (int)(log(node1)/log(2)) != (int)(log(node2) / log(2)));


		Cage.clear();
		for (int i = 0; i < 256; i++)
			Cage.push_back("0");

		Save(NewPopulation[x].Genome, node1);

		Clear(NewPopulation[x].Genome, node1);

		Copy(NewPopulation[x].Genome, NewPopulation[y].Genome, node1, node2);

		Clear(NewPopulation[y].Genome, node2);

		Copy(NewPopulation[y].Genome, Cage, node2, node1);

		NewPopulation[x].Fitness = Fitness(NewPopulation[x].Genome);
		NewPopulation[y].Fitness = Fitness(NewPopulation[y].Genome);
	}
	
}
#endif

#if 1
void Population::Crossover()
{
	int a = 0, b = 0, x, y;
	double Pc_;
	int i = 0;
	sort(NewPopulation.begin(), NewPopulation.end());

	do
	{
		NewPopulation.push_back(NewPopulation[i]);
		i++;
	} while (NewPopulation.size() != POPULATION_SIZE);

	vector<int> cross;
	cross.clear();
	for (int i = 0; i < POPULATION_SIZE / 2; i++)
		cross.push_back(i);
	unsigned seed = (unsigned)chrono::system_clock::now().time_since_epoch().count();
	shuffle(cross.begin(), cross.end(), default_random_engine(seed));

	for (int i = 0; i < POPULATION_SIZE / 2; i += 2)
	{
		Pc_ = (double)rand() / (double)RAND_MAX;
		x = cross[i];
		y = cross[i + 1];
		if(Pc_ < CROSSOVER)
		{

			int node1, node2;
			do
			{
				do
				{
					node1 = rand() % 255 + 1;
				} while (NewPopulation[x].Genome[node1] == "0");

				do
				{
					node2 = rand() % 255 + 1;
				} while (NewPopulation[y].Genome[node2] == "0");


			} while (IsLong(x, node1) + (int)(log(node2) / log(2)) >= 7 || IsLong(y, node2) + (int)(log(node1) / log(2)) >= 7);

			Cage.clear();
			for (int i = 0; i < 256; i++)
				Cage.push_back("0");

			Save(x, node1);

			Clear(x, node1);


			Copy(x, y, node1, node2);
			Clear(y, node2);

			Copy2(y, node2, node1);

			NewPopulation[x].Fitness = Fitness(NewPopulation[x].Genome);
			NewPopulation[y].Fitness = Fitness(NewPopulation[y].Genome);

			
			//Deterministic Crowding 
			if (Distance(NewPopulation[x].Genome, NewPopulation[x + POPULATION_SIZE / 2].Genome) + \
				Distance(NewPopulation[y].Genome, NewPopulation[y + POPULATION_SIZE / 2].Genome) < \
				Distance(NewPopulation[y].Genome, NewPopulation[x + POPULATION_SIZE / 2].Genome) + \
				Distance(NewPopulation[x].Genome, NewPopulation[y + POPULATION_SIZE / 2].Genome))
			{
				if (NewPopulation[x].Fitness < NewPopulation[x + POPULATION_SIZE / 2].Fitness)
				{
					NewPopulation[x + POPULATION_SIZE / 2] = NewPopulation[x];
				}
				if (NewPopulation[y].Fitness < NewPopulation[y + POPULATION_SIZE / 2].Fitness)
				{
					NewPopulation[y + POPULATION_SIZE / 2] = NewPopulation[y];
				}
			}
			else
			{
				if (NewPopulation[x].Fitness < NewPopulation[y + POPULATION_SIZE / 2].Fitness)
				{
					NewPopulation[y + POPULATION_SIZE / 2] = NewPopulation[x];
				}
				if (NewPopulation[y].Fitness < NewPopulation[x + POPULATION_SIZE / 2].Fitness)
				{
					NewPopulation[x + POPULATION_SIZE / 2] = NewPopulation[y];
				}
			}
			

		}
	}
}
#endif

#if 0
void Population::Mutation()
{
	double Pm_ = 0.0;

	sort(NewPopulation.begin(), NewPopulation.end());
	TotalFit2();
	Adaptive();
	int ind;
	for (ind = 0; ind < POPULATION_SIZE; ind++)
	{

		Pm_ = (double)rand() / (double)RAND_MAX;

		if (Pm_ <= PM[ind])
		{
			while (true)
			{
				int node = rand() % 256;
				if (NewPopulation[ind].Genome[node] == "0")
					continue;
				else
				{
					double ran = (double)rand() / (double)RAND_MAX;
					if (ran < 0.9)
					{
						if (NewPopulation[ind].Genome[node] == "x")
						{
							if (node < 128 && node > 0)
							{
								int r = rand() % 8;
								NewPopulation[ind].Genome[node] = operator_dic[r];
								if (r < 4)
								{
									for (int i = node * 2; i <= node * 2 + 1; i++)
									{
										double r1 = (double)rand() / (double)RAND_MAX;
										if (r1 >= 0.5)
											NewPopulation[ind].Genome[i] = "x";

										else
										{
											char buffer[20];
											sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
											string str = buffer;
											NewPopulation[ind].Genome[i] = str;
										}
									}
								}
								else if (r == 4 || r == 5)
								{
									double r1 = (double)rand() / (double)RAND_MAX;
									if (r1 >= 0.5)
										NewPopulation[ind].Genome[node * 2] = "x";

									else
									{
										char buffer[20];
										sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
										string str = buffer;
										NewPopulation[ind].Genome[node * 2] = str;
									}
								}
								else if (r == 6)
									NewPopulation[ind].Genome[node] == "x";
								else if (r == 7)
								{
									char buffer[20];
									sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
									string str = buffer;
									NewPopulation[ind].Genome[node] = str;
								}
							}
						}
						else if ((NewPopulation[ind].Genome[node] == "+" \
							|| NewPopulation[ind].Genome[node] == "-" \
							|| NewPopulation[ind].Genome[node] == "*" \
							|| NewPopulation[ind].Genome[node] == "/" \
							|| NewPopulation[ind].Genome[node] == "sin" \
							|| NewPopulation[ind].Genome[node] == "cos") && node < 128)
						{
							int r2 = rand() % 8;
							if (r2 < 4 && (NewPopulation[ind].Genome[node] == "sin" || NewPopulation[ind].Genome[node] == "cos"))
							{
								double r3 = (double)rand() / (double)RAND_MAX;
								if (r3 >= 0.5)
								{
									char buffer[20];
									sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
									string str = buffer;
									NewPopulation[ind].Genome[2 * node + 1] = str;
								}
								else
									NewPopulation[ind].Genome[2 * node + 1] = "x";
								NewPopulation[ind].Genome[node] = operator_dic[r2];
							}
							else if (r2 >= 4 && r2 <= 5 && (NewPopulation[ind].Genome[node] == "+" || \
								NewPopulation[ind].Genome[node] == "-" || \
								NewPopulation[ind].Genome[node] == "*" || \
								NewPopulation[ind].Genome[node] == "/"))
							{
								Clear(ind, 2 * node + 1);
								NewPopulation[ind].Genome[node] = operator_dic[r2];
							}

							else if (r2 == 6)
							{
								if (node > 1 && node < 256)
								{
									Clear(ind, node);
									NewPopulation[ind].Genome[node] = operator_dic[r2];
								}
							}
							else if (r2 == 7)
							{								
								if (node > 1 && node < 256)
								{
									Clear(ind, node);
									char buffer[20];
									sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
									string str = buffer;
									NewPopulation[ind].Genome[node] = str;
								}
							}
						}
						else if (NewPopulation[ind].Genome[node] == "0")
							continue;
						else
						{
							double r = (double)rand() / (double)RAND_MAX;
							if (r < 0.5)
							{
								double c = atof(NewPopulation[ind].Genome[node].c_str());
								if ((double)rand() / (double)RAND_MAX >= 0.5)
									c += (double)rand() / (double)RAND_MAX;

								else
									c -= (double)rand() / (double)RAND_MAX;

								if (c > 10.0)
									c -= 10.0;
								else if (c < -10.0)
									c += 10.0;

								char buffer[20];
								sprintf_s(buffer, "%.10f", c);
								string str = buffer;
								NewPopulation[ind].Genome[node] = str;
							}
							else
							{
								int r = rand() % 8;
								if (node < 128)
								{
									if (r < 4)
									{
										NewPopulation[ind].Genome[node] = operator_dic[r];
										for (int i = node * 2; i <= node * 2 + 1; i++)
										{
											double r1 = (double)rand() / (double)RAND_MAX;
											if (r1 >= 0.5)
												NewPopulation[ind].Genome[i] = "x";

											else
											{
												char buffer[20];
												sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
												string str = buffer;
												NewPopulation[ind].Genome[i] = str;
											}
										}
									}
									else if (r == 4 || r == 5)
									{

										NewPopulation[ind].Genome[node] = operator_dic[r];
										double r1 = (double)rand() / (double)RAND_MAX;
										if (r1 >= 0.5)
											NewPopulation[ind].Genome[node * 2] = "x";

										else
										{
											char buffer[20];
											sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
											string str = buffer;
											NewPopulation[ind].Genome[node * 2] = str;
										}
									}
									else if (r == 6)
									{
										NewPopulation[ind].Genome[node] = "x";
									}
									else if (r == 7)
										continue;
								}
							}
						}
					}
					else if (ran >= 0.9 && node > 1)
					{
						string tank = "\0";
						Clear(ind, node);
						NewPopulation[ind].Genome[node] = "x";
					}
					break;
				}
			}
		}
		NewPopulation[ind].Fitness = Fitness(NewPopulation[ind].Genome);
	}
}
#endif

#if 1
void Population::Mutation()
{
	double Pm_ = 0.0;

	sort(NewPopulation.begin(), NewPopulation.end());
//	TotalFit2();
//	Adaptive();
	int ind;
	for (ind = 0; ind < POPULATION_SIZE; ind++)
	{
		Pm_ = (double)rand() / (double)RAND_MAX;

//		if (Pm_ <= PM[ind])
		if(Pm_ <= MUTATION)
		{
			while (true)
			{
				int node = rand() % 256;
				if (NewPopulation[ind].Genome[node] == "0")
					continue;
				else
				{
					double x = (double)rand() / (double)RAND_MAX;
					if (x < 1.0)
					{
						if ((NewPopulation[ind].Genome[node] == "x" || \
							NewPopulation[ind].Genome[node] == "+" || \
							NewPopulation[ind].Genome[node] == "-" || \
							NewPopulation[ind].Genome[node] == "*" || \
							NewPopulation[ind].Genome[node] == "/" || \
							NewPopulation[ind].Genome[node] == "sin" || \
							NewPopulation[ind].Genome[node] == "cos") && node < 128)
						{
							int pin = rand() % 3;
							int len = IsLong(ind, node);
							if ((pin == 0) && (log(node) / log(2) + len < 8))
							{
								Cage.clear();
								for (int i = 0; i < 256; i++)
									Cage.push_back("0");
								Save(ind, node);
								Clear(ind, node);
								Copy2(ind, 2 * node, node);
								NewPopulation[ind].Genome[node] = "+";
								char buffer[20];
								sprintf_s(buffer, "%.10f", ((double)rand() / (double)RAND_MAX - 0.5) * 0.1);
								string str = buffer;
								NewPopulation[ind].Genome[2 * node + 1] = str;
							}
							else if ((pin == 1) && (log(node) / log(2) + len < 8))
							{
								Cage.clear();
								for (int i = 0; i < 256; i++)
									Cage.push_back("0");
								Save(ind, node);
								Clear(ind, node);
								Copy2(ind,  2 * node, node);
								NewPopulation[ind].Genome[node] = "*";
								char buffer[20];
								sprintf_s(buffer, "%.10f", ((double)rand() / (double)RAND_MAX - 0.5) * 0.1 + 1);
								string str = buffer;
								NewPopulation[ind].Genome[2 * node + 1] = str;
							}
							else if (pin == 2)
							{
								Clear(ind, node);
								char buffer[20];
								sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
								string str = buffer;
								NewPopulation[ind].Genome[node] = str;
							}
						}
						else
						{
							double x = (double)rand() / (double)RAND_MAX;
							if (x < 0.6)
							{
								double c = atof(NewPopulation[ind].Genome[node].c_str());
								if ((double)rand() / (double)RAND_MAX >= 0.5)
									c += (double)rand() / (double)RAND_MAX * 0.5;

								else
									c -= (double)rand() / (double)RAND_MAX * 0.5;

								if (c > 10.0)
									c -= 10.0;
								else if (c < -10.0)
									c += 10.0;

								char buffer[20];
								sprintf_s(buffer, "%.10f", c);
								string str = buffer;
								NewPopulation[ind].Genome[node] = str;
							}
							else
							{
								NewPopulation[ind].Genome[node] = "x";
							}

						}
					}
					else
					{
						if ((NewPopulation[ind].Genome[node] == "+" \
							|| NewPopulation[ind].Genome[node] == "-" \
							|| NewPopulation[ind].Genome[node] == "*" \
							|| NewPopulation[ind].Genome[node] == "/" \
							|| NewPopulation[ind].Genome[node] == "sin" \
							|| NewPopulation[ind].Genome[node] == "cos") && node < 128)
						{
							int r2 = rand() % 8;
							if (r2 < 4 && (NewPopulation[ind].Genome[node] == "sin" || NewPopulation[ind].Genome[node] == "cos"))
							{
								NewPopulation[ind].Genome[2 * node + 1] = "x";
								NewPopulation[ind].Genome[node] = operator_dic[r2];
							}
							else if (r2 >= 4 && r2 <= 5 && (NewPopulation[ind].Genome[node] == "+" || NewPopulation[ind].Genome[node] == "-" || NewPopulation[ind].Genome[node] == "*" || NewPopulation[ind].Genome[node] == "/")
								&& ((NewPopulation[ind].Genome[(int)(node / 2)] != "sin" && NewPopulation[ind].Genome[(int)(node / 2)] != "cos") || \
								 (NewPopulation[ind].Genome[(int)(node / 4)] != "sin" && NewPopulation[ind].Genome[(int)(node / 4)] != "cos")))
							{
								Clear(ind, 2 * node + 1);
								NewPopulation[ind].Genome[node] = operator_dic[r2];
							}

							else if (r2 == 6)
							{
								if (node > 1 && node < 256)
								{
									Clear(ind, node);
									NewPopulation[ind].Genome[node] = operator_dic[r2];
								}
							}
							/*
							else if (r2 == 7)
							{
								continue;
								
								if (node > 1 && node < 256)
								{
									Clear(ind, node);
									char buffer[20];
									sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
									string str = buffer;
									NewPopulation[ind].Genome[node] = str;
								}
								
							}
							*/
						}
					}
					break;
				}
			}
		}
		NewPopulation[ind].Fitness = Fitness(NewPopulation[ind].Genome);
	}
}
#endif


void Population::CompareFit()
{
	sort(NewPopulation.begin(), NewPopulation.end());
	if (BestIndividual[0].Fitness > NewPopulation[0].Fitness)
		BestIndividual[0] = NewPopulation[0];
}


void Population::Adaptive()
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
			PC.push_back(0.9 - (0.9 - 0.6) * (averagefitness - NewPopulation[i].Fitness) / (averagefitness - NewPopulation[0].Fitness));
			PM.push_back(0.1 - (0.1 - 0.001) * (averagefitness - NewPopulation[i].Fitness) / (averagefitness - NewPopulation[0].Fitness));
		}
	}
}

void Population::Clear(int i, int x)
{
	NewPopulation[i].Genome[x] = "0";
	if (2 * x < 256)
	{
		Clear(i, 2 * x);
		Clear(i, 2 * x + 1);
	}
}

void Population::Save(int i, int x)
{
	Cage[x] = NewPopulation[i].Genome[x];
	if (2 * x < 256)
	{
		Save(i, 2 * x);
		Save(i, 2 * x + 1);
	}
}

void Population::Copy(int x, int y, int node1, int node2)
{
	NewPopulation[x].Genome[node1] = NewPopulation[y].Genome[node2];
	if (2 * max(node1,node2) < 256)
	{
		Copy(x, y, 2 * node1, 2 * node2);
		Copy(x, y, 2 * node1 + 1, 2 * node2 + 1);
	}
}

void Population::Copy2(int x, int node1, int node2)
{
	NewPopulation[x].Genome[node1] = Cage[node2];
	if (2 * max(node1,node2) < 256)
	{
		Copy2(x, 2 * node1, 2 * node2);
		Copy2(x, 2 * node1 + 1, 2 * node2 + 1);
	}
}

string Population::To_string(vector<string> gene, int x)
{
	string function = "\0";
	if ((2 * x > 256) || (gene[2 * x] == (string)"0"))
		function += gene[x];
	else if (gene[x] == tri[0] || gene[x] == tri[1])
		function = function + gene[x] + (string)"(" + To_string(gene, 2 * x) + (string)")";
	else
		function = function + (string)"(" + To_string(gene, 2 * x) + gene[x] + To_string(gene, 2 * x + 1) + (string)")";
	return function;
}

int Population::IsLong(int Ind, int node)
{
	int i = 1;
	int len = 1;
	while (node * pow(2, i) < 256)
	{
		int j = node * (int)pow(2, i);
		int k = 1;
		while (k < (int)pow(2, i))
		{
			if (NewPopulation[Ind].Genome[j] != "0")
			{
				len++;
				break;
			}
			j++;
			k++;
		}
		i++;
	}
	return len;
}

double Population::Distance(vector<string> gene1, vector<string> gene2)
{
	double error = 0.0;
	double fitness[256] = { 0.0 };
	double fitness2[256] = { 0.0 };
	for (int k = 0; k < POINT_NUM; k += 10)
	{
		for (int i = 0; i < 256; i++)
			fitness[i] = 0.0;

		for (int i = 255; i > 0; i--)
		{

			if (gene1[i] == "x")
				fitness[i] = xlist[k];

			else if (gene1[i] == "+" && i < 158)
				fitness[i] = fitness[2 * i] + fitness[2 * i + 1];

			else if (gene1[i] == "-" && i < 158)
				fitness[i] = fitness[2 * i] - fitness[2 * i + 1];

			else if (gene1[i] == "*" && i < 158)
				fitness[i] = fitness[2 * i] * fitness[2 * i + 1];

			else if (gene1[i] == "/" && i < 158)
				fitness[i] = fitness[2 * i] / fitness[2 * i + 1];

			else if (gene1[i] == "sin" && i < 158)
				fitness[i] = sin(fitness[2 * i]);

			else if (gene1[i] == "cos" && i < 158)
				fitness[i] = cos(fitness[2 * i]);

			else if (gene1[i] == "0")
				continue;
			else
				fitness[i] = atof(gene1[i].c_str());
		}

		for (int i = 0; i < 256; i++)
			fitness2[i] = 0.0;


		for (int i = 255; i > 0; i--)
		{

			if (gene2[i] == "x")
				fitness2[i] = xlist[k];

			else if (gene2[i] == "+" && i < 158)
				fitness2[i] = fitness2[2 * i] + fitness2[2 * i + 1];

			else if (gene2[i] == "-" && i < 158)
				fitness2[i] = fitness2[2 * i] - fitness2[2 * i + 1];

			else if (gene2[i] == "*" && i < 158)
				fitness2[i] = fitness2[2 * i] * fitness2[2 * i + 1];

			else if (gene2[i] == "/" && i < 158)
				fitness2[i] = fitness2[2 * i] / fitness2[2 * i + 1];

			else if (gene2[i] == "sin" && i < 158)
				fitness2[i] = sin(fitness2[2 * i]);

			else if (gene2[i] == "cos" && i < 158)
				fitness2[i] = cos(fitness2[2 * i]);

			else if (gene2[i] == "0")
				continue;
			else
				fitness2[i] = atof(gene2[i].c_str());
		}

		error += fabs(fitness[1] - fitness2[1]) / 10;
	}
	return error;
}

void Population::Simplify()
{
	if (BestIndividual[0].Fitness < 0.05)
	{
		for (int i = 0; i < 256; i++)
		{
			if (BestIndividual[0].Genome[i] == "+" && BestIndividual[0].Genome[2 * i].size() > 8 && BestIndividual[0].Genome[2 * i + 1].size() > 8)
			{
				double a = atof(BestIndividual[0].Genome[2 * i].c_str());
				double b = atof(BestIndividual[0].Genome[2 * i + 1].c_str());
				BestIndividual[0].Genome[2 * i] = "0";
				BestIndividual[0].Genome[2 * i + 1] = "0";
				char buffer[20];
				sprintf_s(buffer, "%.10f", a + b);
				string str = buffer;
				BestIndividual[0].Genome[i] = str;
			}
			else if (BestIndividual[0].Genome[i] == "-" && BestIndividual[0].Genome[2 * i].size() > 8 && BestIndividual[0].Genome[2 * i + 1].size() > 8)
			{
				double a = atof(BestIndividual[0].Genome[2 * i].c_str());
				double b = atof(BestIndividual[0].Genome[2 * i + 1].c_str());
				BestIndividual[0].Genome[2 * i] = "0";
				BestIndividual[0].Genome[2 * i + 1] = "0";
				char buffer[20];
				sprintf_s(buffer, "%.10f", a - b);
				string str = buffer;
				BestIndividual[0].Genome[i] = str;
			}
			else if (BestIndividual[0].Genome[i] == "*" && BestIndividual[0].Genome[2 * i].size() > 8 && BestIndividual[0].Genome[2 * i + 1].size() > 8)
			{
				double a = atof(BestIndividual[0].Genome[2 * i].c_str());
				double b = atof(BestIndividual[0].Genome[2 * i + 1].c_str());
				BestIndividual[0].Genome[2 * i] = "0";
				BestIndividual[0].Genome[2 * i + 1] = "0";
				char buffer[20];
				sprintf_s(buffer, "%.10f", a * b);
				string str = buffer;
				BestIndividual[0].Genome[i] = str;
			}
			else if (BestIndividual[0].Genome[i] == "/" && BestIndividual[0].Genome[2 * i].size() > 8 && BestIndividual[0].Genome[2 * i + 1].size() > 8)
			{
				double a = atof(BestIndividual[0].Genome[2 * i].c_str());
				double b = atof(BestIndividual[0].Genome[2 * i + 1].c_str());
				BestIndividual[0].Genome[2 * i] = "0";
				BestIndividual[0].Genome[2 * i + 1] = "0";
				char buffer[20];
				sprintf_s(buffer, "%.10f", a/b);
				string str = buffer;
				BestIndividual[0].Genome[i] = str;
			}
			else if (BestIndividual[0].Genome[i] == "sin" && BestIndividual[0].Genome[2 * i].size() > 8)
			{
				double a = atof(BestIndividual[0].Genome[2 * i].c_str());
				BestIndividual[0].Genome[2 * i] = "0";
				char buffer[20];
				sprintf_s(buffer, "%.10f", sin(a));
				string str = buffer;
				BestIndividual[0].Genome[i] = str;
			}
			else if (BestIndividual[0].Genome[i] == "cos" && BestIndividual[0].Genome[2 * i].size() > 8)
			{
				double a = atof(BestIndividual[0].Genome[2 * i].c_str());
				BestIndividual[0].Genome[2 * i] = "0";
				char buffer[20];
				sprintf_s(buffer, "%.10f", cos(a));
				string str = buffer;
				BestIndividual[0].Genome[i] = str;
			}
			else if (BestIndividual[0].Genome[i] == "-" && (BestIndividual[0].Genome[2 * i] == BestIndividual[0].Genome[2 * i + 1]))
			{
				BestIndividual[0].Genome[2 * i] = "0";
				BestIndividual[0].Genome[2 * i + 1] = "0";
				BestIndividual[0].Genome[i] = "0.000000000";
			}
			else if (BestIndividual[0].Genome[i] == "/" && (BestIndividual[0].Genome[2 * i] == BestIndividual[0].Genome[2 * i + 1]))
			{
				BestIndividual[0].Genome[2 * i] = "0";
				BestIndividual[0].Genome[2 * i + 1] = "0";
				BestIndividual[0].Genome[i] = "1.000000000";
			}
		}
	}
}
