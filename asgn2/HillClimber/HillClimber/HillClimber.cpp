#include <iostream>
#include <fstream>  
#include <vector>
#include <chrono>
#include <random>
#include <ctime>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdio>

using namespace std;

#define POINT_NUM 1000
#define GENERATION 100000

void Read();
void Generate();
string To_string(string gene[], int x);
void Clear(string gene[], int x);
double Operate(string gene[]);
void mutation();

static string operator_dic[8] = { "+","-","*","/","sin","cos","x","a" };
static string calculator[4] = { "+","-","*","/" };
static string tri[2] = { "sin","cos" };
static string cons[2] = { "x","a" };

static double xlist[POINT_NUM] = { 0 };
static double ylist[POINT_NUM] = { 0 };
//static double xtest[POINT_NUM] = { 0 };
//static double ytest[POINT_NUM] = { 0 };
static double BestYlist[POINT_NUM] = { 0 };
static double fitnesslist[GENERATION] = { 0 };

static string tree[256];
static string test[256];

int main()
{
	double fitness = 10000.0;
	double bestfitness = fitness;
	srand(int(time(0)));

	Read();

	bool valid;
	do
	{
		valid = true;
		for (int i = 0; i < 256; i++)
			tree[i] = "0";
		Generate();
		for (int i = 0; i < sizeof(tree); i++)
		{
			if (i < (256 / 2 - 1) && tree[i] == (string)"/" && tree[2 * i + 1] == (string)"0")
			{
				valid = false;
				break;
			}
			if (i < (256 / 2 - 1) && tree[i] == "*" && (tree[2 * i] == "0" || tree[2 * i + 1] == "0"))
			{
				valid = false;
				break;
			}
		}
	} while (!valid);


	for (int loop = 0; loop < GENERATION; loop++)
	{

		mutation();
		double f1 = Operate(tree);
		double f2 = Operate(test);
		if (f1 > f2)
		{
			for (int i = 0; i < 256; i++)
				tree[i] = test[i];
		}
		fitness = Operate(tree);

		string func = To_string(tree, 1);

		if (bestfitness > fitness)
		{
			bestfitness = fitness;
			cout << "--------------------" << endl;
			cout << loop << endl;
			cout << func << endl;
			cout << fitness << endl;
		}
		fitnesslist[loop] = bestfitness;
	}
	ofstream outputfitness("HillClimber5.txt");  //import fitness into a txt file

	for (int i = 0; i < GENERATION; i++)
		outputfitness << setprecision(9) << fitnesslist[i] << endl;

	outputfitness.close();
}

void Read()
{
	int i = 0;
	ifstream input("data.txt");
	double x[1000] = { 0 };
	double y[1000] = { 0 };
	double a, b;
	while (input >> a >> b)
	{
		x[i] = a;
		y[i] = b;
		i++;
	}
	input.close();

	for (int i = 0; i < POINT_NUM; i++)
	{
		xlist[i] = x[i];
		ylist[i] = y[i];
	}

	//	for (int i = 1; i < POINT_NUM; i += 2)
	//	{
	//		xtest[i] = x[i];
	//		ytest[i] = y[i];
	//	}
}

void Generate()
{
	bool flag = false;
	tree[1] = operator_dic[rand() % 6];
	for (int i = 1; i < 7; i++)
	{
		for (int j = (int)pow(2, i); j < (int)pow(2, i + 1); j++)
		{
			for (int k = 0; k < 4; k++)
			{
				if (tree[int(j / 2.0)] == calculator[k])
				{
					tree[j] = operator_dic[rand() % 8];
					flag = true;
					break;
				}
			}

			if (tree[int(j / 2)] == tri[0] || tree[int(j / 2)] == tri[1])
			{
				if (j % 2 == 0)
					tree[j] = operator_dic[rand() % 8];
			}
		}
	}

	for (int i = (int)pow(2, 7); i < (int)pow(2, 8); i++)
	{
		if (tree[int(i / 2)] != (string)"0" && tree[int(i / 2)] != (string)"x" && tree[int(i / 2)] != (string)"a")
		{
			if (tree[int(i / 2)] == (string)"cos" || tree[int(i / 2)] == (string)"sin")
			{
				if (i % 2 == 0)
				{
					tree[i] = operator_dic[6 + rand() % 2];
				}
			}
			else
			{
				tree[i] = operator_dic[6 + rand() % 2];
			}
		}
	}

	for (int i = 0; i < 256; i++)
	{
		if (tree[i] == "a")
		{
			char buffer[20];
			sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
			string str = buffer;
			tree[i] = str;
		}
	}
	//	cout << "Generate Done" << endl;
}


string To_string(string gene[], int x)
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

double Operate(string gene[])
{

	double error = 0.0;
	double fitness[256];
	for (int k = 0; k < POINT_NUM; k++)
	{
		int j = 0;
		for (int i = 0; i < 256; i++)
		{
			fitness[i] = 0.0;
		}

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
			{
				fitness[i] = atof(gene[i].c_str());
			}

		}
		error += fabs(fitness[1] - ylist[k]);
	}
	//	cout << "Operate Done" << endl;
	return error;
}

void mutation()
{

	int node = rand() % 256;

	for (int i = 0; i < 256; i++)
		test[i] = tree[i];

	while (true)
	{
		if (test[node] == "0")
		{
			node = rand() % 256;
			continue;
		}
		else
		{
			double ran = (double)rand() / (double)RAND_MAX;
			if (ran < 0.8)
			{
				if (test[node] == "x")
				{
					if (node < 128)
					{
						int r = rand() % 8;
						test[node] = operator_dic[r];
						if (r < 4)
						{
							//							for (int i = node * 2; i <= node * 2 + 1; i++)
							//							{
							//								double r1 = (double)rand() / (double)RAND_MAX;
							//								if (r1 >= 0.5)
							test[node * 2] = "x";

							//								else
							//								{
							char buffer[20];
							sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
							string str = buffer;
							test[node * 2 + 1] = str;
							//								}
							//							}

						}
						else if (r == 4 || r == 5)
						{
							//							if (test[node / 2] == "sin" || test[node / 2] == "cos")
							//							{
							//								test[node] = "x";
							//							}
							//							else
							//							{
							double r1 = (double)rand() / (double)RAND_MAX;
							if (r1 >= 0.5)
								test[node * 2] = "x";

							else
							{
								char buffer[20];
								sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
								string str = buffer;
								test[node * 2] = str;
							}
							//							}

						}
						else if (r == 6)
							continue;
						else if (r == 7)
						{
							char buffer[20];
							sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
							string str = buffer;
							test[node] = str;
						}

						//						cout << "Mutate X" << endl;
					}
				}
				else if ((test[node] == "+" || test[node] == "-" || test[node] == "*" || test[node] == "/" || test[node] == "sin" || test[node] == "cos") && node < 128)
				{
					int r2 = rand() % 8;
					if (r2 < 4 && (test[node] == "sin" || test[node] == "cos"))
					{
						double r3 = (double)rand() / (double)RAND_MAX;
						if (r3 >= 0.5)
						{
							char buffer[20];
							sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
							string str = buffer;
							test[2 * node + 1] = str;
						}
						else
							test[2 * node + 1] = "x";
						test[node] = operator_dic[r2];
					}
					else if (r2 >= 4 && r2 <= 5 && (test[node] == "+" || test[node] == "-" || test[node] == "*" || test[node] == "/"))
					{
						Clear(test, 2 * node + 1);
						test[node] = operator_dic[r2];
					}

					else if (r2 == 6)
					{
						Clear(test, node);
						if (node > 1)
						{
							test[node] = operator_dic[r2];
						}
					}
					else if (r2 == 7)
					{
						Clear(test, node);
						if (node > 1)
						{
							char buffer[20];
							sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
							string str = buffer;
							test[node] = str;
						}
					}


					//					cout << "Mutate Operator" << endl;
				}
				else if (test[node] == "0")
					continue;
				else
				{
					double r = (double)rand() / (double)RAND_MAX;
					if (r < 0.5)
					{
						double c = atof(test[node].c_str());
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
						test[node] = str;
					}
					else
					{
						if (node < 128)
						{
							int r = rand() % 8;

							if (r < 4)
							{
								test[node] = operator_dic[r];
								//								for (int i = node * 2; i <= node * 2 + 1; i++)
								//								{
								//									double r1 = (double)rand() / (double)RAND_MAX;
								//									if (r1 >= 0.5)
								test[2 * node] = "x";

								//									else
								//									{
								char buffer[20];
								sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
								string str = buffer;
								test[2 * node + 1] = str;
								//									}
								//								}

							}
							else if (r == 4 || r == 5)
							{

								//								if (test[node / 2] != "sin" && test[node / 2] != "cos")
								//								{
								test[node] = operator_dic[r];
								double r1 = (double)rand() / (double)RAND_MAX;
								if (r1 >= 0.5)
									test[node * 2] = "x";

								else
								{
									char buffer[20];
									sprintf_s(buffer, "%.10f", (double)rand() / (double)RAND_MAX * 20.0 - 10.0);
									string str = buffer;
									test[node * 2] = str;
								}

								//								}


							}
							else if (r == 6)
							{
								test[node] = "x";
							}
							else if (r == 7)
								continue;
						}
					}
				}
			}
			else if (ran >= 0.8 && node > 1)
			{
				string tank = "\0";
				Clear(test, node);
				test[node] = "x";
			}
			break;
		}
	}
}

void Clear(string gene[], int x)
{
	gene[x] = "0";
	if (2 * x < 256)
	{
		Clear(gene, 2 * x);
		Clear(gene, 2 * x + 1);
	}
}
