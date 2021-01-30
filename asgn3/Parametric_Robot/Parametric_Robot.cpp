#define _USE_MATH_DEFINES
#include <iostream>
#include <ctime>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <vector>
#include <stdarg.h>
#include <fstream>
#include "omp.h"
#include <algorithm>
#include <math.h>
#include <cmath>
#include <numeric>
#define _MATH_DEFINES_DEFINED
#define GRAPHICS

#define LEN 8192 //  Maximum length of text string


using namespace std;
// Physics Simluator Variables
int mass_number = 16;
int spring_number = 56;
double mass = 0.2;
double length = 0.1;
double gravity = 9.81;
double T = 0;
double timestep = 0.001;
double k_c = 100000;
double k = 6000;
double damping_constant = 0.99;
double friction_coefficient = 0.8;
double W = 2 * M_PI;
double OriLength[16];
int evospring = 12;

int generationNumber = 1;
int robotNumber = 20;
int simulationTime = 4;

float cross_rate = 0.6;
float mutation_rate = 0.1;
double maxpath = 0.0;
vector<double> Fitness;

double asp = 1;     // aspect ratio
int fov = 45;         //  Field of view (for perspective)
double dim = 1.0;  // size of the workd
int moveX, moveY;
int spinX = 0;
int spinY = 0;
int des = 0;
GLfloat world_rotation[16] = { 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1 };
clock_t start = clock();
double duration;
double start_time = 0.0;

double shiny = 1.0;
int mode = 1;
double th = 0.0 * M_PI / 180.0;            //  Azimuth of view angle
double ph = 45.0 * M_PI / 180.0;           //  Elevation of view angle

double view = 1000;
double viewlr = 90 * M_PI / 180.0;

static GLint Frames = 0;
static GLfloat fps = -1;
static GLint T0 = 0;
ofstream bestGene;
ofstream popDis("Best.txt");

struct Mass
{
	double m;
	double p[3];
	double v[3];
	double a[3];
};

struct Spring
{
	double k;
	double l0;
	int m1;
	int m2;
};

struct GENE
{
	// L = L*L_0 + A*sin(W*t+B) + C*sin(W*t+D)
	double A;
	double B;
	double C;
	double D;
	double L;
};

struct Cube
{
	vector<Mass> Masses;
	vector<Spring> Springs;
};


double distance(Mass a, Mass b)
{
	return sqrt(pow(a.p[0] - b.p[0], 2) + pow(a.p[1] - b.p[1], 2) + pow(a.p[2] - b.p[2], 2));
}



vector<int> sort_indexes(const std::vector<double>& v)
{
	// initialize original index locations
	std::vector<int> idx(v.size());
	iota(idx.begin(), idx.end(), 0);
	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),
		[&v](int i1, int i2) {return v[i1] > v[i2]; });
	return idx;
}

class ROBOT
{
private:
	std::vector<GENE> gene;
public:
	double initialLocation[3] = { 0,0,0 };
	std::vector<Mass> Masses;
	std::vector<Spring> Springs;

	void generate_Masses(double X, double Y, double Z)
	{
		Masses.push_back({ mass,{X + 1.5 * length,Y + 1.5 * length,Z},{0,0,0},{0,0,0} });
		Masses.push_back({ mass,{X + 1.5 * length,Y + 1.5 * length,Z + length},{0,0,0},{0,0,0} });
		Masses.push_back({ mass,{X - 1.5 * length,Y + 1.5 * length,Z},{0,0,0},{0,0,0} });
		Masses.push_back({ mass,{X - 1.5 * length,Y + 1.5 * length,Z + length},{0,0,0},{0,0,0} });
		Masses.push_back({ mass,{X - 1.5 * length,Y - 1.5 * length,Z},{0,0,0},{0,0,0} });
		Masses.push_back({ mass,{X - 1.5 * length,Y - 1.5 * length,Z + length},{0,0,0},{0,0,0} });
		Masses.push_back({ mass,{X + 1.5 * length,Y - 1.5 * length,Z},{0,0,0},{0,0,0} });
		Masses.push_back({ mass,{X + 1.5 * length,Y - 1.5 * length,Z + length},{0,0,0},{0,0,0} });

		Masses.push_back({ mass,{X + 1 * length,Y + 1 * length,Z + length},{0,0,0},{0,0,0} });
		Masses.push_back({ mass,{X + 1 * length,Y + 1 * length,Z + 2.0 * length},{0,0,0},{0,0,0} });
		Masses.push_back({ mass,{X - 1 * length,Y + 1 * length,Z + length},{0,0,0},{0,0,0} });
		Masses.push_back({ mass,{X - 1 * length,Y + 1 * length,Z + 2.0 * length},{0,0,0},{0,0,0} });
		Masses.push_back({ mass,{X - 1 * length,Y - 1 * length,Z + length},{0,0,0},{0,0,0} });
		Masses.push_back({ mass,{X - 1 * length,Y - 1 * length,Z + 2.0 * length},{0,0,0},{0,0,0} });
		Masses.push_back({ mass,{X + 1 * length,Y - 1 * length,Z + length},{0,0,0},{0,0,0} });
		Masses.push_back({ mass,{X + 1 * length,Y - 1 * length,Z + 2.0 * length},{0,0,0},{0,0,0} });

		/*
		Masses.push_back({ mass,{X,Y - sqrt(3) / 2.0 * length,Z},{0,0,0},{0,0,0}});

		Masses.push_back({ mass,{X + 1.0 * length, Y + sqrt(3) / 2.0 * length,Z},{0,0,0},{0,0,0} });

		Masses.push_back({ mass,{X - 1.0 * length, Y + sqrt(3) / 2.0 * length,Z},{0,0,0},{0,0,0} });

		Masses.push_back({ mass,{X, Y - sqrt(3) / 2.0 * length,Z + 1.0 * length},{0,0,0},{0,0,0} });

		Masses.push_back({ mass,{X + 1.0 * length, Y + sqrt(3) / 2.0 * length,Z + 1.0 * length},{0,0,0},{0,0,0} });

		Masses.push_back({ mass,{X - 1.0 * length, Y + sqrt(3) / 2.0 * length,Z + 1.0 * length},{0,0,0},{0,0,0} });

		Masses.push_back({ mass,{X,Y + sqrt(3) / 6.0 * length,Z +(sqrt(6) / 3.0 + 1) * length},{0,0,0},{0,0,0}});
		Masses.push_back({ mass,{X + 0.5 * length,Y,Z + length},{0,0,0},{0,0,0}});
		Masses.push_back({ mass,{X - 0.5 * length,Y,Z + length},{0,0,0},{0,0,0}});
		Masses.push_back({ mass,{X,Y + sqrt(3) / 2 * length,Z + length},{0,0,0},{0,0,0}});
		*/
	}

	void generate_Springs()
	{
		int count = 0;
		Springs.push_back({ 2000,distance(Masses[1],Masses[9]), 1,9 });
		Springs.push_back({ 2000,distance(Masses[1],Masses[10]), 1,10 });
		Springs.push_back({ 2000,distance(Masses[1],Masses[14]), 1,14 });
		Springs.push_back({ 2000,distance(Masses[3],Masses[8]), 3,8 });
		Springs.push_back({ 2000,distance(Masses[3],Masses[11]), 3,11 });
		Springs.push_back({ 2000,distance(Masses[3],Masses[12]), 3,12 });
		Springs.push_back({ 2000,distance(Masses[5],Masses[10]), 5,10 });
		Springs.push_back({ 2000,distance(Masses[5],Masses[13]), 5,13 });
		Springs.push_back({ 2000,distance(Masses[5],Masses[14]), 5,14 });
		Springs.push_back({ 2000,distance(Masses[7],Masses[8]), 7,8 });
		Springs.push_back({ 2000,distance(Masses[7],Masses[12]), 7,12 });
		Springs.push_back({ 2000,distance(Masses[7],Masses[15]), 7,15 });

		Springs.push_back({ k,distance(Masses[0],Masses[8]), 0,8 });
		Springs.push_back({ k,distance(Masses[2],Masses[10]), 2,10 });
		Springs.push_back({ k,distance(Masses[4],Masses[12]), 4,12 });
		Springs.push_back({ k,distance(Masses[6],Masses[14]), 6,14 });

		Springs.push_back({ k,distance(Masses[0],Masses[9]), 0,9 });
		Springs.push_back({ k,distance(Masses[2],Masses[11]), 2,11 });
		Springs.push_back({ k,distance(Masses[4],Masses[13]), 4,13 });
		Springs.push_back({ k,distance(Masses[6],Masses[15]), 6,15 });

		Springs.push_back({ k,distance(Masses[0],Masses[1]), 0,1 });
		Springs.push_back({ k,distance(Masses[2],Masses[3]), 2,3 });
		Springs.push_back({ k,distance(Masses[4],Masses[5]), 4,5 });
		Springs.push_back({ k,distance(Masses[6],Masses[7]), 6,7 });
		Springs.push_back({ k,distance(Masses[1],Masses[8]), 1,8 });
		Springs.push_back({ k,distance(Masses[3],Masses[10]), 3,10 });
		Springs.push_back({ k,distance(Masses[5],Masses[12]), 5,12 });
		Springs.push_back({ k,distance(Masses[7],Masses[14]), 7,14 });

		for (int i = 8; i < Masses.size() - 1; i++)
		{
			for (int j = i + 1; j < Masses.size(); j++)
			{
				Springs.push_back({ k,distance(Masses[i],Masses[j]), i,j });
			}
		}
		for (int i = 0; i < 12; i++)
		{
			OriLength[i] = Springs[i].l0;
		}
	}

	ROBOT(double X, double Y, double Z, vector<GENE> legGene)
	{
		// default constructor
		initialLocation[0] = X; initialLocation[1] = Y; initialLocation[2] = Z;
		gene = legGene;
		generate_Masses(X, Y, Z);
		generate_Springs();
	}

	void draw_cube()
	{
		GLfloat color[6][3] = { {1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0},
							 {1.0,1.0,0.0},{1.0,0.0,1.0},{0.0,1.0,1.0} };
		GLUquadric* quad;
		quad = gluNewQuadric();
		for (int i = 0; i < (int)Masses.size(); i++)
		{
			glPushMatrix();
			glMultMatrixf(world_rotation);
			glColor3f(1, 0, 1);
			glTranslated(Masses[i].p[0], Masses[i].p[1], Masses[i].p[2]);
			gluSphere(quad, 0.005, 20, 20);
			glPopMatrix();
		}
		for (int i = 0; i < (int)Springs.size(); i++)
		{
			glPushMatrix();
			glMultMatrixf(world_rotation);
			glLineWidth(2.0f);
			glBegin(GL_LINES);
			glColor3d(1, 1, 1);
			glVertex3f(Masses[Springs[i].m1].p[0], Masses[Springs[i].m1].p[1], Masses[Springs[i].m1].p[2]);
			glVertex3f(Masses[Springs[i].m2].p[0], Masses[Springs[i].m2].p[1], Masses[Springs[i].m2].p[2]);
			glEnd();
			glPopMatrix();
		}


		glPushMatrix();
		glMultMatrixf(world_rotation);
		glBegin(GL_QUADS);
		glColor3f(0.2, 0.1, 0.3);
		glVertex3f(GLfloat(Masses[8].p[0]), GLfloat(Masses[8].p[1]), GLfloat(Masses[8].p[2]));
		glVertex3f(GLfloat(Masses[9].p[0]), GLfloat(Masses[9].p[1]), GLfloat(Masses[9].p[2]));
		glVertex3f(GLfloat(Masses[10].p[0]), GLfloat(Masses[10].p[1]), GLfloat(Masses[10].p[2]));
		glVertex3f(GLfloat(Masses[11].p[0]), GLfloat(Masses[11].p[1]), GLfloat(Masses[11].p[2]));
		glEnd();
		glBegin(GL_QUADS);
		glColor3f(0.2, 0.1, 0.3);
		glVertex3f(GLfloat(Masses[10].p[0]), GLfloat(Masses[10].p[1]), GLfloat(Masses[10].p[2]));
		glVertex3f(GLfloat(Masses[11].p[0]), GLfloat(Masses[11].p[1]), GLfloat(Masses[11].p[2]));
		glVertex3f(GLfloat(Masses[12].p[0]), GLfloat(Masses[12].p[1]), GLfloat(Masses[12].p[2]));
		glVertex3f(GLfloat(Masses[13].p[0]), GLfloat(Masses[13].p[1]), GLfloat(Masses[13].p[2]));
		glEnd();
		glBegin(GL_QUADS);
		glColor3f(0.2, 0.1, 0.3);
		glVertex3f(GLfloat(Masses[12].p[0]), GLfloat(Masses[12].p[1]), GLfloat(Masses[12].p[2]));
		glVertex3f(GLfloat(Masses[13].p[0]), GLfloat(Masses[13].p[1]), GLfloat(Masses[13].p[2]));
		glVertex3f(GLfloat(Masses[14].p[0]), GLfloat(Masses[14].p[1]), GLfloat(Masses[14].p[2]));
		glVertex3f(GLfloat(Masses[15].p[0]), GLfloat(Masses[15].p[1]), GLfloat(Masses[15].p[2]));
		glEnd();
		glBegin(GL_QUADS);
		glColor3f(0.2, 0.1, 0.3);
		glVertex3f(GLfloat(Masses[8].p[0]), GLfloat(Masses[8].p[1]), GLfloat(Masses[8].p[2]));
		glVertex3f(GLfloat(Masses[10].p[0]), GLfloat(Masses[10].p[1]), GLfloat(Masses[10].p[2]));
		glVertex3f(GLfloat(Masses[12].p[0]), GLfloat(Masses[12].p[1]), GLfloat(Masses[12].p[2]));
		glVertex3f(GLfloat(Masses[14].p[0]), GLfloat(Masses[14].p[1]), GLfloat(Masses[14].p[2]));
		glEnd();
		glBegin(GL_QUADS);
		glColor3f(0.2, 0.1, 0.3);
		glVertex3f(GLfloat(Masses[9].p[0]), GLfloat(Masses[9].p[1]), GLfloat(Masses[9].p[2]));
		glVertex3f(GLfloat(Masses[11].p[0]), GLfloat(Masses[11].p[1]), GLfloat(Masses[11].p[2]));
		glVertex3f(GLfloat(Masses[13].p[0]), GLfloat(Masses[13].p[1]), GLfloat(Masses[13].p[2]));
		glVertex3f(GLfloat(Masses[15].p[0]), GLfloat(Masses[15].p[1]), GLfloat(Masses[15].p[2]));
		glEnd();
		glPopMatrix();


		// draw line between middle point and initial position
		double x = 0; double y = 0; double z = 0;
		for (int j = 0; j < mass_number; j++) {
			x = x + Masses[j].p[0];
			y = y + Masses[j].p[1];
			z = z + Masses[j].p[2];
		}
		x = x / mass_number;
		y = y / mass_number;
		z = z / mass_number;
		glPushMatrix();
		glMultMatrixf(world_rotation);
		glBegin(GL_LINES);
		glColor3f(0.0, 1.0, 0.0);
		glVertex3f(GLfloat(x), GLfloat(y), 0.0);
		glVertex3f(GLfloat(initialLocation[0]), GLfloat(initialLocation[1]), GLfloat(0.0));
		glEnd();
		glPopMatrix();
		double distance = sqrt(pow(x, 2) + pow(y, 2));
		//printf("Time: %f, Distance: %f\n", T, distance);
	}

	void robotUpdate()
	{
		bool move = true;

		vector<vector<double>> forces(mass_number, vector<double>(3));


		for (int i = 0; i < 12; i++)
		{
			Springs[i].l0 = OriLength[i] + gene[i].A * OriLength[i] * sin(W * T + gene[i].B) +
				gene[i].C * OriLength[i] * sin(2 * W * T + gene[i].D);
			if (Springs[i].l0 < 0)
			{
				cout << Springs[i].l0 << endl;
			}
		}

		for (int i = 0; i < spring_number; i++)
		{
			Mass mass1 = Masses[Springs[i].m1];
			Mass mass2 = Masses[Springs[i].m2];
			double position_distance[3] = { mass2.p[0] - mass1.p[0],mass2.p[1] - mass1.p[1],mass2.p[2] - mass1.p[2] };
			double l_now = sqrt(pow(position_distance[0], 2) + pow(position_distance[1], 2) + pow(position_distance[2], 2));
			double this_force = Springs[i].k * fabs(Springs[i].l0 - l_now);
			int flag = 1;

			if (l_now > Springs[i].l0) {
				flag = -1;
			}
			forces[Springs[i].m1][0] = forces[Springs[i].m1][0] - flag * this_force * position_distance[0] / l_now;
			forces[Springs[i].m1][1] = forces[Springs[i].m1][1] - flag * this_force * position_distance[1] / l_now;
			forces[Springs[i].m1][2] = forces[Springs[i].m1][2] - flag * this_force * position_distance[2] / l_now;
			forces[Springs[i].m2][0] = forces[Springs[i].m2][0] + flag * this_force * position_distance[0] / l_now;
			forces[Springs[i].m2][1] = forces[Springs[i].m2][1] + flag * this_force * position_distance[1] / l_now;
			forces[Springs[i].m2][2] = forces[Springs[i].m2][2] + flag * this_force * position_distance[2] / l_now;
		}
		for (int i = 0; i < mass_number; i++)
		{
			forces[i][2] = forces[i][2] - Masses[i].m * gravity;
			if (Masses[i].p[2] <= 0)
			{
				forces[i][2] = forces[i][2] + k_c * fabs(Masses[i].p[2]);
				double F_H = sqrt(pow(forces[i][0], 2) + pow(forces[i][1], 2));
				double F_V = forces[i][2];
				if (abs(F_H) < abs(F_V) * friction_coefficient)
				{
					forces[i][0] = 0;
					forces[i][1] = 0;
					Masses[i].v[0] = 0;
					Masses[i].v[1] = 0;
				}

				else
				{
					for (int j = 0; j < 2; j++)
					{
						if (forces[i][j] < 0)
						{
							forces[i][j] = forces[i][j] + abs(F_V) * friction_coefficient * abs(forces[i][j] / F_H);
							if (forces[i][j] > 0)
								forces[i][j] = 0;
						}
						else
						{
							forces[i][j] = forces[i][j] - abs(F_V) * friction_coefficient * abs(forces[i][j] / F_H);
							if (forces[i][j] < 0)
								forces[i][j] = 0;
						}
					}
				}

			}
			for (int j = 0; j < 3; j++) {
				Masses[i].a[j] = forces[i][j] / Masses[i].m;
				Masses[i].v[j] = damping_constant * (Masses[i].v[j] + Masses[i].a[j] * timestep);
				Masses[i].p[j] = Masses[i].p[j] + Masses[i].v[j] * timestep;
			}
			//double velocity = sqrt(pow(Masses[i].v[0], 2) + pow(Masses[i].v[1], 2) + pow(Masses[i].v[2], 2));
			//cout << "v:" << velocity << endl;

		};
	}
};

class Simulation {
private:
	int populationSize;
	vector<double> populationDistance;
	vector<vector<GENE>> populationGene;
	vector<vector<GENE>> newPopulationGene;
	vector<vector<GENE>> maxgene;
	vector<ROBOT> robots;
public:
	double averageDistance;
	double maxDistance;

	Simulation(int popSize)
	{
		populationSize = popSize;
		generateGenes(populationSize);
		generateRobots();
	}

	void startSim(double time) {
		if (T < time) {
			simUpdate();
#ifdef  GRAPHICS
			simDraw();
#endif
		}
		else {
			printf("### Generation %d ###", generationNumber);
			double time = omp_get_wtime() - start_time;
			printf(" Time: %f ###\n", time);
			calculatePopulationDistance();
			selection();
			crossOver();
			mutation();
			populationGene.clear();
			populationGene.shrink_to_fit();
			populationGene = newPopulationGene;
			generationNumber++;
			robots.clear(); robots.shrink_to_fit();
			generateRobots();
			T = 0;
			start_time = omp_get_wtime();
		}
	}

	void selection() {
		vector<int> index = sort_indexes(populationDistance);
		newPopulationGene.clear();
		newPopulationGene.shrink_to_fit();
		popDis << maxpath / 8.0 << endl;
		if (populationDistance[index[0]] > maxpath)
		{
			maxpath = populationDistance[index[0]];
			maxgene.clear();
			maxgene.push_back(populationGene[index[0]]);
		}
		newPopulationGene.push_back(maxgene[0]);

		for (int i = 0; i < index.size() / 2 - 1; i++) {
			newPopulationGene.push_back(populationGene[index[i]]);
		}
		for (int i = 0; i < newPopulationGene[0].size(); i++) {
			bestGene << newPopulationGene[0][i].A << " " << newPopulationGene[0][i].B << " " << newPopulationGene[0][i].C << " " << newPopulationGene[0][i].D;
		}
		bestGene << "\n";
	}

	void crossOver() {
		for (int n = 0; n < populationGene.size() / 4; n++)
		{
			int parentIndex1 = rand() % (populationGene.size() / 2);
			int parentIndex2 = rand() % (populationGene.size() / 2);

			newPopulationGene.push_back(newPopulationGene[parentIndex1]);
			newPopulationGene.push_back(newPopulationGene[parentIndex2]);
			double cross_Probability = (double)rand() / (double)RAND_MAX;
			if (cross_Probability < cross_rate) {
				int crossOverPoint1, crossOverPoint2, crossOverPoint;
				while (true)
				{
					crossOverPoint1 = rand() % evospring;
					crossOverPoint2 = rand() % evospring;
					if (crossOverPoint1 < crossOverPoint2)
					{
						break;
					}
				}
				int l = crossOverPoint2 - crossOverPoint1;
				while (true)
				{
					crossOverPoint = rand() % evospring;
					if (crossOverPoint + l < evospring)
					{
						break;
					}
				}
				GENE tank;
				for (int i = crossOverPoint1; i <= crossOverPoint2; i++)
				{
					tank = newPopulationGene[parentIndex1][i];
					newPopulationGene[parentIndex1][i] = newPopulationGene[parentIndex2][crossOverPoint];
					newPopulationGene[parentIndex2][crossOverPoint] = tank;
					crossOverPoint += 1;
				}
			}
		}
	}

	void mutation()
	{
		for (int i = 0; i < newPopulationGene.size(); i++)
		{
			for (int j = 0; j < evospring; j++)
			{
				double mutationProbability = (double)rand() / (double)RAND_MAX;
				if (mutationProbability < mutation_rate)
				{
					newPopulationGene[i][j].A = newPopulationGene[i][j].A * (1.0 + (0.1 * (double)rand() / ((double)RAND_MAX) - 0.05));
					if (abs(newPopulationGene[i][j].A) > 0.5)
						newPopulationGene[i][j].A = 0.25;

					newPopulationGene[i][j].B = newPopulationGene[i][j].B * (1.0 + (0.1 * (double)rand() / ((double)RAND_MAX) - 0.05));

					newPopulationGene[i][j].C = newPopulationGene[i][j].C * (1.0 + (0.1 * (double)rand() / ((double)RAND_MAX) - 0.05));
					if (abs(newPopulationGene[i][j].C) > 0.5)
						newPopulationGene[i][j].C = 0.25;

					newPopulationGene[i][j].D = newPopulationGene[i][j].D * (1.0 + (0.1 * (double)rand() / ((double)RAND_MAX) - 0.05));
				}
			}
		}

	}

	void generateGenes(int number)
	{
		srand(time(0));
		for (int i = 0; i < number; i++) {
			vector<GENE> tempVec;
			for (int j = 0; j < evospring; j++) {
				double A = 0.5 * (double)rand() / ((double)RAND_MAX);;
				double B = -2 * M_PI + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (4 * M_PI)));
				double C = 0.5 * (double)rand() / ((double)RAND_MAX);
				double D = -2 * M_PI + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (4 * M_PI)));
				GENE temp{ A,B,C,D };
				tempVec.push_back(temp);
			}
			populationGene.push_back(tempVec);
		}
	}

	void generateRobots()
	{
		for (int i = 0; i < populationSize; i++) {
			double X = 2.0 * (double)rand() / (double)RAND_MAX;
			double Y = 2.0 * (double)rand() / (double)RAND_MAX;
			robots.push_back(ROBOT(X, Y, 0.0, populationGene[i]));
		}
	}

	void simUpdate()
	{
#pragma omp parallel for num_threads(16)
		for (int i = 0; i < populationSize; i++) {
			robots[i].robotUpdate();
		}
	}

	void simDraw()
	{
		for (int i = 0; i < populationSize; i++)
			robots[i].draw_cube();
	}

	void calculatePopulationDistance()
	{
		populationDistance.clear();
		populationDistance.shrink_to_fit();
		for (int i = 0; i < populationSize; i++) {
			double x = 0; double y = 0;
			for (int j = 0; j < mass_number; j++) {
				x = x + robots[i].Masses[j].p[0];
				y = y + robots[i].Masses[j].p[1];
			}
			x = x / mass_number;
			y = y / mass_number;
			//double distance = pow(fabs(x - robots[i].initialLocation[0]),2)+pow(fabs(y - robots[i].initialLocation[1]),2 );
			//populationDistance.push_back(sqrt(distance));
			double distance = fabs(x - robots[i].initialLocation[0]);
			populationDistance.push_back(distance);
		}
		averageDistance = 0;
		maxDistance = 0;
		int m = 0;
		for (int i = 0; i < populationSize; i++)
		{
			averageDistance = averageDistance + populationDistance[i];
			if (maxDistance < populationDistance[i])
			{
				maxDistance = populationDistance[i];
			}
		}

		for (int i = 0; i < evospring; i++)
		{
			cout << populationGene[m][i].A << ",";
			cout << populationGene[m][i].B << ",";
			cout << populationGene[m][i].C << ",";
			cout << populationGene[m][i].D << "," << endl;
		}
		averageDistance = averageDistance / populationSize;
		cout << "Maximum Velocity: " << maxDistance / 8.0 << endl;
		cout << "Average Velocity: " << averageDistance / 8.0 << endl;
	}
};

Simulation sim1(robotNumber);


#if 1

void Print(const char* format, ...)
{
	char buf[LEN];
	char* ch = buf;
	va_list args;
	//  Turn the parameters into a character string
	va_start(args, format);
	vsnprintf(buf, LEN, format, args);
	va_end(args);
	//  Display the characters one at a time at the current raster position
	while (*ch)
		glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, *ch++);
}

void drawGrid() {
	for (int i = -dim / 2; i < dim / 2 + 1; i++) {
		for (int j = -dim / 2; j < dim / 2 + 1; j++) {
			float white[] = { 1,1,1,1 };
			float black[] = { 0,0,0,1 };
			glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shiny);
			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, white);
			glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, black);

			glPushMatrix();
			glTranslatef(i * 2, 0, j * 2);
			glBegin(GL_QUADS);
			//
			glNormal3f(0, 1, 0);
			glColor3f(0.45, 0.45, 0.45);
			glVertex3f(+0, -0.01, +0);
			glVertex3f(+1, -0.01, +0);
			glVertex3f(+1, -0.01, +1);
			glVertex3f(+0, -0.01, +1);
			//
			glNormal3f(0, 1, 0);
			glColor3f(0.5, 0.5, 0.5);
			glVertex3f(-1, -0.01, +0);
			glVertex3f(+0, -0.01, +0);
			glVertex3f(+0, -0.01, +1);
			glVertex3f(-1, -0.01, +1);
			//
			glNormal3f(0, 1, 0);
			glColor3f(0.45, 0.45, 0.45);
			glVertex3f(-1, -0.01, -1);
			glVertex3f(+0, -0.01, -1);
			glVertex3f(+0, -0.01, +0);
			glVertex3f(-1, -0.01, +0);
			//
			glNormal3f(0, 1, 0);
			glColor3f(0.5, 0.5, 0.5);
			glVertex3f(-0, -0.01, -1);
			glVertex3f(+1, -0.01, -1);
			glVertex3f(+1, -0.01, +0);
			glVertex3f(-0, -0.01, +0);
			glEnd();
			glPopMatrix();
		}
	}

}


void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glLoadIdentity();
	const double len = 2.0;  //  Length of axes

	double Ex = -2 * dim * sin(th) * cos(ph);
	double Ey = +2 * dim * sin(ph);
	double Ez = +2 * dim * cos(th) * cos(ph);
	gluLookAt(Ex, Ey, Ez, 0, 0, 0, 0, cos(ph), 0);
	glPushMatrix();
	glRotated(spinX, 0, 1, 0);
	glRotated(spinY, 1, 0, 0);
	glTranslated(0, 0, des);

	glDisable(GL_LIGHTING);

	drawGrid();
#ifdef GRAPHICS
	glColor3f(1, 1, 1);
	//glWindowPos2i(0, 0);
	Print("Generation: %d", generationNumber);
	//glWindowPos2i(850, 0);
	Print("Max: %4.2f", sim1.maxDistance);
	//glWindowPos2i(850, 50);
	Print("Average %4.2f", sim1.averageDistance);
#endif 1
	sim1.startSim(simulationTime);
	T = T + timestep;

	Frames++;
	GLint t = glutGet(GLUT_ELAPSED_TIME);
	if (t - T0 >= 1000) {
		GLfloat seconds = ((double)t - T0) / 1000.0;
		fps = Frames / seconds;
		//printf("%d frames in %6.3f seconds = %6.3f FPS\n", Frames, seconds, fps);
		T0 = t;
		Frames = 0;
	}
	glRasterPos3d(0.0, 2, 0.0);

	glColor3f(1, 1, 1);
	glPopMatrix();
	glutSwapBuffers();
}

#endif

#if 0
void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glLoadIdentity();
	gluLookAt(-1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	glPushMatrix();
	glRotated(spinX, 0, 1, 0);
	glRotated(spinY, 1, 0, 0);
	glTranslated(0, 0, des);
	for (int i = 0; i < cube_numbers; i++) {
		Update_cube(cubes[i].Masses, cubes[i].Springs);
	}
	glPopMatrix();
	glutSwapBuffers();
}
#endif
void Project(double fov, double asp, double dim)
{

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (fov)
		gluPerspective(fov, asp, dim / 1, 16 * dim);
	else
		glOrtho(-asp * dim, asp * dim, -dim, +dim, -dim, +dim);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void mouseMove(int x, int y) {
	int dx = x - moveX;
	int dy = y - moveY;
	printf("dx;%dx,dy:%dy\n", dx, dy);
	spinX += dx;
	spinY += dy;
	glutPostRedisplay();
	moveX = x;
	moveY = y;
}
void key_pressed(unsigned char ch, int x, int y)
{
	//  Exit on ESC
	if (ch == 27)
		exit(0);
	else if (ch == '0')
		th = ph = 0;
	else if (ch == '-' && ch > 1)
		fov++;
	else if (ch == '=' && ch < 179)
		fov--;
	else if (ch == GLUT_KEY_PAGE_DOWN)
		dim += 0.1;
	else if (ch == GLUT_KEY_PAGE_UP && dim > 1)
		dim -= 0.1;
	else if (ch == 'a' || ch == 'A')
		spinX += 5;
	else if (ch == 'd' || ch == 'D')
		spinX -= 5;
	else if (ch == 'w' || ch == 'W')
		spinY += 5;
	else if (ch == 's' || ch == 'S')
		spinY -= 5;
	else if (ch == 'q' || ch == 'Q')
		viewlr -= 15;
	else if (ch == 'e' || ch == 'E')
		viewlr += 15;
	//  Reproject
	Project(fov, asp, dim);
	//  Tell GLUT it is necessary to redisplay the scene
	glutPostRedisplay();
}
void reshape(int width, int height)
{
	//  Ratio of the width to the height of the window
	asp = (double)width / height;
	glViewport(0, 0, width, height);
	Project(fov, asp, dim);
}
void idle()
{
	glutPostRedisplay();
}

int main(int argc, char* argv[])
{
#ifdef GRAPHICS
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(1000, 1200);
	glutCreateWindow("cubes");
	glutIdleFunc(idle);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMotionFunc(mouseMove);
	glutKeyboardFunc(key_pressed);
	glutMainLoop();
	return 0;
#endif

#ifndef GRAPHICS
	while (1) {
		double start_time = omp_get_wtime();
		sim1.startSim(simulationTime);
		T = T + timestep;
		//printf("%f\n", T);
		double time = omp_get_wtime() - start_time;
	}
#endif
};