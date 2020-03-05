#pragma once

#include "Evaluator.h"

#include <random>
#include <vector>

#define DEFAULT_CROSS_PROB 0.9
#define ALT_CROSS_PROB_MIN 0.05
#define ALT_CROSS_PROB_MAX 0.3
#define DEFAULT_POPULATION_NUMBER 50
#define DEFAULT_DIFF_WEIGHT 0.5
#define NUMBER_OF_POPULATIONS 3
#define MIGRATION_INTERVAL 1024
#define LARGE_PERIOD_INTERVAL 10
#define ALT_PROB_TRIGGER_POINT 5
#define CALAMITY_POINT 100
#define MINIMAL_EXPECTED_IMPROVEMENT_FRACTION 1.001

using namespace std;

class COptimizer
{
public:
	COptimizer(CEvaluator &cEvaluator);
	~COptimizer();

	void vInitialize();
	void vRunIteration();

	vector<double> *pvGetCurrentBest() { return &v_current_best; }

private:
	void v_fill_randomly(vector<double> &vSolution);

	CEvaluator &c_evaluator;

	double d_current_best_fitness;
	vector<double> v_current_best;

	mt19937 c_rand_engine;
	
	double crossProb;
	double diffWeight;
	int populationNumber;
	int genotypeSize;

	int runCounter;
	int bestPopulation;
	int migrationPause;

	double stuckPoint;
	int stuckCount;
	
	vector<vector<double>*>* populations;
	vector<double>* populationsQuality;
	double lastIterationBestQuality[NUMBER_OF_POPULATIONS];

	void initPopulation();
	void randomizePopulation(int populationOffset);
	void iterateForPopulation(int populationOffset);
	bool indexesAreDifferent(int first, int second, int third, int fourth);
	int getBestSolutionOffset(int populationOffset);

	double getRandomDouble(double lowerBound, double upperBound);
	int getRandomInt(int lowerBound, int upperBound);

};