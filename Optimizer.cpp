#include "Optimizer.h"

#include <cfloat>
#include <iostream>
#include <windows.h>

using namespace std;

COptimizer::COptimizer(CEvaluator &cEvaluator)
	: c_evaluator(cEvaluator)
{
	random_device c_seed_generator;
	c_rand_engine.seed(c_seed_generator());
	
	d_current_best_fitness = 0;

	vector<vector<double>*>* populations = new vector<vector<double>*>[NUMBER_OF_POPULATIONS];
	vector<double>* populationsQuality = new vector<double>[NUMBER_OF_POPULATIONS];

}

void COptimizer::vInitialize() {

	d_current_best_fitness = DBL_MAX;
	v_current_best.clear();

	diffWeight = DEFAULT_DIFF_WEIGHT;
	crossProb = DEFAULT_CROSS_PROB;
	populationNumber = DEFAULT_POPULATION_NUMBER;
	genotypeSize = c_evaluator.iGetNumberOfDimensions();
	runCounter = 0;
	migrationPause = false;

	initPopulation();
	stuckPoint = d_current_best_fitness;
	stuckCount = 0;

}

void COptimizer::vRunIteration()
{

	for (int i = 0; i < NUMBER_OF_POPULATIONS; i++)
	{
		iterateForPopulation(i);
	}

	runCounter++;

	// wszystkie optymalizacje po danym okresie
	if (runCounter % MIGRATION_INTERVAL == 0) 
	{

		// przygotowanie najlepszych rozwi¹zañ z ka¿dej populacji
		int bestIndexes[NUMBER_OF_POPULATIONS];

		for (int i = 0; i < NUMBER_OF_POPULATIONS; i++)
		{
			bestIndexes[i] = getBestSolutionOffset(i);

			cout << populationsQuality[i][bestIndexes[i]] << " ";
		}

		cout << endl;

		if (runCounter % (MIGRATION_INTERVAL * LARGE_PERIOD_INTERVAL / 2) == 0)
		{
			migrationPause = false;

			if (runCounter % (MIGRATION_INTERVAL * LARGE_PERIOD_INTERVAL) == 0) 
			{

				if (runCounter > MIGRATION_INTERVAL * LARGE_PERIOD_INTERVAL)
				{
					// zerowanie populacji która utknê³a
					for (int i = 0; i < NUMBER_OF_POPULATIONS; i++)
					{
						if (lastIterationBestQuality[i] < populationsQuality[i][bestIndexes[i]] * MINIMAL_EXPECTED_IMPROVEMENT_FRACTION &&
							i != bestPopulation)
						{
							randomizePopulation(i);
							migrationPause = true;
							cout << "RANDOMIZING: " << i << endl;
						}
					}

				}

				// dane do pó¿niejszej oceny
				for (int i = 0; i < NUMBER_OF_POPULATIONS; i++)
					lastIterationBestQuality[i] = populationsQuality[i][bestIndexes[i]];
			}
		}
		
		// zmniejszenie prawdopodobieñstwa krzy¿owania w przypadku utkniêcia
		if (stuckPoint < d_current_best_fitness * MINIMAL_EXPECTED_IMPROVEMENT_FRACTION) 
		{
			stuckCount++;

			if (stuckCount % ALT_PROB_TRIGGER_POINT == 0)
			{
				crossProb = getRandomDouble(ALT_CROSS_PROB_MIN, ALT_CROSS_PROB_MAX);

				cout << "ALT PROB " << crossProb << endl;
			}
			
		} 
		else
		{
			crossProb = DEFAULT_CROSS_PROB;
			stuckCount = 0;
			stuckPoint = d_current_best_fitness;
		}

		// wprowadzenie losowego rozwi¹zania 
		for (int i = 0; i < NUMBER_OF_POPULATIONS; i++) 
		{
			if (bestIndexes[i] != populationNumber - 1)
				v_fill_randomly(*populations[i][populationNumber - 1]);
			else
				v_fill_randomly(*populations[i][populationNumber - 2]);
		}

		// migracje 
		if (!migrationPause || runCounter % (MIGRATION_INTERVAL * LARGE_PERIOD_INTERVAL * LARGE_PERIOD_INTERVAL) == 0)
		{
			vector<double>* temp = populations[0][bestIndexes[0]];

			for (int i = 0; i < NUMBER_OF_POPULATIONS - 1; i++)
				populations[i][bestIndexes[i]] = populations[i + 1][bestIndexes[i + 1]];

			populations[NUMBER_OF_POPULATIONS - 1][bestIndexes[NUMBER_OF_POPULATIONS - 1]] = temp;
		}

		cout << "Runcounter: " << runCounter << endl;
	} 
}

void COptimizer::v_fill_randomly(vector<double> &vSolution)
{
	vSolution.resize(genotypeSize);

	double d_lower_bound, d_upper_bound;

	for (int i = 0; i < genotypeSize; i++)
	{
		c_evaluator.bGetLowerBound(i, d_lower_bound);
		c_evaluator.bGetUpperBound(i, d_upper_bound);

		uniform_real_distribution<double> c_uniform_real_distribution(d_lower_bound, d_upper_bound);

		vSolution[(size_t)i] = c_uniform_real_distribution(c_rand_engine);
	}
}

COptimizer::~COptimizer()
{
	if (populations != NULL)
	{
		delete [] populations;
	}

	if (populationsQuality != NULL)
	{
		delete [] populationsQuality;
	}

}

void COptimizer::initPopulation() 
{

	for (int i = 0; i < NUMBER_OF_POPULATIONS; i++) 
	{
		populations[i].resize(populationNumber);
		populationsQuality[i].resize(populationNumber);

		for (int j = 0; j < populationNumber; j++)
			populations[i][j] = new vector<double>[populationNumber];

		randomizePopulation(i);
	}
}

void COptimizer::randomizePopulation(int populationOffset)
{
	for (int i = 0; i < populationNumber; i++) 
	{
		v_fill_randomly(*populations[populationOffset][i]);
		populationsQuality[populationOffset][i] = c_evaluator.dEvaluate(populations[populationOffset][i]);
	}

	int bestIndex = getBestSolutionOffset(populationOffset);

	if (d_current_best_fitness > populationsQuality[populationOffset][bestIndex]) 
	{
		v_current_best = *populations[populationOffset][bestIndex];
		d_current_best_fitness = populationsQuality[populationOffset][bestIndex];
	}
}

void COptimizer::iterateForPopulation(int populationOffset) 
{

	double oldDiffWieght = diffWeight;

	diffWeight += (populationOffset - NUMBER_OF_POPULATIONS / 2) / (NUMBER_OF_POPULATIONS - 1.) * diffWeight;

	for (int individualOffset = 0; individualOffset < populationNumber; ++individualOffset) 
	{

		int baseIndividualOffset = getRandomInt(0, populationNumber - 1);
		int add0IndividualOffset = getRandomInt(0, populationNumber - 1);
		int add1IndividualOffset = getRandomInt(0, populationNumber - 1);

		if (indexesAreDifferent(individualOffset, -1, add0IndividualOffset, add1IndividualOffset)) 
		{
			vector<double>* newIndividual = new vector<double>;
			newIndividual->resize(genotypeSize);

			for (int geneOffset = 0; geneOffset < genotypeSize; ++geneOffset) 
			{
				if (getRandomDouble(0, 1) < crossProb) 
				{

					(*newIndividual)[geneOffset] =
							(*(populations[populationOffset])[baseIndividualOffset])[geneOffset] +
							diffWeight * ((*(populations[populationOffset])[add0IndividualOffset])[geneOffset] -
							(*(populations[populationOffset])[add1IndividualOffset])[geneOffset]);
				}
				else (*newIndividual)[geneOffset] = (*(populations[populationOffset])[individualOffset])[geneOffset];
			}

			double fitness = c_evaluator.dEvaluate(newIndividual);

			if (fitness < populationsQuality[populationOffset][individualOffset]) 
			{
				populations[populationOffset][individualOffset]->assign(newIndividual->begin(), newIndividual->end());
				populationsQuality[populationOffset][individualOffset] = fitness;
			} 
			
			delete newIndividual;
		}
		else individualOffset--;
	}

	int bestIndex = getBestSolutionOffset(populationOffset);

	if (populationsQuality[populationOffset][bestIndex] < d_current_best_fitness) 
	{
		d_current_best_fitness = populationsQuality[populationOffset][bestIndex];
		v_current_best = *(populations[populationOffset])[bestIndex];
		bestPopulation = populationOffset;

		lastIterationBestQuality[populationOffset]++;

		cout << "Population: " << populationOffset << " " << d_current_best_fitness << endl;
	}

	diffWeight = oldDiffWieght;
}

bool COptimizer::indexesAreDifferent(int first, int second, int third, int fourth) 
{
	return first != second && first != third && first != fourth && second != third && second != fourth && third != fourth;
}

int COptimizer::getBestSolutionOffset(int populationOffset)
{
	int bestIndex = 0;
	for (int i = 1; i < populationNumber; ++i) 
	{
		if (populationsQuality[populationOffset][bestIndex] > populationsQuality[populationOffset][i])
			bestIndex = i;
	}
	return bestIndex;
}

double COptimizer::getRandomDouble(double lowerBound, double upperBound)
{
	uniform_real_distribution<double> c_uniform_real_distribution(lowerBound, upperBound);

	return c_uniform_real_distribution(c_rand_engine);
}

int COptimizer::getRandomInt(int lowerBound, int upperBound)
{
	uniform_int_distribution<int> c_uniform_int_distribution(lowerBound, upperBound);

	return c_uniform_int_distribution(c_rand_engine);
}
