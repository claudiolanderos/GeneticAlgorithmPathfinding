#pragma once
#include <string>
#include <vector>
#include <random>

struct Location
{
	std::string mName;
	double mLatitude;
	double mLongitude;
};

struct Population
{
	std::vector<std::vector<int>> mMembers;
};

std::vector<Location> GetLocations(std::string inputfile);

Population GeneratePopulation(std::vector<Location> locations, int popsize, std::mt19937& randomGenerator);

std::vector<std::pair<int, double> > CalculatePopulationFitness(Population population, std::vector<Location> locations);
