#include <iostream>
#include <random>
#include "TSP.h"
#include <fstream>
#include <algorithm>
#include <sstream>

int main(int argc, const char* argv[])
{
    std::string path = argv[1];
    std::string inputfile, parsedInt;
    int popsize, generations, mutationchance, seed;
    
    inputfile = argv[1];
    parsedInt = argv[2];
    popsize = std::stoi(parsedInt);
    parsedInt = argv[3];
    generations = std::stoi(parsedInt);
    parsedInt = argv[4];
    mutationchance = std::stoi(parsedInt);
    parsedInt = argv[5];
    seed = std::stoi(parsedInt);
    
    std::cout << "Input file: " << inputfile <<std::endl;
    std::cout << "Pop size: " << popsize <<std::endl;
    std::cout << "Generations: " << generations <<std::endl;
    std::cout << "Mutation chance: " << mutationchance <<std::endl;
    std::cout << "Seed: " << seed <<std::endl;

    std::mt19937 random(seed);
    
    std::vector<Location> locations = GetLocations(inputfile);
    Population population = GeneratePopulation(locations, popsize, random);
    std::vector<std::pair<int, double> > fitnessPairs = CalculatePopulationFitness(population, locations);
    
    std::vector<std::pair<int, int> > parents = ParentSelection(fitnessPairs, random);

	return 0;
}
