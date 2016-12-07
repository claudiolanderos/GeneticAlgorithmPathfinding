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
    int popsize, generations, seed;
    double mutationchance;
    
    inputfile = argv[1];
    parsedInt = argv[2];
    popsize = std::stoi(parsedInt);
    parsedInt = argv[3];
    generations = std::stoi(parsedInt);
    parsedInt = argv[4];
    mutationchance = std::stod(parsedInt);
    mutationchance /= 100.0;
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

    for (int i = 1; i <= generations; i++) {
        population = GenerateNewPopulation(parents, population, static_cast<int>(locations.size()), mutationchance, random);
        
        std::ofstream outputfile;
        outputfile.open("log.txt", std::fstream::app);
        outputfile << "GENERATION: ";
        outputfile << i;
        outputfile << "\n";
        std::for_each(population.mMembers.begin(), population.mMembers.end(), [&outputfile](std::vector<int>& member){
            std::for_each(member.begin(), member.end()-1, [&outputfile](int& n){
                outputfile << n;
                outputfile << ",";
            });
            outputfile << *(member.end()-1);
            outputfile << "\n";
        });
        outputfile.close();
        fitnessPairs = CalculatePopulationFitness(population, locations);
        if(i+1 <= generations)
        {
            parents = ParentSelection(fitnessPairs, random);
        }
    }
    
    std::sort(fitnessPairs.begin(), fitnessPairs.end(), [](std::pair<int, double>& first, std::pair<int, double>& second){
        return first.second < second.second;
    });
    
    std::ofstream outputfile;
    outputfile << "SOLUTION: \n";
    outputfile.open("log.txt", std::fstream::app);
    std::for_each(population.mMembers[fitnessPairs[0].first].begin(), population.mMembers[fitnessPairs[0].first].end(), [&outputfile, &locations](int& n){
        outputfile << locations[n].mName;
        outputfile << "\n";
    });
    outputfile << locations[0].mName;
    outputfile << "\n";
    outputfile << "DISTANCE: ";
    outputfile << fitnessPairs[0].second;
    outputfile << " miles";
    outputfile.close();
    

	return 0;
}
