#include "TSP.h"
#include <fstream>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <iostream>
#include <sstream>

std::vector<Location> GetLocations(std::string inputfile)
{
    std::ifstream::pos_type size;
    std::vector<Location> locations;
    // Open the file for input and start ATE (at the end)
    std::ifstream file (inputfile, std::ios::in|std::ios::ate);
    if (file.is_open())
    {
        size = file.tellg(); // Save the size of the file
        double length = size;
        file.seekg (0, std::ios::beg); // Seek back to start of file
        
        std::cout << "length: " << length << std::endl;
        
        std::string line;
        std::string token;
        while (std::getline(file, line)) {
            Location location;
            std::istringstream ss(line);
            std::getline(ss, token, ',');
            location.mName = token;
            std::cout << "mName: " << location.mName << std::endl;
            std::getline(ss, token, ',');
            location.mLatitude = std::stod(token);
            std::cout << "mLatitude: " << location.mLatitude << std::endl;
            std::getline(ss, token);
            location.mLongitude = std::stod(token);
            std::cout << "mLongitude: " << location.mLongitude << std::endl;
            locations.push_back(location);
        }
    }
    return locations;
}

Population GeneratePopulation(std::vector<Location> locations, int popsize, std::mt19937& randomGenerator)
{
    Population population;
    population.mMembers.resize(popsize);
    std::vector<int> initialLocations;

    initialLocations.resize(locations.size());
    int n = 0;
    std::for_each(initialLocations.begin(), initialLocations.end(), [&n](int& i){
        i = n;
        n++;
    });
    
    std::generate(population.mMembers.begin(), population.mMembers.end(), [&initialLocations, &randomGenerator] {
        std::vector<int> randomizedLocations;
        randomizedLocations = initialLocations;
        std::shuffle(randomizedLocations.begin()+1, randomizedLocations.end(), randomGenerator);
        return randomizedLocations;
    });
    
    std::ofstream outputfile;
    outputfile.open("log.txt");
    outputfile << "INITIAL POPULATION:\n";
    std::for_each(population.mMembers.begin(), population.mMembers.end(), [&outputfile](std::vector<int>& member){
        std::for_each(member.begin(), member.end()-1, [&outputfile](int& n){
            outputfile << n;
            outputfile << ",";
        });
        outputfile << *(member.end()-1);
        outputfile << "\n";
    });
    
    return population;
}

double CalculateHaversineDistance(int loc1, int loc2, const std::vector<Location>& locations)
{
    Location location1 = locations[loc1];
    Location location2 = locations[loc2];
    double latitude1 = location1.mLatitude* 0.0174533;
    double longitude1 = location1.mLongitude * 0.0174533;
    double latitude2 = location2.mLatitude* 0.0174533;
    double longitude2 = location2.mLongitude * 0.0174533;
    
    double dlon = (longitude2) - (longitude1);
    double dlat = (latitude2) - (latitude1);
    double a = std::pow((sin(dlat/2.0)),2) + cos(latitude1) * cos(latitude2) * std::pow((sin(dlon/2.0)),2);
    double c = 2.0 * atan2(sqrt(a), sqrt(1.0-a));
    return 3961.0 * c;
}

std::vector<std::pair<int, double> > CalculatePopulationFitness(Population population, std::vector<Location> locations)
{
    std::vector<std::pair<int, double> > pairs;

    std::for_each(population.mMembers.begin(), population.mMembers.end(), [&](std::vector<int>& member){
        std::vector<double> stopsDistance;
        stopsDistance.resize(locations.size());
        std::pair<int, double> pair;
        double fitness = 0;
        
        std::adjacent_difference(member.begin(), member.end(), stopsDistance.begin(), [&locations] (int loc1, int loc2){
            return CalculateHaversineDistance(loc2, loc1, locations);
        });
        
        stopsDistance[0] = CalculateHaversineDistance(*(member.end()-1), *member.begin(), locations);
        
        fitness = std::accumulate(stopsDistance.begin(), stopsDistance.end(), 0.0f);
        pair = std::make_pair(pairs.size(), fitness);
        pairs.push_back(pair);
    });
    
    std::ofstream outputfile;
    outputfile.open("log.txt", std::fstream::app);
    outputfile << "FITNESS:\n";
    std::for_each(pairs.begin(), pairs.end(), [&outputfile](std::pair<int, double>& pair){
        outputfile << pair.first;
        outputfile << ":";
        outputfile << pair.second;
        outputfile << "\n";
    });
    return pairs;
}

int GetIndex(int n, std::vector<double> probabilityVector, double runningSum, double probability)
{
    runningSum += probabilityVector[n];
    if(runningSum >= probability)
    {
        return n;
    }
    else {
        n++;
        return GetIndex(n, probabilityVector, runningSum, probability);
    }
}

std::vector<std::pair<int, int> > ParentSelection(std::vector<std::pair<int, double> > fitnessPairs, std::mt19937& randomGenerator)
{
    std::sort(fitnessPairs.begin(), fitnessPairs.end(), [](std::pair<int, double>& first, std::pair<int, double>& second){
        return first.second < second.second;
    });
    
    std::vector<double> probabilityVector;
    probabilityVector.resize(fitnessPairs.size());
    std::generate(probabilityVector.begin(), probabilityVector.end(), [&probabilityVector]{
        return 1.0/probabilityVector.size();
    });
    
    probabilityVector[fitnessPairs[0].first] *= 6.0;
    probabilityVector[fitnessPairs[1].first] *= 6.0;

    std::for_each(fitnessPairs.begin()+2, fitnessPairs.begin()+(fitnessPairs.size()/2-1), [&probabilityVector](std::pair<int, double>& pair){
        probabilityVector[pair.first] *= 3.0;
    });
    
    double sum = std::accumulate(probabilityVector.begin(), probabilityVector.end(), 0.0f);
    std::transform(probabilityVector.begin(), probabilityVector.end(), probabilityVector.begin(), [sum] (double& n){
        return n / sum;
    });
    
    std::vector<std::pair<int, int> > pairs;
    pairs.resize(fitnessPairs.size());
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    std::for_each(pairs.begin(), pairs.end(), [&probabilityVector, &distribution, &randomGenerator] (std::pair<int, int>& pair){
        double probability1 = distribution(randomGenerator);
        int index1 = GetIndex(0, probabilityVector, 0, probability1);
        
        double probability2 = distribution(randomGenerator);
        int index2 = GetIndex(0, probabilityVector, 0, probability2);
        
        pair = std::make_pair(index1, index2);
    });
    
    std::ofstream outputfile;
    outputfile.open("log.txt", std::fstream::app);
    outputfile << "SELECTED PAIRS:\n";
    std::for_each(pairs.begin(), pairs.end(), [&outputfile](std::pair<int, int>& pair){
        outputfile << "(";
        outputfile << pair.first;
        outputfile << ",";
        outputfile << pair.second;
        outputfile << ")";
        outputfile << "\n";
    });
    
    return pairs;
}

std::vector<int> CrossoverPairs(std::vector<int> parentA, std::vector<int> parentB, int size, int mutationChance, std::mt19937& randomGenerator)
{
    std::uniform_int_distribution<int> distribution1(1, size-2);
    int crossoverIndex1 = distribution1(randomGenerator);
    std::uniform_int_distribution<int> distribution2(0, 1);
    int crossoverIndex2 = distribution2(randomGenerator);
    
    std::vector<int> child;
    if(crossoverIndex2 == 1)
    {
        std::copy_n(parentA.begin(), crossoverIndex1, std::back_inserter(child));
        std::copy_if(parentB.begin(), parentB.end(), std::back_inserter(child), [&child](int elem){
            return (std::find(child.begin(), child.end(), elem) == child.end());
        });
    }
    else {
        std::copy_n(parentB.begin(), crossoverIndex1, std::back_inserter(child));
        std::copy_if(parentA.begin(), parentA.end(), std::back_inserter(child), [&child](int elem){
            return (std::find(child.begin(), child.end(), elem) == child.end());
        });
    }
    
    std::uniform_real_distribution<double> distribution3(0.0, 1.0);
    double mutationValue = distribution3(randomGenerator);
    if (mutationValue <= mutationChance) {
        std::uniform_int_distribution<int> mutationDistribution(1, size-1);
        int index1 = mutationDistribution(randomGenerator);
        int index2 = mutationDistribution(randomGenerator);
        std::swap(child[index1], child[index2]);
    }
    
    return child;
}

Population GenerateNewPopulation(std::vector<std::pair<int, int> > parents, Population population, int size, int mutationChance, std::mt19937&randomGenerator)
{
    Population newPopulation;
    population.mMembers.resize(population.mMembers.size());
    std::for_each(parents.begin(), parents.end(), [&] (std::pair<int, int> pair){
        newPopulation.mMembers.push_back(CrossoverPairs(population.mMembers[pair.first], population.mMembers[pair.second], size, mutationChance, randomGenerator));
    });
    
    
    return newPopulation;
}
