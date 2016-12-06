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
    double dlon = (location2.mLongitude * 0.0174533) - (location1.mLongitude * 0.0174533);
    double dlat = (location2.mLatitude * 0.0174533) - (location1.mLatitude * 0.0174533);
    double a = std::pow((sin(dlat/2)),2) + cos(location1.mLatitude) * cos(location2.mLatitude) * std::pow((sin(dlon/2)),2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    return 3961 * c;
}

double SumLocationDistances(double fitness, Location loc)
{
    return fitness + loc.mLatitude;
}

std::vector<std::pair<int, double> > CalculatePopulationFitness(Population population, std::vector<Location> locations)
{
    std::vector<std::pair<int, double> > pairs;

    std::for_each(population.mMembers.begin(), population.mMembers.end(), [&](std::vector<int>& member){
        std::vector<double> stopsDistance;
        stopsDistance.resize(locations.size()+1);
        std::pair<int, double> pair;
        double fitness = 0;
        
        /*std::vector<Location> orderedLocations;
        std::for_each(member.begin(), member.end(), [&orderedLocations, &locations] (int n){
            orderedLocations.push_back(locations[n]);
        });*/
        
        std::adjacent_difference(member.begin(), member.end(), stopsDistance.begin(), [&locations] (int loc1, int loc2){
            return CalculateHaversineDistance(loc2, loc1, locations);
        });
        
        //std::transform(orderedLocations.begin(), orderedLocations.end(), orderedLocations.begin()+1, stopsDistance.begin(), CalculateHaversineDistance);
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
