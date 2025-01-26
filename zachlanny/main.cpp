#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>

using namespace std;
using namespace std::chrono;

struct City {
    int id;
    int x, y;
    bool visited = false;
};

double distance(const City &a, const City &b) {
    return sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2));
}

int findNearestCity(const vector<City> &cities, int currentCity) {
    double minDist = 1e9;
    int nearestCity = -1;
    for (int i = 0; i < cities.size(); ++i) {
        if (!cities[i].visited && i != currentCity) {
            double dist = distance(cities[currentCity], cities[i]);
            if (dist < minDist) {
                minDist = dist;
                nearestCity = i;
            }
        }
    }
    return nearestCity;
}

vector<int> nearestNeighborTSP(vector<City> &cities) {
    vector<int> tour;
    int currentCity = 0;
    cities[currentCity].visited = true;
    tour.push_back(currentCity);

    for (int i = 1; i < cities.size(); ++i) {
        int nearestCity = findNearestCity(cities, currentCity);
        if (nearestCity != -1) {
            cities[nearestCity].visited = true;
            tour.push_back(nearestCity);
            currentCity = nearestCity;
        }
    }

    return tour;
}

double calculateTourDistance(const vector<City> &cities, const vector<int> &tour) {
    double totalDistance = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        totalDistance += distance(cities[tour[i]], cities[tour[i + 1]]);
    }
    totalDistance += distance(cities[tour.back()], cities[tour[0]]);
    return totalDistance;
}

vector<City> readCitiesFromFile(const string &filename) {
    vector<City> cities;
    ifstream file(filename);
    if (file.is_open()) {
        int numCities;
        file >> numCities;
        for (int i = 0; i < numCities; ++i) {
            City city;
            file >> city.id >> city.x >> city.y;
            cities.push_back(city);
        }
        file.close();
    } else {
        cerr << "Unable to open file " << filename << endl;
    }
    return cities;
}

int main() {
    string filename = "cities100.txt";
    vector<City> cities = readCitiesFromFile(filename);

    cout << "Read cities from file:" << endl;
    for (size_t i = 0; i < cities.size(); ++i) {
        cout << cities[i].id << " (" << cities[i].x << ", " << cities[i].y << ")" << endl;
    }

    auto start = high_resolution_clock::now();
    vector<int> tour = nearestNeighborTSP(cities);
    auto stop = high_resolution_clock::now();

    cout << "Tour: ";
    for (size_t i = 0; i < tour.size(); ++i) {
        cout << cities[tour[i]].id << " ";
    }
    cout << endl;

    double totalDistance = calculateTourDistance(cities, tour);
    cout << "Total tour distance: " << totalDistance << endl;

    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by TSP: " << duration.count() << " microseconds" << endl;

    return 0;
}