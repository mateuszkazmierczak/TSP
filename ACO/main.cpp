#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <limits>
#include <random>
#include <algorithm>
#include <chrono>
#include <atomic>
#include <thread>

using namespace std;
using namespace std::chrono;

struct City {
    int id;
    double x, y;
};

// Calculate Euclidean distance between two cities
double calculateDistance(const City &a, const City &b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

// Calculates total distance of a given tour
double calculateTourDistance(const vector<City> &cities, const vector<int> &tour) {
    double totalDistance = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        totalDistance += calculateDistance(cities[tour[i]], cities[tour[i + 1]]);
    }
    // Close the loop: last city back to first city
    totalDistance += calculateDistance(cities[tour.back()], cities[tour[0]]);
    return totalDistance;
}

// Read cities from file
vector<City> readCitiesFromFile(const string &filename) {
    vector<City> cities;
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Unable to open file " << filename << endl;
        return cities;
    }

    int numCities;
    if (file >> numCities) {
        for (int i = 0; i < numCities; ++i) {
            City city;
            file >> city.id >> city.x >> city.y;
            cities.push_back(city);
        }
    } else {
        // Fallback: read cities until EOF
        City city;
        while (file >> city.id >> city.x >> city.y) {
            cities.push_back(city);
        }
    }
    
    file.close();
    return cities;
}

// ACO parameters
const int NUM_ANTS = 50;
const double ALPHA = 0.5;  // Pheromone importance
const double BETA = 5.0;   // Distance (heuristic) importance
const double EVAPORATION_RATE = 0.5;
const double Q = 100.0;    // Pheromone deposit factor
const int MAX_EXECUTION_TIME_MS = 0.5 * 60 * 1000; // 3 minutes

// Atomic flag for time monitoring
atomic<bool> timeExceeded(false);

// Time monitoring thread function
void monitorTime() {
    this_thread::sleep_for(milliseconds(MAX_EXECUTION_TIME_MS));
    timeExceeded = true;
}

// Main ACO TSP function
vector<int> acoTSP(const vector<City> &cities) {
    // Reset time exceeded flag
    timeExceeded = false;
    
    // Start time monitoring thread
    thread timeMonitor(monitorTime);
    
    int numCities = static_cast<int>(cities.size());
    
    // Initialize pheromone matrix
    vector<vector<double>> pheromones(numCities, vector<double>(numCities, 1.0));
    
    // Precompute distances
    vector<vector<double>> distances(numCities, vector<double>(numCities, 0.0));
    for (int i = 0; i < numCities; ++i) {
        for (int j = 0; j < numCities; ++j) {
            distances[i][j] = calculateDistance(cities[i], cities[j]);
        }
    }

    // Random engine
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> rand01(0.0, 1.0);
    uniform_int_distribution<int> cityDist(0, numCities - 1);

    vector<int> bestTour;
    double bestTourLength = numeric_limits<double>::infinity();

    int iteration = 0;
    while (!timeExceeded) {
        // Each ant's tour and tour length
        vector<vector<int>> antsTours(NUM_ANTS, vector<int>(numCities));
        vector<double> antsTourLengths(NUM_ANTS, 0.0);

        // Construct tours for each ant
        for (int ant = 0; ant < NUM_ANTS; ++ant) {
            // Pick a random start city
            int currentCity = cityDist(gen);
            vector<bool> visited(numCities, false);

            antsTours[ant][0] = currentCity;
            visited[currentCity] = true;

            // Build the rest of the tour
            for (int step = 1; step < numCities; ++step) {
                vector<double> probabilities(numCities, 0.0);
                double sumProbabilities = 0.0;

                // Calculate selection probabilities
                for (int city = 0; city < numCities; ++city) {
                    if (!visited[city]) {
                        double tau = pow(pheromones[currentCity][city], ALPHA);
                        // Add small constant to avoid division by zero
                        double eta = pow(1.0 / (distances[currentCity][city] + 1e-12), BETA);
                        probabilities[city] = tau * eta;
                        sumProbabilities += probabilities[city];
                    }
                }

                int nextCity = -1;
                if (sumProbabilities > 0.0) {
                    double threshold = rand01(gen) * sumProbabilities;
                    double runningSum = 0.0;
                    for (int city = 0; city < numCities; ++city) {
                        if (!visited[city]) {
                            runningSum += probabilities[city];
                            if (runningSum >= threshold) {
                                nextCity = city;
                                break;
                            }
                        }
                    }
                }

                // If no city was chosen, pick a random unvisited city
                if (nextCity == -1) {
                    vector<int> unvisitedCities;
                    for (int city = 0; city < numCities; ++city) {
                        if (!visited[city]) {
                            unvisitedCities.push_back(city);
                        }
                    }
                    if (!unvisitedCities.empty()) {
                        nextCity = unvisitedCities[cityDist(gen) % unvisitedCities.size()];
                    }
                }

                // Move ant
                antsTours[ant][step] = nextCity;
                visited[nextCity] = true;
                antsTourLengths[ant] += distances[currentCity][nextCity];
                currentCity = nextCity;
            }

            // Close the loop for the tour
            antsTourLengths[ant] += distances[currentCity][antsTours[ant][0]];

            // Update best tour if improved
            if (antsTourLengths[ant] < bestTourLength) {
                bestTourLength = antsTourLengths[ant];
                bestTour = antsTours[ant];
            }

            // Break if time exceeded
            if (timeExceeded) break;
        }

        // Break if time exceeded
        if (timeExceeded) break;

        // Pheromone evaporation
        for (int i = 0; i < numCities; ++i) {
            for (int j = 0; j < numCities; ++j) {
                pheromones[i][j] *= (1.0 - EVAPORATION_RATE);
                // Keep pheromone above a tiny threshold (optional)
                if (pheromones[i][j] < 1e-12) {
                    pheromones[i][j] = 1e-12;
                }
            }
        }

        // Pheromone update (deposit)
        for (int ant = 0; ant < NUM_ANTS; ++ant) {
            for (int step = 0; step < numCities - 1; ++step) {
                int city1 = antsTours[ant][step];
                int city2 = antsTours[ant][step + 1];
                pheromones[city1][city2] += (Q / antsTourLengths[ant]);
                pheromones[city2][city1] += (Q / antsTourLengths[ant]);
            }
            // Close the loop
            int lastCity = antsTours[ant][numCities - 1];
            int firstCity = antsTours[ant][0];
            pheromones[lastCity][firstCity] += (Q / antsTourLengths[ant]);
            pheromones[firstCity][lastCity] += (Q / antsTourLengths[ant]);
        }

        iteration++;
    }

    // Detach the monitoring thread
    timeMonitor.detach();

    cout << "Total iterations before time limit: " << iteration << endl;
    return bestTour;
}

int main(int argc, char *argv[]) {
    // Expect filename as command-line argument
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <filename>" << endl;
        return 1;
    }

    string filename = argv[1];
    // Read cities from file
    vector<City> cities = readCitiesFromFile(filename);
    if (cities.empty()) {
        cerr << "Error: No cities loaded. Check file format or path." << endl;
        return 1;
    }

    cout << "Loaded " << cities.size() << " cities from: " << filename << "\n\n";
    for (const auto &c : cities) {
        cout << "ID: " << c.id << " at (" << c.x << ", " << c.y << ")\n";
    }
    cout << endl;

    auto start = high_resolution_clock::now();
    // Run ACO
    vector<int> bestTour = acoTSP(cities);
    auto stop = high_resolution_clock::now();

    // Print best tour
    cout << "Best tour found:" << endl;
    for (int idx : bestTour) {
        cout << idx << " ";
    }
    cout << endl;

    double bestDistance = calculateTourDistance(cities, bestTour);
    cout << "Best tour distance: " << bestDistance << endl;

    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "ACO execution time: " << duration.count() << " ms" << endl;

    return 0;
}