#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <limits>
#include <random>
#include <algorithm>
#include <chrono>
#include <atomic>
#include <thread>

using namespace std;
using namespace std::chrono;

struct ACOConfig {
    int numAnts = 50;
    double alpha = 1.0;
    double beta = 5.0;
    double evaporationRate = 0.5;
    double q = 100.0;
    int maxExecutionTimeMs = 3 * 60 * 1000; // 3 minutes default
};

struct City {
    int id;
    double x, y;
};

// Read ACO configuration from file
ACOConfig readConfigFromFile(const string& filename) {
    ACOConfig config;
    ifstream configFile(filename);
    
    if (!configFile.is_open()) {
        cerr << "Warning: Could not open config file. Using default parameters." << endl;
        return config;
    }

    string line;
    while (getline(configFile, line)) {
        istringstream iss(line);
        string key;
        double value;
        
        if (getline(iss, key, '=') && iss >> value) {
            if (key == "num_ants") config.numAnts = static_cast<int>(value);
            else if (key == "alpha") config.alpha = value;
            else if (key == "beta") config.beta = value;
            else if (key == "evaporation_rate") config.evaporationRate = value;
            else if (key == "q") config.q = value;
            else if (key == "max_execution_time_ms") config.maxExecutionTimeMs = static_cast<int>(value);
        }
    }

    return config;
}

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

// Atomic flag for time monitoring
atomic<bool> timeExceeded(false);

// Time monitoring thread function
void monitorTime(int maxExecutionTimeMs) {
    this_thread::sleep_for(milliseconds(maxExecutionTimeMs));
    timeExceeded = true;
}

// Main ACO TSP function
vector<int> acoTSP(const vector<City> &cities, const ACOConfig& config) {
    // Reset time exceeded flag
    timeExceeded = false;
    
    // Start time monitoring thread
    thread timeMonitor(monitorTime, config.maxExecutionTimeMs);
    
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
        vector<vector<int>> antsTours(config.numAnts, vector<int>(numCities));
        vector<double> antsTourLengths(config.numAnts, 0.0);

        // Construct tours for each ant
        for (int ant = 0; ant < config.numAnts; ++ant) {
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
                        double tau = pow(pheromones[currentCity][city], config.alpha);
                        // Add small constant to avoid division by zero
                        double eta = pow(1.0 / (distances[currentCity][city] + 1e-12), config.beta);
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
                pheromones[i][j] *= (1.0 - config.evaporationRate);
                // Keep pheromone above a tiny threshold (optional)
                if (pheromones[i][j] < 1e-12) {
                    pheromones[i][j] = 1e-12;
                }
            }
        }

        // Pheromone update (deposit)
        for (int ant = 0; ant < config.numAnts; ++ant) {
            for (int step = 0; step < numCities - 1; ++step) {
                int city1 = antsTours[ant][step];
                int city2 = antsTours[ant][step + 1];
                pheromones[city1][city2] += (config.q / antsTourLengths[ant]);
                pheromones[city2][city1] += (config.q / antsTourLengths[ant]);
            }
            // Close the loop
            int lastCity = antsTours[ant][numCities - 1];
            int firstCity = antsTours[ant][0];
            pheromones[lastCity][firstCity] += (config.q / antsTourLengths[ant]);
            pheromones[firstCity][lastCity] += (config.q / antsTourLengths[ant]);
        }

        iteration++;
    }

    // Detach the monitoring thread
    timeMonitor.detach();

    cout << "Total iterations before time limit: " << iteration << endl;
    return bestTour;
}

int main(int argc, char *argv[]) {
    bool showCities = false;
    string citiesFilename;
    string configFilename;
    
    // -c for showing every city
    vector<string> args;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-c") {
            showCities = true;
        } else {
            args.push_back(arg);
        }
    }

    // After parsing, the first non-flag argument is the cities file
    if (args.empty()) {
        cerr << "Usage: " << argv[0] << " [-c] <cities_file> [config_file]" << endl;
        return 1;
    }

    citiesFilename = args[0];
    ACOConfig config;

    // If a second non-flag argument is provided, it's the config file
    if (args.size() > 1) {
        configFilename = args[1];
        config = readConfigFromFile(configFilename);
    }

    // Read cities from file
    vector<City> cities = readCitiesFromFile(citiesFilename);
    if (cities.empty()) {
        cerr << "Error: No cities loaded. Check file format or path." << endl;
        return 1;
    }

    cout << "Loaded " << cities.size() << " cities from: " << citiesFilename << "\n\n";

    // Conditionally display cities based on the -c flag
    if (showCities) {
        for (const auto &c : cities) {
            cout << "ID: " << c.id << " at (" << c.x << ", " << c.y << ")\n";
        }
        cout << endl;
    }

    auto start = high_resolution_clock::now();
    vector<int> bestTour = acoTSP(cities, config);
    auto stop = high_resolution_clock::now();
    
    // Print best tour using the actual city IDs from the input file
    cout << "Best tour found (city IDs):" << endl;
    for (int idx : bestTour) {
        cout << cities[idx].id << " ";
    }
    cout << endl;

    double bestDistance = calculateTourDistance(cities, bestTour);
    cout << "Best tour distance: " << bestDistance << endl;

    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "ACO execution time: " << duration.count() << " ms" << endl;

    return 0;
}
