#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>

using namespace std;

struct City {
    int id;
    int x, y;
};

vector<City> generateCities(int num, int maxX, int maxY) {
    vector<City> cities;
    srand(time(0));

    for (int i = 0; i < num; ++i) {
        int x = rand() % (maxX + 1);
        int y = rand() % (maxY + 1);
        cities.push_back({i + 1, x, y});
    }

    return cities;
}

void writeCitiesToFile(const vector<City> &cities, const string &filename) {
    ofstream file(filename);

    if (file.is_open()) {
        file << cities.size() << endl;
        for (const City &city : cities) {
            file << city.id << " " << city.x << " " << city.y << endl;
        }
        file.close();
        cout << "Cities written to " << filename << endl;
    } else {
        cerr << "Unable to open file " << filename << endl;
    }
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        cerr << "Usage: " << argv[0] << " <num_cities> <max_x> <max_y> <output_file>" << endl;
        return 1;
    }

    int numCities = atoi(argv[1]);
    int maxX = atoi(argv[2]);
    int maxY = atoi(argv[3]);
    string outputFile = argv[4];

    vector<City> cities = generateCities(numCities, maxX, maxY);
    writeCitiesToFile(cities, outputFile);

    return 0;
}