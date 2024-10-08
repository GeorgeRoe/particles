#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// reads particle data from the particle data file
void read_file(std::vector<double>& posx, std::vector<double>& posy, std::vector<double>& posz) {
  /*
   * files are formatted as such:
   * the first line contains how many rows worth of data the file stores
   * the subsequent lines store data in columns for the x, y and z positions of the particles
   */

  // open the file
  std::ifstream ParticleData("particle_data.txt");
  std::string line;

  // read the first line to find how many particles should be read
  std::getline(ParticleData, line);
  int file_length = std::stoi(line);

  // loop over the lines of the file and parse the text as double precision positions
  int i = 0;
  while (i < file_length && std::getline(ParticleData, line)) {
    posx.push_back(std::stod(line.substr(1, 19)));
    posy.push_back(std::stod(line.substr(25, 19)));
    posz.push_back(std::stod(line.substr(49, 19)));
    i++;
  }

  // close the file
  ParticleData.close();
}


// finds the distance between two points on an axis through a PBC boundary
double PBC_distance(double a,  double b, double lower_boundary, double upper_boundary) {
  // find the distance between each point and its closest boundary wall
  double min_a = std::min(upper_boundary - a, a - lower_boundary);
  double min_b = std::min(upper_boundary - b, b - lower_boundary);

  // sum the two distances to find the PBC distance
  return min_a + min_b;
}

// finds the distance between two points ignoring PBC
double distance(double a, double b) {
  return std::abs(a - b);
}

// finds the shortest distance between two points on an axis using PBC
double shortest_distance(double a, double b, double lower_boundary, double upper_boundary) {
  // find the distances between the particles using both methods and return the shortest
  double pbc = PBC_distance(a, b, lower_boundary, upper_boundary);
  double standard = distance(a, b);

  return std::min(pbc, standard);
}

// uses three dimensional pythagoras to find the magnitude of a vector
double magnitude(double vector[3]) {
  double result = 0;

  // add the square of each axis' value
  for (int axis = 0; axis < 3; axis++) {
    result += vector[axis] * vector[axis];
  }

  // square root the result to finish pythag
  return sqrt(result);
}

// counts the number of particles considered pairs given a cutoff
int count_pairs(
  std::vector<double>& posx, std::vector<double>& posy, std::vector<double>& posz,
  double lower_boundary[3], double upper_boundary[3], double cutoff) {

  // variables used for iterating over data
  double i_pos[3], j_pos[3], difference[3];

  int pairs = 0;

  // triangular nested loop to find every combination of pairs without double counting
  for (int i = 0; i < posx.size(); i++) {
    i_pos[0] = posx[i];
    i_pos[1] = posy[i];
    i_pos[2] = posz[i];

    for (int j = i + 1; j < posx.size(); j++) {
      j_pos[0] = posx[j];
      j_pos[1] = posy[j];
      j_pos[2] = posz[j];

      // find the difference in position between the two particles on each axis
      for (int axis = 0; axis < 3; axis++) {
        difference[axis] = shortest_distance(i_pos[axis], j_pos[axis], lower_boundary[axis], upper_boundary[axis]);
      }

      // calculate the magnitude of the resulting vector, if it is less than the cutoff then the particles will be considered a pair
      if (magnitude(difference) < cutoff) {
        pairs += 1;
      }
    }
  }

  return pairs;
}

int main() {
  // the bounds of the simulation (max position and minimu position)
  double lower_boundary[3] = { 0, 0, 0 }, upper_boundary[3] = { 1, 1, 1 };

  // the maximum distance between two particles before they arent considered a pair
  double cutoff = 0.05;

  // lists containing particle positions
  std::vector<double> posx, posy, posz; 
  read_file(posx, posy, posz);

  // count the pairs and print the result
  int pairs = count_pairs(posx, posy, posz, lower_boundary, upper_boundary, cutoff);
  std::cout << pairs << "\n";

  return 0;
}
