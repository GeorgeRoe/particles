#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

void read_file(std::vector<double>& posx, std::vector<double>& posy, std::vector<double>& posz) {
  std::ifstream ParticleData("particle_data.txt");
  std::string line;

  std::getline(ParticleData, line);

  int file_length = std::stoi(line);

  int i = 0;

  while (i < file_length && std::getline(ParticleData, line)) {
    posx.push_back(std::stod(line.substr(1, 19)));
    posy.push_back(std::stod(line.substr(25, 19)));
    posz.push_back(std::stod(line.substr(49, 19)));
    i++;
  }

  ParticleData.close();
}

double PBC_distance(double a,  double b, double lower_boundary, double upper_boundary) {
  double min_a = std::min(upper_boundary - a, a - lower_boundary);
  double min_b = std::min(upper_boundary - b, b - lower_boundary);

  return min_a + min_b;
}

double distance(double a, double b) {
  return std::abs(a - b);
}

double shortest_distance(double a, double b, double lower_boundary, double upper_boundary) {
  double pbc = PBC_distance(a, b, lower_boundary, upper_boundary);
  double standard = distance(a, b);

  return std::min(pbc, standard);
}

double magnitude(double vector[3]) {
  double result = 0;

  for (int axis = 0; axis < 3; axis++) {
    result += vector[axis] * vector[axis];
  }

  return sqrt(result);
}

int count_pairs(
  std::vector<double>& posx, std::vector<double>& posy, std::vector<double>& posz,
  double lower_boundary[3], double upper_boundary[3], double cutoff) {

  double i_pos[3], j_pos[3], difference[3];

  int pairs = 0;

  for (int i = 0; i < posx.size(); i++) {
    i_pos[0] = posx[i];
    i_pos[1] = posy[i];
    i_pos[2] = posz[i];

    for (int j = i + 1; j < posx.size(); j++) {
      j_pos[0] = posx[j];
      j_pos[1] = posy[j];
      j_pos[2] = posz[j];

      for (int axis = 0; axis < 3; axis++) {
        difference[axis] = shortest_distance(i_pos[axis], j_pos[axis], lower_boundary[axis], upper_boundary[axis]);
      }

      if (magnitude(difference) < cutoff) {
        pairs += 1;
      }
    }
  }

  return pairs;
}

int main() {
  double lower_boundary[3] = { 0, 0, 0 }, upper_boundary[3] = { 1, 1, 1 };
  double cutoff = 0.05;

  std::vector<double> posx, posy, posz; 

  read_file(posx, posy, posz);

  int pairs = count_pairs(posx, posy, posz, lower_boundary, upper_boundary, cutoff);

  std::cout << pairs << "\n";

  return 0;
}
