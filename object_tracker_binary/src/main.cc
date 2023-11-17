#include <object_tracker/object_tracker.h>

#include <iostream>
#include <sstream>
#include <vector>

int main() {
  /// Read number of objects and time steps from user interface
  std::string input;
  char comma;  // to store the comma character

  size_t num_obj{0};
  size_t times{0};
  {
    do {
      std::cout << "Enter <number of objects>,<time steps>: ";
      std::getline(std::cin, input);  // Read a line of input from the user
      std::stringstream ss(input);
      ss >> num_obj >> comma >> times;

      if ((num_obj == 0) || (times == 0)) {
        fprintf(
            stderr,
            "Invalid input: number of objects and time steps must be non zero");
      }
    } while ((num_obj == 0) || (times == 0));
  }

  std::vector<tracker::ObjectTracker> trackers;

  for (size_t i = 0; i < num_obj; i++) {
    float x_top, y_top, x_bot, y_bot;

    /// Read initial from user interface
    std::cout
        << "Enter initial state of bound <x_top>,<y_top>,<x_bot>,<y_bot>: ";
    std::getline(std::cin, input);  // Read a line of input from the user
    std::stringstream ss(input);
    ss >> x_top >> comma >> y_top >> comma >> x_bot >> comma >> y_bot;

    tracker::ObjectTracker obj(x_top, y_top, x_bot, y_bot);
    trackers.emplace_back(obj);
  }

  for (size_t i = 0; i < times; i++) {
    for (auto &tracker : trackers) {
      float x_top, y_top, x_bot, y_bot;

      /// Read measurements from user interface
      std::cout << "Enter measurement <x_top>,<y_top>,<x_bot>,<y_bot>: ";
      std::getline(std::cin, input);  // Read a line of input from the user
      std::stringstream ss(input);
      ss >> x_top >> comma >> y_top >> comma >> x_bot >> comma >> y_bot;

      const mat_t x_meas = {{(x_top + x_bot) / 2.0f}};
      const mat_t y_meas = {{(y_top + y_bot) / 2.0f}};

      tracker.estimate(x_meas, y_meas);

      float xhat, yhat;
      tracker.get_position(xhat, yhat);
      fprintf(stdout, "%f,%f\n", xhat, yhat);
    }
  }
  exit(EXIT_SUCCESS);
}
