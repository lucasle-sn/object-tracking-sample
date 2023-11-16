#include <cassert>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>

constexpr float cSampleTime = 0.1f;
constexpr float cGamma = 0.01f;

template <typename T>
void mprint(std::vector<std::vector<T>> A) {
  for (const auto &row : A) {
    for (auto &val : row) {
      std::cout << val << " ";
    }
    std::cout << "\n";
  }
}

template <typename T>
std::vector<std::vector<T>> mmultiply(const std::vector<std::vector<T>> &A,
                                      const std::vector<std::vector<T>> &B) {
  const size_t nA_row = A.size();
  const size_t nA_col = A[0].size();
  const size_t nB_row = B.size();
  const size_t nB_col = B[0].size();

  if (nA_col != nB_row) {
    fprintf(stderr,
            "Number of columns of matrix A must be equal to number of rows of "
            "matrix B\n");
    assert(false);
  }

  std::vector<std::vector<T>> ret(nA_row, std::vector<T>(nB_col, 0));
  for (size_t i = 0; i < nA_row; i++) {
    for (size_t j = 0; j < nB_col; j++) {
      for (size_t k = 0; k < nA_col; k++) {
        ret[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return ret;
}

template <typename T, typename U>
std::vector<std::vector<T>> mmultiply(const std::vector<std::vector<T>> &A,
                                      const U coff) {
  const size_t nA_row = A.size();
  const size_t nA_col = A[0].size();

  std::vector<std::vector<T>> ret(nA_row, std::vector<T>(nA_col, 0));
  for (size_t i = 0; i < nA_row; i++) {
    for (size_t j = 0; j < nA_col; j++) {
      ret[i][j] = A[i][j] * coff;
    }
  }
  return ret;
}

template <typename T, typename... Args>
std::vector<std::vector<T>> mmultiply(const std::vector<std::vector<T>> &A,
                                      const std::vector<std::vector<T>> &B,
                                      Args... args) {
  return mmultiply(A, mmultiply(B, args...));
}

template <typename T>
std::vector<std::vector<T>> madd(const std::vector<std::vector<T>> &A,
                                 const std::vector<std::vector<T>> &B) {
  const size_t nA_row = A.size();
  const size_t nA_col = A[0].size();
  const size_t nB_row = B.size();
  const size_t nB_col = B[0].size();

  if ((nA_row != nB_row) || (nA_col != nB_col)) {
    fprintf(stderr, "Two matrices must have the same size\n");
    assert(false);
  }

  std::vector<std::vector<T>> ret(nA_row, std::vector<T>(nA_col, 0));
  for (size_t i = 0; i < nA_row; i++) {
    for (size_t j = 0; j < nA_col; j++) {
      ret[i][j] = A[i][j] + B[i][j];
    }
  }
  return ret;
}

template <typename T, typename... Args>
std::vector<std::vector<T>> madd(const std::vector<std::vector<T>> &A,
                                 const std::vector<std::vector<T>> &B,
                                 Args... args) {
  return madd(A, madd(B, args...));
}

template <typename T>
std::vector<std::vector<T>> msubtract(const std::vector<std::vector<T>> &A,
                                      const std::vector<std::vector<T>> &B) {
  return madd(A, mmultiply(B, -1.0f));
}

template <typename T>
std::vector<std::vector<T>> mtranspose(const std::vector<std::vector<T>> &A) {
  const size_t nA_row = A.size();
  const size_t nA_col = A[0].size();

  std::vector<std::vector<T>> ret(nA_col, std::vector<T>(nA_row, 0));
  for (size_t i = 0; i < nA_row; i++) {
    for (size_t j = 0; j < nA_col; j++) {
      ret[j][i] = A[i][j];
    }
  }
  return ret;
}

template <typename T = float>
std::vector<std::vector<T>> midentity(size_t size) {
  std::vector<std::vector<T>> ret(size, std::vector<T>(size, 0));
  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      if (i == j) {
        ret[i][j] = static_cast<T>(1);
      }
    }
  }
  return ret;
}

typedef std::vector<std::vector<float>> mat_t;

/// Continuous model (observer c.f)
const struct {
  mat_t A = {{-cGamma, 1}, {0, 0}};
  mat_t B = {{0}, {1}};
  mat_t C = {{1, 0}};
  mat_t D = {{0}};
} cConModel;

const mat_t Psi =
    madd(midentity<float>(2), mmultiply(cConModel.A, cSampleTime / 2.0f),
         mmultiply(cConModel.A, cConModel.A, cSampleTime *cSampleTime / 6));

// Having x[k+1] = A*x[k] + B*u[k]; y[k] = C*x[k] + D*u[k]
const struct {
  mat_t A = madd(midentity<float>(2), mmultiply(cConModel.A, Psi, cSampleTime));
  mat_t B = mmultiply(Psi, cConModel.B, cSampleTime);
  mat_t C = cConModel.C;
  mat_t D = {{0}};
} cDisModel;

// Covariance matrices V (model) & W (measurement)
const struct {
  mat_t model = {{1, 0}, {0, 0.1f}};
  mat_t meas = {{1}};
} cCov;

class ObjectTracker {
  struct object_cord {
    mat_t state = {{0}, {0}};
    mat_t P = {{0, 0}, {0, 0}};
    //    mat_t meas = {{0}};
    struct {
      mat_t prior = {{0}, {0}};
      mat_t P_prior = {{0, 0}, {0, 0}};
    } internal;
  };

 public:
  ObjectTracker() noexcept = delete;
  explicit ObjectTracker(float x_top, float y_top, float x_bot, float y_bot,
                         float u_x = 0.0f, float u_y = 0.0f) noexcept {
    x.state = {{(x_top + x_bot) / 2.0f}, {u_x}};
    y.state = {{(y_top + y_bot) / 2.0f}, {u_y}};
  }

  bool estimate(const mat_t &x_measure, const mat_t &y_measure,
                float u_x = 0.0f, float u_y = 0.0f) {
    return (estimate(x.state, x_measure, x.P, u_x) &&
            estimate(y.state, y_measure, y.P, u_y));
  }

  void predict(float u_x = 0.0f, float u_y = 0.0f) {
    predict(x.state, x.P, x.internal.prior, x.internal.P_prior, u_x);
    predict(y.state, y.P, y.internal.prior, y.internal.P_prior, u_y);
  }

  bool update(const mat_t &x_measure, const mat_t &y_measure) {
    return (
        update(x.internal.prior, x.internal.P_prior, x.state, x_measure, x.P) &&
        update(y.internal.prior, y.internal.P_prior, y.state, y_measure, y.P));
  }

  void get_position(float &x_pos, float &y_pos) noexcept {
    x_pos = x.state[0][0];
    y_pos = y.state[0][0];
  }

 private:
  void predict(const mat_t &xhat, const mat_t &P, mat_t &xhat_prior,
               mat_t &P_prior, float u = 0.0f) {
    // Prediction
    xhat_prior = madd(mmultiply(cDisModel.A, xhat), mmultiply(cDisModel.B, u));
    P_prior =
        madd(mmultiply(cDisModel.A, P, mtranspose(cDisModel.A)), cCov.model);
  }

  bool update(const mat_t &xhat_prior, const mat_t &P_prior, mat_t &xhat,
              const mat_t &measure, mat_t &P) {
    // Update
    auto tmp = madd(mmultiply(cDisModel.C, P_prior, mtranspose(cDisModel.C)),
                    cCov.meas);
    if (std::abs(tmp[0][0]) < 1e-5) {
      fprintf(stderr, "C*P_prior*transpose(C) = 0 - Can't be inverted\n");
      return false;
    }

    auto K = mmultiply(P_prior, mtranspose(cDisModel.C), 1.0f / tmp[0][0]);
    xhat = madd(
        xhat_prior,
        mmultiply(K, msubtract(measure, mmultiply(cDisModel.C, xhat_prior))));
    P = mmultiply(msubtract(midentity<float>(2), mmultiply(K, cDisModel.C)),
                  P_prior);
    return true;
  }

  bool estimate(mat_t &xhat, const mat_t &measure, mat_t &P, float u = 0.0f) {
    // Prediction
    auto xhat_prior =
        madd(mmultiply(cDisModel.A, xhat), mmultiply(cDisModel.B, u));
    auto P_prior =
        madd(mmultiply(cDisModel.A, P, mtranspose(cDisModel.A)), cCov.model);

    // Update
    auto tmp = madd(mmultiply(cDisModel.C, P_prior, mtranspose(cDisModel.C)),
                    cCov.meas);
    if (std::abs(tmp[0][0]) < 1e-5) {
      fprintf(stderr, "C*P_prior*transpose(C) = 0 - Can't be inverted\n");
      return false;
    }

    auto K = mmultiply(P_prior, mtranspose(cDisModel.C), 1.0f / tmp[0][0]);
    xhat = madd(
        xhat_prior,
        mmultiply(K, msubtract(measure, mmultiply(cDisModel.C, xhat_prior))));
    P = mmultiply(msubtract(midentity<float>(2), mmultiply(K, cDisModel.C)),
                  P_prior);
    return true;
  }

  struct object_cord x;
  struct object_cord y;
};

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

  std::vector<ObjectTracker> trackers;

  for (size_t i = 0; i < num_obj; i++) {
    float x_top, y_top, x_bot, y_bot;

    /// Read initial from user interface
    std::cout
        << "Enter initial state of bound <x_top>,<y_top>,<x_bot>,<y_bot>: ";
    std::getline(std::cin, input);  // Read a line of input from the user
    std::stringstream ss(input);
    ss >> x_top >> comma >> y_top >> comma >> x_bot >> comma >> y_bot;

    ObjectTracker obj(x_top, y_top, x_bot, y_bot);
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
