#include <object_tracker/internal/maths.h>
#include <object_tracker/object_tracker.h>

namespace tracker {

constexpr float cSampleTime = 0.1f;
constexpr float cGamma = 0.01f;

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

ObjectTracker::ObjectTracker(float x_top, float y_top, float x_bot, float y_bot,
                             float u_x, float u_y) noexcept {
  x.state = {{(x_top + x_bot) / 2.0f}, {u_x}};
  y.state = {{(y_top + y_bot) / 2.0f}, {u_y}};
}

bool ObjectTracker::estimate(const mat_t &x_measure, const mat_t &y_measure,
                             float u_x, float u_y) {
  return (estimate(x.state, x_measure, x.P, u_x) &&
          estimate(y.state, y_measure, y.P, u_y));
}

void ObjectTracker::predict(float u_x, float u_y) {
  predict(x.state, x.P, x.internal.prior, x.internal.P_prior, u_x);
  predict(y.state, y.P, y.internal.prior, y.internal.P_prior, u_y);
}

bool ObjectTracker::update(const mat_t &x_measure, const mat_t &y_measure) {
  return (
      update(x.internal.prior, x.internal.P_prior, x.state, x_measure, x.P) &&
      update(y.internal.prior, y.internal.P_prior, y.state, y_measure, y.P));
}

void ObjectTracker::get_position(float &x_pos, float &y_pos) noexcept {
  x_pos = x.state[0][0];
  y_pos = y.state[0][0];
}

void ObjectTracker::predict(const mat_t &xhat, const mat_t &P,
                            mat_t &xhat_prior, mat_t &P_prior, float u) {
  // Prediction
  xhat_prior = madd(mmultiply(cDisModel.A, xhat), mmultiply(cDisModel.B, u));
  P_prior =
      madd(mmultiply(cDisModel.A, P, mtranspose(cDisModel.A)), cCov.model);
}

bool ObjectTracker::update(const mat_t &xhat_prior, const mat_t &P_prior,
                           mat_t &xhat, const mat_t &measure, mat_t &P) {
  // Update
  auto tmp =
      madd(mmultiply(cDisModel.C, P_prior, mtranspose(cDisModel.C)), cCov.meas);
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

bool ObjectTracker::estimate(mat_t &xhat, const mat_t &measure, mat_t &P,
                             float u) {
  // Prediction
  auto xhat_prior =
      madd(mmultiply(cDisModel.A, xhat), mmultiply(cDisModel.B, u));
  auto P_prior =
      madd(mmultiply(cDisModel.A, P, mtranspose(cDisModel.A)), cCov.model);

  // Update
  auto tmp =
      madd(mmultiply(cDisModel.C, P_prior, mtranspose(cDisModel.C)), cCov.meas);
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

}  // namespace tracker
