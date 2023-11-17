#ifndef OBJECT_TRACKER_OBJECT_TRACKER_H
#define OBJECT_TRACKER_OBJECT_TRACKER_H

#include <vector>

using mat_t = std::vector<std::vector<float>>;

namespace tracker {

class ObjectTracker {
  struct object_cord {
    mat_t state = {{0}, {0}};
    mat_t P = {{0, 0}, {0, 0}};
    struct {
      mat_t prior = {{0}, {0}};
      mat_t P_prior = {{0, 0}, {0, 0}};
    } internal;
  };

 public:
  ObjectTracker() noexcept = delete;
  explicit ObjectTracker(float x_top, float y_top, float x_bot, float y_bot,
                         float u_x = 0.0f, float u_y = 0.0f) noexcept;

  bool estimate(const mat_t &x_measure, const mat_t &y_measure,
                float u_x = 0.0f, float u_y = 0.0f);

  void predict(float u_x = 0.0f, float u_y = 0.0f);

  bool update(const mat_t &x_measure, const mat_t &y_measure);

  void get_position(float &x_pos, float &y_pos) noexcept;

 private:
  void predict(const mat_t &xhat, const mat_t &P, mat_t &xhat_prior,
               mat_t &P_prior, float u = 0.0f);

  bool update(const mat_t &xhat_prior, const mat_t &P_prior, mat_t &xhat,
              const mat_t &measure, mat_t &P);

  bool estimate(mat_t &xhat, const mat_t &measure, mat_t &P, float u);

  struct object_cord x;
  struct object_cord y;
};

}  // namespace tracker

#endif  // OBJECT_TRACKER_OBJECT_TRACKER_H
