#ifndef PSTP_DEPENDENCE_LOGGERS_H
#define PSTP_DEPENDENCE_LOGGERS_H

#include <cstddef>

namespace pstp {

struct KNNLogger {
  size_t vis_node_num = 0;
  size_t generate_box_num = 0;
  size_t check_box_num = 0;
  size_t skip_box_num = 0;
};

struct RangeQueryLogger {
  size_t vis_node_num = 0;
  size_t generate_box_num = 0;
  size_t full_box_num = 0;
  size_t skip_box_num = 0;
};

}  // namespace pstp

#endif  // PSTP_DEPENDENCE_LOGGERS_H_
