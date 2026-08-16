#ifndef PSI_DEPENDENCE_LOGGERS_H
#define PSI_DEPENDENCE_LOGGERS_H

#include <cstddef>

namespace psi
{

struct knn_logger {
	size_t vis_leaf_num = 0;
	size_t vis_interior_num = 0;
	size_t generate_box_num = 0;
	size_t check_box_num = 0;
	size_t skip_box_num = 0;
};

struct range_query_logger {
	size_t vis_leaf_num = 0;
	size_t vis_interior_num = 0;
	size_t generate_box_num = 0;
	size_t full_box_num = 0;
	size_t skip_box_num = 0;
};

} // namespace psi

#endif // PSI_DEPENDENCE_LOGGERS_H_
