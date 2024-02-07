#ifndef CARTOCROW_STENOMAP_MEDIAL_AXIS_H
#define CARTOCROW_STENOMAP_MEDIAL_AXIS_H

#include "../core/core.h"

// TODO: properly define all of this
namespace cartocrow::medial_axis {

class MedialAxis {
  public:
	// Constructs a new medial axis given single polygon
	MedialAxis(const Polygon<Inexact>& shape);
};

} // namespace cartocrow::necklace_map

#endif //CARTOCROW_STENOMAP_MEDIAL_AXIS_H
