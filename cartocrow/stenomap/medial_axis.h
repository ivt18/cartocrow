#ifndef CARTOCROW_STENOMAP_MEDIAL_AXIS_H
#define CARTOCROW_STENOMAP_MEDIAL_AXIS_H

#include "../core/core.h"
#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_2/IO/print.h>
#include <CGAL/draw_straight_skeleton_2.h>

// TODO: properly define all of this
namespace cartocrow::medial_axis {
typedef CGAL::Straight_skeleton_2<Inexact> Ss;
typedef boost::shared_ptr<Ss> SsPtr;

class MedialAxis {
  private:
    // pointer to interior straight skeleton
    SsPtr iss = NULL;

  public:
	// Constructs a new medial axis given single polygon
	MedialAxis(const Polygon<Inexact>& shape);
};

} // namespace cartocrow::necklace_map

#endif //CARTOCROW_STENOMAP_MEDIAL_AXIS_H
