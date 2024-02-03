#include "medial_axis.h"

namespace cartocrow::medial_axis {

    MedialAxis::MedialAxis(const Polygon<Inexact>& shape) {
        assert(shape.is_counter_clockwise_oriented());
        
        // create interior straight skeleton
        // TODO: maybe make a separate function to just compute this skeleton,
        // and call it in the constructor
        // TODO: maybe even put this under the Polygon class?
        iss = CGAL::create_interior_straight_skeleton_2(shape.vertices_begin(), shape.vertices_end());

        // print it/draw it
        // TODO: find a way of actually returning this as an object
        // TODO: also maybe make this a separate function?
        CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*iss);
        CGAL::draw(*iss);
    }

} // namespace cartocrow::medial_axis
