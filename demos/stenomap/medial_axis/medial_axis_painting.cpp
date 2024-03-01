#include "medial_axis_painting.h"

#include "cartocrow/core/core.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_renderer.h"

#include <CGAL/create_straight_skeleton_2.h>

using namespace cartocrow;
typedef CGAL::Straight_skeleton_2<Inexact> Ss;
typedef boost::shared_ptr<Ss> SsPtr;

MedialAxisPainting::MedialAxisPainting(SsPtr medial_axis) 
    : medialAxis(medial_axis) {}

void MedialAxisPainting::paint(renderer::GeometryRenderer& renderer) const {
    renderer.setMode(renderer::GeometryRenderer::stroke);
    renderer.setStroke(Color{100,100,100}, 4);
    std::cout<<"OK";
    // TODO: implement this when drawing medial axis is added (technically a straight skeleton)
    /* renderer.draw(medialAxis); */
    // Iterate through halfedges to draw segments
    for(auto edge = medialAxis->halfedges_begin(); edge != medialAxis->halfedges_end(); ++edge) {
        if (!edge->is_border()) {
            auto start = edge->prev()->vertex()->point();
            auto end = edge->vertex()->point();
            // Create a segment from start to end point and draw it
            Segment<Inexact> segment(start, end);
            renderer.draw(segment);
        }
    }
}

