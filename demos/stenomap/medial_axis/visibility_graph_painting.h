#ifndef CARTOCROW_STENOMAP_VISIBILITY_GRAPH_PAINTING_H
#define CARTOCROW_STENOMAP_VISIBILITY_GRAPH_PAINTING_H

#include "cartocrow/stenomap/medial_axis.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_renderer.h"

using namespace cartocrow;

// The \ref renderer::GeometryPainting "GeometryPainting" for a \ref medialAxis.
class VisibilityGraphPainting : public renderer::GeometryPainting {

    public:
        // Creates a new painting with the given medial axis
        VisibilityGraphPainting(medial_axis::MedialAxis medial_axis);

    protected:
        void paint(renderer::GeometryRenderer& renderer) const override;

    private:
        medial_axis::MedialAxis _medial_axis;
};

#endif // CARTOCROW_STENOMAP_MEDIAL_AXIS_PAINTING_H
