#ifndef CARTOCROW_STENOMAP_MEDIAL_AXIS_PAINTING_H
#define CARTOCROW_STENOMAP_MEDIAL_AXIS_PAINTING_H

#include "cartocrow/stenomap/medial_axis.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_renderer.h"

using namespace cartocrow;

// The \ref renderer::GeometryPainting "GeometryPainting" for a \ref medialAxis.
class MedialAxisPainting : public renderer::GeometryPainting {

    public:
        // Creates a new painting with the given medial axis
        MedialAxisPainting(medial_axis::MedialAxisData medial_axis_data);

    protected:
        void paint(renderer::GeometryRenderer& renderer) const override;

    private:
        medial_axis::MedialAxisData medial_axis_data;
};

#endif // CARTOCROW_STENOMAP_MEDIAL_AXIS_PAINTING_H
