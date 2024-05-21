#include "medial_axis_painting.h"

using namespace cartocrow;

MedialAxisPainting::MedialAxisPainting(medial_axis::MedialAxisData medial_axis_data)
    : medial_axis_data(medial_axis_data) {};

void MedialAxisPainting::paint(renderer::GeometryRenderer& renderer) const {
    renderer.setMode(renderer::GeometryRenderer::stroke);
    renderer.setStroke(Color{255, 0, 0}, 4);

    for (auto edge : medial_axis_data.edges) {
        Segment<Inexact> s(edge.first, edge.second);
        renderer.draw(s);
    }
}

