#include "pruned_medial_axis_painting.h"

using namespace cartocrow;

PrunedMedialAxisPainting::PrunedMedialAxisPainting(medial_axis::MedialAxis medial_axis)
    : _medial_axis(medial_axis) {};

void PrunedMedialAxisPainting::paint(renderer::GeometryRenderer& renderer) const {
    renderer.setMode(renderer::GeometryRenderer::stroke);
    renderer.setStroke(Color{0, 255, 0}, 4);
    // Iterate over adjacency list to draw medial axis
    for (auto vertex : _medial_axis.get_graph()) {
        for (auto end : vertex.second) {
            // Create a segment from start to end point and draw it
            Segment<Inexact> segment(vertex.first, end);
            renderer.draw(segment);
        }
    }
}

