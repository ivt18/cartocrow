#include "backbone_painting.h"

using namespace cartocrow;

BackbonePainting::BackbonePainting(medial_axis::MedialAxis medial_axis)
    : _medial_axis(medial_axis) {};

void BackbonePainting::paint(renderer::GeometryRenderer& renderer) const {
    renderer.setMode(renderer::GeometryRenderer::stroke);
    renderer.setStroke(Color{255, 0, 255}, 4);
    // Iterate over adjacency list to draw medial axis
    std::vector<Point<Inexact>> backbone = _medial_axis.get_backbone();
    for (int i = 0; i < backbone.size() - 1; i++) {
        Segment<Inexact> segment(backbone[i], backbone[i+1]);
        renderer.draw(segment);
    }
}

