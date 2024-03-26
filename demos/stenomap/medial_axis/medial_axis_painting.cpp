#include "medial_axis_painting.h"

using namespace cartocrow;

MedialAxisPainting::MedialAxisPainting(medial_axis::MedialAxis medial_axis)
    : _medial_axis(medial_axis) {};

void MedialAxisPainting::paint(renderer::GeometryRenderer& renderer) const {
    renderer.setMode(renderer::GeometryRenderer::stroke);
    renderer.setStroke(Color{255, 0, 0}, 4);
    // // Iterate over adjacency list to draw medial axis
    // for (auto vertex : _medial_axis.get_graph()) {
    //     for (auto end : vertex.second) {
    //         // Create a segment from start to end point and draw it
    //         Segment<Inexact> segment(vertex.first, end);
    //         renderer.draw(segment);
    //     }
    // }
    // Iterate over branches to draw medial axis
    for (auto branch : _medial_axis.get_branches()) {
        for (int i = 0; i < branch.size() - 1; i++) {
            // Create a segment from start to end point and draw it
            Segment<Inexact> segment(branch[i], branch[i + 1]);
            renderer.draw(segment);
        }
    }
}

