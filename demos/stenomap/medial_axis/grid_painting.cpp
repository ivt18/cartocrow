#include "grid_painting.h"

using namespace cartocrow;

GridPainting::GridPainting(medial_axis::MedialAxis medial_axis)
    : _medial_axis(medial_axis) {};

void GridPainting::paint(renderer::GeometryRenderer& renderer) const {
    renderer.setMode(renderer::GeometryRenderer::vertices);
    renderer.setStroke(Color{0,0,255}, 1);
    for (auto row : _medial_axis.get_grid()) {
        for (auto point : row) {
            renderer.draw(point);
        }
    }
}

