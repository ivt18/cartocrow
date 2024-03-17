#include "pruned_grid_painting.h"

using namespace cartocrow;

PrunedGridPainting::PrunedGridPainting(medial_axis::MedialAxis medial_axis)
    : _medial_axis(medial_axis) {};

void PrunedGridPainting::paint(renderer::GeometryRenderer& renderer) const {
    renderer.setMode(renderer::GeometryRenderer::vertices);
    renderer.setStroke(Color{255,255,0}, 1);
    for (auto point : _medial_axis.get_pruned_grid()) {
        renderer.draw(point);
    }
}

