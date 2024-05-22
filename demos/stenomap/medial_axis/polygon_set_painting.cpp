#include "polygon_set_painting.h"

#include "cartocrow/core/core.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_renderer.h"

using namespace cartocrow;

PolygonSetPainting::PolygonSetPainting(PolygonSet<Inexact> polySet) 
    : polySet(polySet) {}

void PolygonSetPainting::paint(renderer::GeometryRenderer& renderer) const {
    renderer.setMode(renderer::GeometryRenderer::stroke);
    renderer.setStroke(Color{200,200,100}, 4);
    renderer.draw(polySet);
}

