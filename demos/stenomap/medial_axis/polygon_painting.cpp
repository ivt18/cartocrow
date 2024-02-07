#include "polygon_painting.h"

#include "cartocrow/core/core.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_renderer.h"

using namespace cartocrow;

PolygonPainting::PolygonPainting(Polygon<Inexact> polygon) 
    : poly(polygon) {}

void PolygonPainting::paint(renderer::GeometryRenderer& renderer) const {
    renderer.setMode(renderer::GeometryRenderer::stroke);
    renderer.setStroke(Color{100,100,100}, 4);
    renderer.draw(poly);
}

