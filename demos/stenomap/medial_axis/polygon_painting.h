#include "cartocrow/core/core.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_renderer.h"

using namespace cartocrow;
using namespace cartocrow::renderer;

/// The \ref renderer::GeometryPainting "GeometryPainting" for a \ref medialAxis.
class PolygonPainting : public GeometryPainting {

  public:
	    /// Creates a new painting with the given polygon and medial axis
	    PolygonPainting(Polygon<Inexact> poly);

  protected:
	    void paint(GeometryRenderer& renderer) const override;

  private:
      Polygon<Inexact> poly;
};
