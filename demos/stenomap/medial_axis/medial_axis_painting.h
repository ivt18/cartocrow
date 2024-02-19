#include "cartocrow/core/core.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_renderer.h"

#include <CGAL/create_straight_skeleton_2.h>

using namespace cartocrow;
typedef CGAL::Straight_skeleton_2<Inexact> Ss;
typedef boost::shared_ptr<Ss> SsPtr;

/// The \ref renderer::GeometryPainting "GeometryPainting" for a \ref medialAxis.
class MedialAxisPainting : public renderer::GeometryPainting {

  public:
	      /// Creates a new painting with the given polygon and medial axis
	      MedialAxisPainting(SsPtr medial_axis);

  protected:
	      void paint(renderer::GeometryRenderer& renderer) const override;

  private:
        SsPtr medialAxis;
};
