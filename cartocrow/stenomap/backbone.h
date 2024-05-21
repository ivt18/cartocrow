#ifndef CARTOCROW_STENOMAP_BACKBONE_H
#define CARTOCROW_STENOMAP_BACKBONE_H

#include "../core/core.h"
#include <CGAL/Polygon_2.h>


namespace cartocrow::backbone{

template <typename K> using AdjacencyList = std::map<Point<K>, std::list<Point<K>>>;
template <typename K> using WeightList = std::map< std::pair<Point<K>, Point<K>>, double>;

class Backbone {
  private:
    Polygon<Inexact> pgon;
    AdjacencyList<Inexact> graph;
    WeightList<Inexact> distances;
	

  public:
    Backbone(const Polygon<Inexact>& shape);
    void build_visibility_graph();
    void add_vertex(const Point<Inexact>& s);
    void add_edge(const Point<Inexact>& s, const Point<Inexact>& t, double dist);
	

};

} // namespace cartocrow::backbone

#endif //CARTOCROW_STENOMAP_BACKBONE_H
