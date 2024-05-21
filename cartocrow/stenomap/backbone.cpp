#include "backbone.h"
#include "../core/core.h"
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_naive_point_location.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel               Kernel;
typedef Kernel::Point_2                                                 Point_2;
typedef Kernel::Segment_2                                               Segment_2;
typedef CGAL::Arr_segment_traits_2<cartocrow::Inexact>                              Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                                   Arrangement_2;
typedef Arrangement_2::Face_handle                                      Face_handle;
typedef Arrangement_2::Edge_const_iterator                              Edge_const_iterator;
typedef Arrangement_2::Ccb_halfedge_circulator                          Ccb_halfedge_circulator;



namespace cartocrow::backbone {

void Backbone::add_vertex(const Point<Inexact>& s) {
    graph[s];
}

void Backbone::add_edge(const Point<Inexact>& s, const Point<Inexact>& t, double dist) {
    // Add both vertices to the adjacency list
    add_vertex(s);
    add_vertex(t);

    // Add the edge to the adjacency list in both directions
    graph[s].push_back(t);
    graph[t].push_back(s);
    distances[std::make_pair(t,s)] = dist;
}


void Backbone::build_visibility_graph() {
    for (const auto& src : pgon) {
        for (const auto& tgt : pgon) {
            double dist = CGAL::squared_distance(src, tgt);
            add_edge(src, tgt, dist);
        }
    }
}

Backbone::Backbone(const Polygon<Inexact>& shape) {
    assert(shape.is_counterclockwise_oriented());
    pgon = shape;
    /* build_visibility_graph(); */

    /* for (const auto& pair : distances) { */
    /*     std::cout << "(" << pair.first.first  << "," << pair.first.second << ")" << ": " << pair.second << std::endl; */
    /* } */
      // create a polygon and put some points in it

    Arrangement_2 env;
    CGAL::insert_non_intersecting_curves(env,pgon.edges().begin(),pgon.edges().end());

    // find the face of the query point
    // (usually you may know that by other means)
    Point<Inexact> q(20,0);
    Arrangement_2::Face_const_handle * face;
    CGAL::Arr_naive_point_location<Arrangement_2> pl(env);
    CGAL::Arr_point_location_result<Arrangement_2>::Type obj = pl.locate(q);

    // The query point locates in the interior of a face
    face = boost::get<Arrangement_2::Face_const_handle> (&obj);

    // compute regularized visibility area
    // Define visibility object type that computes regularized visibility area
    typedef CGAL::Simple_polygon_visibility_2<Arrangement_2, CGAL::Tag_true> RSPV;
    Arrangement_2 regular_output;
    RSPV regular_visibility(env);
    regular_visibility.compute_visibility(q, *face, regular_output);
    std::cout << "Regularized visibility region of " << q << " has "
        << regular_output.number_of_edges()
        << " edges:" << std::endl;
    for (Edge_const_iterator eit = regular_output.edges_begin(); eit != regular_output.edges_end(); ++eit)
        std::cout << "[" << eit->source()->point() << " -> " << eit->target()->point() << "]" << std::endl;
}

} // namespace cartocrow::medial_axis
