#include "medial_axis.h"

/* namespace cartocrow::tree {

    TreeNode::TreeNode(const Point<Inexact>& v) {
        vertex = v;
    }

    void TreeNode::add_child(TreeNode *n) {
        children.push_back(n);
    }

    std::vector<TreeNode*>::const_iterator TreeNode::children_begin() {
        return children.begin();
    }

    std::vector<TreeNode*>::const_iterator TreeNode::children_end() {
        return children.end();
    }

    TreeNode Tree::get_root() {
        return root;
    }
} // namespace cartocrow::tree */

namespace cartocrow::medial_axis {

    void MedialAxis::add_vertex(const Point<Inexact>& s) {
        graph[s];
    }

    void MedialAxis::add_edge(const Point<Inexact>& s, const Point<Inexact>& t) {
        // Add both vertices to the adjacency list
        add_vertex(s);
        add_vertex(t);

        // Add the edge to the adjacency list in both directions
        graph[s].push_back(t);
        graph[t].push_back(s);
    }

    void MedialAxis::print_adjacency_list() {
        for (const auto& pair : graph) {
            std::cout << "Vertex " << pair.first << " adjacency list:\n";
            for (const auto& neighbor : pair.second) {
                std::cout << "\t" << neighbor << "\n";
            }
            std::cout << std::endl;
        }
    }

    MedialAxis::MedialAxis(const Polygon<Inexact>& shape) {
        assert(shape.is_counter_clockwise_oriented());
        
        // create interior straight skeleton
        // TODO: maybe make a separate function to just compute this skeleton,
        // and call it in the constructor
        // TODO: maybe even put this under the Polygon class?
        iss = CGAL::create_interior_straight_skeleton_2(shape.vertices_begin(), shape.vertices_end());
        for (auto halfedge = iss->halfedges_begin(); halfedge != iss->halfedges_end(); ++halfedge) {
            if (halfedge->is_inner_bisector()) {
                Point<Inexact> t = halfedge->vertex()->point();
                Point<Inexact> s = halfedge->opposite()->vertex()->point();
                add_edge(s, t);
            }
        }
    }

    void MedialAxis::PrintMedialAxis() {
        assert(iss != NULL);
        CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*iss);
    }

} // namespace cartocrow::medial_axis
