#ifndef CARTOCROW_STENOMAP_MEDIAL_AXIS_H
#define CARTOCROW_STENOMAP_MEDIAL_AXIS_H

#include "../core/core.h"
#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_2/IO/print.h>
#include <CGAL/draw_straight_skeleton_2.h>

/* namespace cartocrow::tree {

// TODO: put this somewhere else
class TreeNode {
    private:        
        Point<Inexact> vertex;
        std::vector<TreeNode*> children;

    public:
        // Constructs a new tree node with vertex v
        TreeNode(const Point<Inexact>& v);
        // Adds node n as a child
        void add_child(TreeNode *n);
        // Returns const iterator to beginning of children of node
        std::vector<TreeNode*>::const_iterator children_begin();
        // Returns const iterator to end of children of node
        std::vector<TreeNode*>::const_iterator children_end();
};

class Tree {
    private:
        TreeNode root;

    public:
        Tree();
        TreeNode get_root();
};

} */

// TODO: properly define all of this
namespace cartocrow::medial_axis {
typedef CGAL::Straight_skeleton_2<Inexact> Ss;
typedef boost::shared_ptr<Ss> SsPtr;

template<typename K> using AdjacencyList = std::map<Point<K>, std::list<Point<K>>>;

class MedialAxis {
  private:
    // pointer to interior straight skeleton
    SsPtr iss;
    // adjacency list
    AdjacencyList<Inexact> graph;

  public:
	// Constructs a new medial axis given single polygon
	MedialAxis(const Polygon<Inexact>& shape);
    // Prints medial axis to IO
    void print_medial_axis();
    // Adds vertex to the adjacency list
    void add_vertex(const Point<Inexact>& s);
    // Adds edge to the adjacency list
    void add_edge(const Point<Inexact>& s, const Point<Inexact>& t);
    // Prints the adjacency list
    void print_adjacency_list();
};

} // namespace cartocrow::necklace_map

#endif //CARTOCROW_STENOMAP_MEDIAL_AXIS_H
