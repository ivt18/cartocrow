#ifndef CARTOCROW_STENOMAP_MEDIAL_AXIS_H
#define CARTOCROW_STENOMAP_MEDIAL_AXIS_H

#include "../core/core.h"
#include <CGAL/create_straight_skeleton_2.h>

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

template <typename K> using AdjacencyList = std::map<Point<K>, std::list<Point<K>>>;
template <typename K> using Grid = std::vector<std::vector<Point<K>>>;
template <typename K> using Branch = std::vector<Point<K>>;

template <typename K> using RadiusList = std::map<Point<K>, double>;

class MedialAxis {
  private:
	// pointer to interior straight skeleton
	SsPtr iss;
	// adjacency list

    AdjacencyList<Inexact> graph;
    
    // grid of points spanning the whole polygon
    Grid<Inexact> grid;

	// branch list
	std::vector<Branch<Inexact>> branches;

  public:
	Polygon<Inexact> polygon;
	RadiusList<Inexact> radius_list;
	// Constructs a new medial axis given single polygon
	MedialAxis(const Polygon<Inexact>& shape);
	// Adds vertex to the adjacency list
	void add_vertex(const Point<Inexact>& s);
	// Removes vertex from the adjacency list
	void remove_vertex(const Point<Inexact>& s);
	// Adds edge to the adjacency list
	void add_edge(const Point<Inexact>& s, const Point<Inexact>& t);
	// Prints the adjacency list
	void print_adjacency_list();
	// Calculates the weight function for each point
	void calculate_weight_function();
	// Compute branches
	void compute_branches();
	// Get weight of branch
	double get_branch_weight(int index);
	// Remove branch from graph
	void remove_branch(int index);
	// Get copy of graph with branch passed as parameter removed
	AdjacencyList<Inexact> temporary_remove_branch(int index);
	// Prunes the points of the medial axis, using parameter t as threshold
	void prune_points(double t);
	SsPtr getIss() const {
		return iss;
	}
    AdjacencyList<Inexact> get_graph() const {
        return graph;
    }
    Grid<Inexact> get_grid() const {
        return grid;
    }
    
    // Computes the grid of points containing the polygon
    void compute_grid(unsigned int cells_x, unsigned int cells_y, unsigned int points_per_cell_edge);
};

} // namespace cartocrow::medial_axis

#endif //CARTOCROW_STENOMAP_MEDIAL_AXIS_H
