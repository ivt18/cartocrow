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
template <typename K> using SimpleGrid = std::vector<Point<K>>;
template <typename K> using Branch = std::vector<Point<K>>;
template <typename K> using BranchClosestGridPoints = std::map<int, std::vector<int>>;
template <typename K> using GridClosestBranches = std::map<int, int>;
template <typename K> using RemovedGridPoints = std::unordered_set<Point<K>*>;

template <typename K> using RadiusList = std::map<Point<K>, double>;
template <typename K> using AreaLostList = std::map<Point<K>, int>;

class MedialAxis {
  private:
	// pointer to interior straight skeleton
	SsPtr iss;
	// adjacency list

    AdjacencyList<Inexact> graph;
	
	// map between centroid and points in its radius
	AdjacencyList<Inexact> centroid_neighborhoods;

	// map between centroid and the points that consider the centroid as the closest
	AdjacencyList<Inexact> centroid_closest_points;
    
    // grid of points spanning the whole polygon
    Grid<Inexact> grid;

	// grid of points inside the polygon with pruning
	SimpleGrid<Inexact> grid_pruned;

    // closest branches for each grid point
    GridClosestBranches<Inexact> grid_closest_branches;

    BranchClosestGridPoints<Inexact> branch_closest_grid_points;

    RemovedGridPoints<Inexact> removed_grid_points;

	// branch list
	std::vector<Branch<Inexact>> branches;
    std::vector<Segment<Inexact>> noBranchSegments; // To store segments not in any branch

	// Map between a centroid and the area lost
	AreaLostList<Inexact> centroid_area_lost;
	

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
	int get_branch_weight(int index);
	// Remove branch from graph
	void remove_branch(int index);
	// Get copy of graph with branch passed as parameter removed
	AdjacencyList<Inexact> temporary_remove_branch(int index);
	// Prunes the points of the medial axis, using parameter t as threshold
	void prune_points(double t);
    // Computes the (points_x * points_y) grid of points containing the polygon
    void compute_grid(unsigned int points_x, unsigned int points_y);
	void compute_centroid_neighborhoods();
	void compute_centroid_closest_points();

    void compute_grid_closest_branches();
    void print_grid_closest_branches();

	void prune_grid();

	std::vector<Point<Inexact>> get_points_not_in_branch(int i);

    // Getters
    AdjacencyList<Inexact> get_graph() const {
        return graph;
    }
    Grid<Inexact> get_grid() const {
        return grid;
    }

	SimpleGrid<Inexact> get_pruned_grid() const {
        return grid_pruned;
    }
};

} // namespace cartocrow::medial_axis

#endif //CARTOCROW_STENOMAP_MEDIAL_AXIS_H
