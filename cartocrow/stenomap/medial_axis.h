#ifndef CARTOCROW_STENOMAP_MEDIAL_AXIS_H
#define CARTOCROW_STENOMAP_MEDIAL_AXIS_H

#include "../core/core.h"
#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_with_holes_2.h>



// TODO: properly define all of this
namespace cartocrow::medial_axis {
typedef CGAL::Straight_skeleton_2<Inexact> Ss;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_Inexact;
typedef CGAL::Polygon_2<K> Polygon_2;
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
	int detail_level;
    PolygonSet<Inexact> region;
	std::vector<Polygon_2> circle_polygons;

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
	void compute_radius();
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

	void retract_end_branches(double retraction_percentage);
    void apply_modified_negative_offset(double constant_offset, double min_radius);
    Point<Inexact> interpolate(const Point<Inexact>& start, const Point<Inexact>& end, double ratio);

	Polygon_2 approximate_circle(const Point_Inexact& center, double radius, int n_sides);

	void store_points_on_medial_axis();
	std::vector<Point<Inexact>> points_on_medial_axis;
	void compute_negatively_offset_polygon();
	std::vector<Point<Inexact>> get_points_not_in_branch(int i);

    // Getters
    AdjacencyList<Inexact> get_graph() const {
        return graph;
    }
	// Return branches of the medial axis
	std::vector<Branch<Inexact>> get_branches() const {
		return branches;
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
