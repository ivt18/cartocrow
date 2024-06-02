#include "medial_axis.h"
#include "../core/core.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <cmath>
#include <unordered_map>
#include <numbers>

namespace cartocrow::medial_axis {

    MedialAxis::MedialAxis(const Polygon<Inexact>& shape) {
        assert(shape.is_counterclockwise_oriented());
        polygon = shape;

        // compute the medial axis
        compute_voronoi_diagram();
        filter_voronoi_diagram_to_medial_axis();
    }

    bool equals(Point<Inexact> p, Point<Inexact> t) {
        return (CGAL::squared_distance(p,t) < 0.01);
    }

    void MedialAxis::calculate_weight_function() {
        for (const auto& vertex: medial_axis_graph) {
            const Point<Inexact>& current_point = vertex.first;
            double radius = INFINITY;
            Segment<Inexact> segment;
            Segment<Inexact> segment2;
            for (auto edgeIt = polygon.edges_begin(); edgeIt != polygon.edges_end(); ++edgeIt) {
                Segment<Inexact> edge = *edgeIt;
                Inexact::FT dist = CGAL::squared_distance(current_point, edge);
                if (dist < radius) radius = dist;
            }
            radius_list.insert(std::make_pair(current_point, radius));
            
        }
    }

    void MedialAxis::compute_radius() {
        radius_list.clear();
        for (const auto& cur_point: points_on_medial_axis) {
            double radius = INFINITY;
            Segment<Inexact> segment;
            Segment<Inexact> segment2;
            for (auto edgeIt = polygon.edges_begin(); edgeIt != polygon.edges_end(); ++edgeIt) {
                Segment<Inexact> edge = *edgeIt;
                Inexact::FT dist = CGAL::squared_distance(cur_point, edge);
                if (dist < radius) radius = dist;
            }
            radius_list.insert(std::make_pair(cur_point, radius));
        }
    }

    void MedialAxis::compute_centroid_neighborhoods() {
        for (const auto& row : grid) {
            for (const auto& point : row) {
                for (const auto& vertex: medial_axis_graph) {
                    const Point<Inexact>& current_point = vertex.first;
                    centroid_neighborhoods[current_point];
                    double radius = radius_list[current_point];
                    if (CGAL::squared_distance(current_point, point) <= radius) {
                        centroid_neighborhoods[current_point].push_back(point);
                    }
                }
            }
        }
    }

    // Iterate through all grid points and find the closest centroid
    void MedialAxis::compute_centroid_closest_points() {
        for (const auto& row : grid) {
            for (const auto& point : row) {
                Point<Inexact> closest_centroid;
                double min_distance = INFINITY;
                for (const auto& vertex: medial_axis_graph) {
                    const Point<Inexact>& current_point = vertex.first;
                    centroid_closest_points[current_point];
                    double distance = CGAL::squared_distance(current_point, point);
                    if (distance < min_distance) {
                        min_distance = distance;
                        closest_centroid = current_point;
                    }
                }
                centroid_closest_points[closest_centroid].push_back(point);
            }
        }
    }
    // Remove all points from the grid that are not in the polygon
    void MedialAxis::prune_grid() {
        for (const auto& row : grid) {
            std::vector<Point<Inexact>> pruned_row;
            for (const auto& point : row) {
                if (CGAL::bounded_side_2(polygon.vertices_begin(), polygon.vertices_end(), point) != CGAL::ON_UNBOUNDED_SIDE) {
                    grid_pruned.push_back(point);
                }
            }
        }
    }

    void MedialAxis::add_vertex(AdjacencyList<Inexact>& graph, const Point<Inexact>& s) {
        graph[s];
        centroid_area_lost[s];
    }

    void MedialAxis::remove_vertex(AdjacencyList<Inexact>& graph, const Point<Inexact>& s) {
        // Remove s from the adjacency lists of its neighbors.
        for (const auto& neighbor : graph[s]) {
            auto& neighborsList = graph[neighbor];
            neighborsList.remove(s);
        }

        // Finally, remove the point from the graph.
        graph.erase(s);
    }

    void MedialAxis::add_edge(AdjacencyList<Inexact>& graph, const Point<Inexact>& s, const Point<Inexact>& t) {
        // Add both vertices to the adjacency list
        add_vertex(graph, s);
        add_vertex(graph, t);

        // Add the edge to the adjacency list in both directions
        graph[s].push_back(t);
        graph[t].push_back(s);
    }

    void MedialAxis::print_adjacency_list() {
        std::cout << "Printing adjacency list..." << std::endl;
        for (const auto& pair : medial_axis_graph) {
            std::cout << "Vertex (" << pair.first << ") adjacency list:\n";
            for (const auto& neighbor : pair.second) {
                std::cout << "\t(" << neighbor << ")\n";
            }
            std::cout << std::endl;
        }
        std::cout << "Finished printing adjacency list!" << std::endl;
    }

    void MedialAxis::compute_branches() {
        std::set<Point<Inexact>> visited;
        std::map<Point<Inexact>, Point<Inexact>> parents; // To reconstruct paths
        branches.clear();
  
        // Identify endpoints and junction points
        std::set<Point<Inexact>> endpoints, junctions;
        endpoints.clear();
        junctions.clear();
        for (const auto& entry : medial_axis_graph) {
            size_t degree = entry.second.size();
            if (degree == 2) { // Endpoint
                endpoints.insert(entry.first);
            } else if (degree > 4) { // Junction point
                junctions.insert(entry.first);
            }
        }

        // Ensure graph is not empty and has endpoints and junctions
        if (medial_axis_graph.empty() || endpoints.empty() || junctions.empty()) {
            std::cout << "Graph, endpoints, or junctions are empty." << std::endl;
            return;
        }

        // For each endpoint, find the nearest junction point using BFS
        for (const auto& endpoint : endpoints) {
            std::queue<Point<Inexact>> queue;
            visited.clear();
            parents.clear();
            bool found = false;

            queue.push(endpoint);
            visited.insert(endpoint);
            parents[endpoint] = endpoint; // Marking the start

            while (!queue.empty() && !found) {
                Point<Inexact> current = queue.front();
                queue.pop();

                if (junctions.find(current) != junctions.end()) {
                    // Junction found, reconstruct the path
                    std::vector<Point<Inexact>> branch;
                    for (Point<Inexact> at = current; at != endpoint; at = parents[at]) {
                        branch.push_back(at);
                    }
                    branch.push_back(endpoint); // Add the endpoint itself
                    std::reverse(branch.begin(), branch.end());
                    branches.push_back(branch);
                    found = true;
                    break;
                }

                for (const auto& adj : medial_axis_graph.at(current)) {
                    if (visited.find(adj) == visited.end()) {
                        queue.push(adj);
                        visited.insert(adj);
                        parents[adj] = current;
                    }
                }
            }
        }

        // Identify segments not in any branch
        std::set<std::pair<Point<Inexact>, Point<Inexact>>> allSegments;
        for (const auto& branch : branches) {
            for (size_t i = 0; i < branch.size() - 1; ++i) {
                allSegments.insert({branch[i], branch[i + 1]});
                allSegments.insert({branch[i + 1], branch[i]}); // Assuming undirected graph
            }
        }

        noBranchSegments.clear();
        for (const auto& entry : medial_axis_graph) {
            for (const auto& adj : entry.second) {
                if (allSegments.find({entry.first, adj}) == allSegments.end()) {
                    noBranchSegments.push_back(Segment<Inexact>(entry.first, adj));
                }
            }
        }
        std::cout << "branches size: " << branches.size() << std::endl;
    }

    void MedialAxis::remove_branch(int index) {
        Branch<Inexact> branch = branches[index];

        // If the branch only contains the junction point or is empty, nothing to remove.
        if (branch.size() <= 1) {
            return;
        }

        // Due to implementation of compute_branches(), the junction point will be the last point 
        // in the branch vector, which we should not remove when removing the branch
        for (int i = 0; i < branch.size() - 1; ++i) {
            remove_vertex(medial_axis_graph, branch[i]);
        }

        for (int j = 0; j < branch_closest_grid_points[index].size(); j++) {
            // since indices go in decreasing order we can simply call erase without worrying
            // about the indices shifting down
            grid_pruned.erase(grid_pruned.begin() + branch_closest_grid_points[index][j]);
        }
    }

     std::vector<Point<Inexact>> MedialAxis::get_points_not_in_branch(int i) {
        Branch<Inexact> branch = branches[i];
        std::vector<Point<Inexact>> points_not_in_branch;
        for (const auto& point: medial_axis_graph) {
            Point<Inexact> cur = point.first;
            for (int j = 0; j < branch.size() - 1; j++) {
                if (cur == branch[j]) {
                    break;
                }
                if (j == branch.size() - 2) {
                    points_not_in_branch.push_back(cur);
                }
            }
        }
        return points_not_in_branch;
    }
    

    // For each branch, count number of points of the grid that have the branch as the closest centroid and are not in the neighbourhood of any other centroid
    int MedialAxis::get_branch_weight(int i) {
        int weight = branch_closest_grid_points[i].size();
        for (int k = 0; k < branches[i].size()-1; k++) {
            weight += centroid_area_lost[branches[i][k]];
        }
        return weight;
    }

    void MedialAxis::prune_points(double t) {
        calculate_weight_function();
        compute_branches();
        compute_grid_closest_branches();
        int grid_size = grid_pruned.size();
        int removedArea = 0;
        while (true) {
            double minWeight = INFINITY;
            int minIndex = -1;

            std::cout << "number of branches = " << branches.size() << std::endl;
            // Find the branch with the minimum weight
            for (int i = 0; i < branches.size(); ++i) {
                double weight = get_branch_weight(i);
                if (weight < minWeight) {
                    minWeight = weight;
                    minIndex = i;
                }
            }
            if (minIndex == -1 || minWeight >= t * grid_size) {
                std::cout << "Pruning finished" << std::endl;
                return;
            }
            remove_branch(minIndex);
            compute_branches();
            compute_grid_closest_branches();
            centroid_area_lost[branches[minIndex].back()] += minWeight;
            removedArea += minWeight;
        }
    }

    void MedialAxis::compute_grid(unsigned int points_x, unsigned int points_y) {
        Point<Inexact> top_left, bottom_right;
        CGAL::Bbox_2 bbox = polygon.bbox();
        top_left = Point<Inexact>(bbox.xmin(), bbox.ymax());
        bottom_right = Point<Inexact>(bbox.xmax(), bbox.ymin());
        auto height = (top_left.y() - bottom_right.y()) / (points_y - 1);
        auto width = (bottom_right.x() - top_left.x()) / (points_x - 1);
        for (unsigned int y = 0; y < points_y; y++) {
            std::vector<Point<Inexact>> row(points_x);
            auto delta_y = y * height;
            for (unsigned int x = 0; x < points_x; x++) {
                row[x] = Point<Inexact>(top_left.x() + width * x, top_left.y() - delta_y);
            }
            grid.push_back(row);
        }
    }

    void MedialAxis::compute_grid_closest_branches() {
        grid_closest_branches.clear();
        branch_closest_grid_points.clear();
        // Make sure points inside radius of one of the inner points 
        // are not counted for any branch
        for (int p = 0; p < grid_pruned.size(); p++) {
            double rad_sqr_distance;
            for (int j = 0; j < branches.size(); j++) {
                rad_sqr_distance = CGAL::squared_distance(grid_pruned[p], branches[j].back());
                if (rad_sqr_distance <= radius_list[branches[j].back()]) {
                    grid_closest_branches[p] = -1;
                }
            }
        }
        for (int p = 0; p < grid_pruned.size(); p++) {
            double min_sqr_distance = INFINITY;
            double cur_sqr_distance;
            // Segments in a branch
            for (int j = 0; j < branches.size(); j++) {
                for (auto i = branches[j].begin(); i != branches[j].end() - 1; i++) {
                    Segment<Inexact> cur_segment(*i, *(i+1));
                    cur_sqr_distance = CGAL::squared_distance(grid_pruned[p], cur_segment);
                    if (cur_sqr_distance < min_sqr_distance && grid_closest_branches[p] != -1) {
                        min_sqr_distance = cur_sqr_distance;
                        grid_closest_branches[p] = j;
                    }
                }
            }
            // Segments in no branch
            for (int j = 0; j < noBranchSegments.size(); j++) {
                cur_sqr_distance = CGAL::squared_distance(grid_pruned[p], noBranchSegments[j]);
                if (cur_sqr_distance < min_sqr_distance && grid_closest_branches[p] != -1) {
                    grid_closest_branches[p] = -1;
                }
            }
        }
        // It is crucial that these go in decreasing order for removing the points from 
        // the grid in remove_branch()
        for (int p = grid_pruned.size()-1; p >= 0; p--) {
            int index = grid_closest_branches[p];
            if (index != -1) 
                branch_closest_grid_points[index].push_back(p);
        }
    }

    void MedialAxis::print_grid_closest_branches() {
        std::cout << "Printing grid closest branches" << std::endl;
            std::cout << "branch size:" << branch_closest_grid_points.size() << std::endl;
        for (int i = 0; i < branches.size(); i++) {
            std::cout << "branch index :" << i << std::endl;
            std::cout << "branch size:" << branch_closest_grid_points[i].size() << std::endl;
            for (int p : branch_closest_grid_points[i]) {
                std::cout << grid_pruned[p] << std::endl;
            }
            std::cout << std::endl;
        }
    }

void MedialAxis::retract_end_branches(double retraction_percentage) {
    for (int j = 0; j < branches.size(); j++) {
        // Function to retract branches
        auto retract = [&](double percentage) {
            double total_length = 0.0;
            for (size_t i = 0; i < branches[j].size() - 1; i++) {
                total_length += std::sqrt(CGAL::squared_distance(branches[j][i], branches[j][i + 1]));
            }

            double retract_length = total_length * percentage;

            double retracted_length = 0.0;
            size_t last_index = branches[j].size() - 1; 
            for (size_t i = 0; i < branches[j].size() - 1; i++) {
                double segment_length = std::sqrt(CGAL::squared_distance(branches[j][i], branches[j][i + 1]));

                if (retracted_length + segment_length >= retract_length) {
                    // Find the new point to insert based on the remaining length to retract
                    double remaining_length = retract_length - retracted_length;
                    double ratio = remaining_length / segment_length;
                    
                    Point<Inexact> new_point = interpolate(branches[j][i], branches[j][i + 1], ratio);
                    remove_vertex(medial_axis_graph, branches[j][i + 1]);
                    remove_vertex(medial_axis_graph, branches[j][i]);
                    add_edge(medial_axis_graph, new_point, branches[j][i + 1]);
                    branches[j][i] = new_point; // Replace the point at i-1 with the new point
                    last_index = i; 
                    retracted_length += segment_length;
                     
                    break;
                } else {
                    retracted_length += segment_length;
                    last_index = i; 
                }
            }

            // Remove the retracted part of the branch
            branches[j].erase(branches[j].begin(), branches[j].begin() + last_index);
        };

        // First retraction
        retract(retraction_percentage);

        // Second retraction based on the new length
        retract(retraction_percentage);
    }
}

    void MedialAxis::compute_voronoi_diagram() {
        vd.clear();

        // Add polygon boundary to voronoi diagram
        for (std::size_t i = 0; i < polygon.size(); i++) {
            vd.insert(Site_2<Inexact>::construct_site_2(polygon[i], polygon[(i + 1) % polygon.size()]));
        }

        // Ensure voronoi diagram is valid
        assert(vd.is_valid());
    }


//compute the last new leaf
Point<Inexact> MedialAxis::interpolate(const Point<Inexact>& start, const Point<Inexact>& end, double ratio) {
    auto x = start.x() + ratio * (end.x() - start.x());
    auto y = start.y() + ratio * (end.y() - start.y());
    return Point<Inexact>(x, y);
}


void MedialAxis::apply_modified_negative_offset(double constant_offset, double min_radius) {
	compute_radius();
	// iterate through each vertex-radius pair in the radius_list.
	for (auto& vertex_radius_pair : radius_list) {
		// reduce the radius by the constant offset.
		//std :: cout << "Applying modified negative offset..." << std::endl;
		double new_radius = vertex_radius_pair.second - (constant_offset * constant_offset);

		// ensure the new radius does not fall under the min threshold.
		// if it does, set it to the mini theshold.
		if (new_radius < min_radius) {
			new_radius = min_radius;
			//new_radius = 0;
		}
		// print old vs new radious
		//std::cout << "Old radius: " << vertex_radius_pair.second << " New radius: " << new_radius << std::endl;
		vertex_radius_pair.second = new_radius;
		//std::cout << new_radius << std::endl;
	}
	compute_negatively_offset_polygon();
}

void MedialAxis::store_points_on_medial_axis() {
	double threshold = 1e-1; // make this smaller maybe ;1e-5 is too small
	points_on_medial_axis.clear(); 

	for (const auto& branch : branches) {
    for (size_t i = 0; i < branch.size() - 1; ++i) {
      Point<Inexact> start = branch[i];
      Point<Inexact> end = branch[i + 1];
      points_on_medial_axis.push_back(start);
      points_on_medial_axis.push_back(end);
      for (double i = 1; i < 50; i++) {
        Point<Inexact> new_point = interpolate(start, end, i / 50.0);
        points_on_medial_axis.push_back(new_point);
      }
    }
  }
}

Polygon_2 MedialAxis::approximate_circle(const Point_Inexact& center, double radius, int n_sides) {
	Polygon_2 poly;
	for (int i = 0; i < n_sides; ++i) {
		double angle = 2 * std::numbers::pi * i / n_sides;
		poly.push_back(Point_Inexact(center.x() + radius * std::cos(angle), center.y() + radius * std::sin(angle)));
	}
	return poly;
}

void MedialAxis::compute_negatively_offset_polygon() {
	// Generate approximated circle polygons for each point with the modified radius
	for (const auto& vertex_radius_pair : radius_list) {
		Polygon_2 poly = approximate_circle(vertex_radius_pair.first, std::sqrt(vertex_radius_pair.second), 20);
		circle_polygons.push_back(poly);
	}

	for (const auto& poly : circle_polygons) {
		region.join(poly);
	}
}

    bool MedialAxis::vertex_in_polygon(const VertexHandle<Inexact>& vh) {
        // Make sure we do not check vertices more than once
        std::set<VertexHandle<Inexact>> considered;

        // Skip vertices which have already been considered, since a vertex may be connected to multiple halfedges
        if(considered.count(vh) != 0) {
            return false;
        }
        // Ensure we don't look at a vertex twice
        considered.insert(vh);

        // Determine if the vertex is inside the polygon
        const auto orientation = CGAL::oriented_side(vh->point(), polygon);
        bool inside_of_polygon = orientation == CGAL::ON_ORIENTED_BOUNDARY || orientation == CGAL::POSITIVE;

        // If the vertex was inside the polygon make a note of it
        if(inside_of_polygon) {
            return true;
        }

        return false;
    }

    std::set<VertexHandle<Inexact>> MedialAxis::identify_vertices_inside_polygon() {
        // Ensure voronoi diagram is valid
        assert(vd.is_valid());
        // Keep track of vertices inside polygon
        std::set<VertexHandle<Inexact>> inside;

        for (auto it = vd.bounded_halfedges_begin(); it != vd.bounded_halfedges_end(); it++) {
            if (vertex_in_polygon(it->source()) && !inside.contains(it->source()))
                inside.insert(it->source());
            if (vertex_in_polygon(it->target()) && !inside.contains(it->target()))
                inside.insert(it->target());
        }

        return inside;
    }

    double MedialAxis::z_cross_product(Point<Inexact> p, Point<Inexact> q, Point<Inexact> r) {
        const auto dx1 = q.x() - p.x();
        const auto dy1 = q.y() - p.y();
        const auto dx2 = r.x() - q.x();
        const auto dy2 = r.y() - q.y();
        return dx1 * dy2 - dy1 * dx2;
    }

    std::set<Point<Inexact>> MedialAxis::identify_concave_vertices_polygon() {
        // Ensure voronoi diagram is valid
        assert(vd.is_valid());
        // Keep track of concave vertices of polygon
        std::set<Point<Inexact>> concave;

        for (size_t i = 0; i < polygon.size(); i++) {
            const auto cross_product = z_cross_product(
                    polygon[(i + 0) % polygon.size()],
                    polygon[(i + 1) % polygon.size()],
                    polygon[(i + 2) % polygon.size()]
            );
            
            if (cross_product < 0)
                concave.insert(polygon[(i + 1) % polygon.size()]);
        }

        return concave;
    }

    void MedialAxis::filter_voronoi_diagram_to_medial_axis() {
        auto inside = identify_vertices_inside_polygon();
        auto concave = identify_concave_vertices_polygon();

        for (auto it = vd.bounded_halfedges_begin(); it != vd.bounded_halfedges_end(); it++) {
            const VertexHandle<Inexact> p = it->source();
            const VertexHandle<Inexact> q = it->target();

            // filter identity edges
            if (p->point() == q->point())
                continue;

            // filter voronoi diagram to only those vertices inside the polygon
            if (!inside.contains(p) || !inside.contains(q))
                continue;

            // drop those edges which are not part of the medial axis
            if (concave.contains(p->point()) || concave.contains(q->point()))
                continue;

            // add the edge and points to the data
            add_edge(medial_axis_graph, p->point(), q->point());
        }
    }
	
void MedialAxis::compute_visibility_graph(std::list<Point<Inexact>> feature_points) {
    // extract the polygons from the polygon set
    /* std::list<PolygonWithHoles<Inexact>> polygons;
    region.polygons_with_holes(std::back_inserter(polygons)); */
    
    for(auto p : feature_points) {
      for (auto q : feature_points) {
        if (p == q) continue;
        bool visible = true;
        Segment<Inexact> seg(p,q);
        /* for (const PolygonWithHoles<Inexact>& pwh : polygons) {
            const Polygon<Inexact> outer_boundary = pwh.outer_boundary();

            // iterate over edges in outer boundary
            for (auto edgeIt = outer_boundary.edges_begin(); edgeIt != outer_boundary.edges_end(); ++edgeIt) {
              Segment<Inexact> edge = *edgeIt;
              if (CGAL::do_intersect(seg, edge)) {
                visible = false;
                break;
              }
            }

            // iterate over the holes
            for (auto holeIt = pwh.holes_begin(); holeIt != pwh.holes_end(); ++holeIt) {
                for (auto edgeIt = holeIt->edges_begin(); edgeIt != holeIt->edges_end(); ++edgeIt) {
                  Segment<Inexact> edge = *edgeIt;
                  if (CGAL::do_intersect(seg, edge)) {
                    visible = false;
                    break;
                  }
                }
            }
        }
        if (visible) {
          add_edge(visibility_graph, p, q);
        } */
      }
    }
  }

} // namespace cartocrow::medial_axis
