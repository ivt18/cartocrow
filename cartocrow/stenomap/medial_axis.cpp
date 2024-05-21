#include "medial_axis.h"
#include "../core/core.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <cmath>

namespace cartocrow::medial_axis {

    void MedialAxis::add_vertex(const Point<Inexact>& s) {
        graph[s];
        centroid_area_lost[s];
    }

    void MedialAxis::calculate_weight_function() {
        for (const auto& vertex: graph) {
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

    void MedialAxis::compute_centroid_neighborhoods() {
        for (const auto& row : grid) {
            for (const auto& point : row) {
                for (const auto& vertex: graph) {
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
                for (const auto& vertex: graph) {
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

    void MedialAxis::remove_vertex(const Point<Inexact>& s) {
        // Remove s from the adjacency lists of its neighbors.
        for (const auto& neighbor : graph[s]) {
            auto& neighborsList = graph[neighbor];
            neighborsList.remove(s);
        }

        // Finally, remove the point from the graph.
        graph.erase(s);
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
        std::cout << "Printing adjacency list..." << std::endl;
        for (const auto& pair : graph) {
            std::cout << "Vertex (" << pair.first << ") adjacency list:\n";
            for (const auto& neighbor : pair.second) {
                std::cout << "\t(" << neighbor << ")\n";
            }
            std::cout << std::endl;
        }
        std::cout << "Finished printing adjacency list!" << std::endl;
    }

    void MedialAxis::compute_branches() {
        std::cout << "computing branches" << std::endl;
        std::set<Point<Inexact>> visited;
        std::map<Point<Inexact>, Point<Inexact>> parents; // To reconstruct paths
        branches.clear();
  
        // Identify endpoints and junction points
        std::set<Point<Inexact>> endpoints, junctions;
        for (const auto& entry : graph) {
            size_t degree = entry.second.size();
            if (degree == 2) { // Endpoint
                endpoints.insert(entry.first);
            } else if (degree > 3) { // Junction point
                junctions.insert(entry.first);
            }
        }

        // Ensure graph is not empty and has endpoints and junctions
        if (graph.empty() || endpoints.empty() || junctions.empty()) {
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

                for (const auto& adj : graph.at(current)) {
                    if (visited.find(adj) == visited.end()) {
                        queue.push(adj);
                        visited.insert(adj);
                        parents[adj] = current;
                    }
                }
            }
        }
        std::cout << "done computing branches" << std::endl;
    }

    void MedialAxis::remove_branch(int index) {
        Branch<Inexact> branch = branches[index];
        // Due to implementation of compute_branches(), the junction point will be the last point 
        // in the branch vector, which we should not remove when removing the branch
        if (branch.size() <= 1) {
            // If the branch only contains the junction point or is empty, nothing to remove.
            return;
        }

        for (size_t i = 0; i < branch.size() - 1; ++i) {
            remove_vertex(branch[i]);
        }
        for (auto p : branch_closest_grid_points[branches[index]]) {
            for (auto it = grid_pruned.begin(); it != grid_pruned.end(); it++) {
                double dist = CGAL::squared_distance(*p, *it);
                if (dist < 0.01) {
                    grid_pruned.erase(it);
                    break;
                }
            }
        }
    }

     std::vector<Point<Inexact>> MedialAxis::get_points_not_in_branch(int i) {
        Branch<Inexact> branch = branches[i];
        std::vector<Point<Inexact>> points_not_in_branch;
        for (const auto& point: graph) {
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
        int weight = branch_closest_grid_points[branches[i]].size();
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
        while (true) {
            double minWeight = INFINITY;
            int minIndex = -1;

            std::cout << "number of branches = " << branches.size() << std::endl;
            // Find the branch with the minimum weight
            for (size_t i = 0; i < branches.size(); ++i) {
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
            std::cout << "removing branch " << minIndex << std::endl;
            remove_branch(minIndex);
            std::cout << "graph size = " << graph.size() << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            compute_branches();
            compute_grid_closest_branches();
            centroid_area_lost[branches[minIndex].back()] += minWeight;
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
        // TODO: optimize this to not go over all branches
        for (auto point: grid_pruned) {
            double min_sqr_distance = INFINITY;
            double cur_sqr_distance;
            double rad_sqr_distance;
            for (auto branch : branches) {
                for (auto i = branch.begin(); i != branch.end() && (i+1) != branch.end(); i++) {
                    Segment<Inexact> cur_segment(*i, *(i+1));
                    cur_sqr_distance = CGAL::squared_distance(point, cur_segment);
                    rad_sqr_distance = CGAL::squared_distance(point, branch.back());
                    if (cur_sqr_distance < min_sqr_distance && rad_sqr_distance >= radius_list[branch.back()]) {
                        min_sqr_distance = cur_sqr_distance;
                        grid_closest_branches[point] = branch;
                    }
                }
            }
        }
        for (auto& point : grid_pruned) {
            branch_closest_grid_points[grid_closest_branches[point]].push_back(&point);
        }
    }

    void MedialAxis::print_grid_closest_branches() {
        std::cout << "Printing grid closest branches" << std::endl;
            std::cout << "branch size:" << branch_closest_grid_points.size() << std::endl;
        for (int i = 0; i < branches.size(); i++) {
            std::cout << "branch index :" << i << std::endl;
            std::cout << "branch size:" << branch_closest_grid_points[branches[i]].size() << std::endl;
            for (auto point : branch_closest_grid_points[branches[i]]) {
                std::cout << point->x() << "," << point->y() << std::endl;
            }
            std::cout << std::endl;
        }
    }

    MedialAxis::MedialAxis(const Polygon<Inexact>& shape) {
        assert(shape.is_counter_clockwise_oriented());
        polygon = shape;

        // compute the medial axis
        compute_voronoi_diagram();
        filter_voronoi_diagram_to_medial_axis();
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
            add_edge(p->point(), q->point());
        }
    }

} // namespace cartocrow::medial_axis
