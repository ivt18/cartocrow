#include "medial_axis.h"
#include "../core/core.h"
#include <CGAL/CORE/ExprRep.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <cmath>
#include <unordered_map>

namespace cartocrow::medial_axis {

    void MedialAxis::add_vertex(const Point<Inexact>& s) {
        graph[s];
        centroid_area_lost[s];
    }

    void MedialAxis::calculate_weight_function() {
        radius_list.clear();
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

        // Identify segments not in any branch
        std::set<std::pair<Point<Inexact>, Point<Inexact>>> allSegments;
        for (const auto& branch : branches) {
            for (size_t i = 0; i < branch.size() - 1; ++i) {
                allSegments.insert({branch[i], branch[i + 1]});
                allSegments.insert({branch[i + 1], branch[i]}); // Assuming undirected graph
            }
        }

        for (const auto& entry : graph) {
            for (const auto& adj : entry.second) {
                if (allSegments.find({entry.first, adj}) == allSegments.end()) {
                    noBranchSegments.push_back(Segment<Inexact>(entry.first, adj));
                }
            }
        }
    }

    void MedialAxis::remove_branch(int index) {
        Branch<Inexact> branch = branches[index];
        // Due to implementation of compute_branches(), the junction point will be the last point 
        // in the branch vector, which we should not remove when removing the branch
        if (branch.size() <= 1) {
            // If the branch only contains the junction point or is empty, nothing to remove.
            return;
        }

        for (int i = 0; i < branch.size() - 1; ++i) {
            remove_vertex(branch[i]);
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
        while (true) {
            double minWeight = INFINITY;
            int minIndex = -1;

            /* std::cout << "number of branches = " << branches.size() << std::endl; */
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
            std::cout << "removing branch " << minIndex << " with endpoint " << branches[minIndex][0] << std::endl;
            remove_branch(minIndex);
            /* std::cout << "graph size = " << graph.size() << std::endl; */
            std::cout << std::endl;
            std::cout << std::endl;
            /* std::cout << std::endl; */
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
        // Make sure points inside radius of one of the inner points 
        // are not counted for any branch
        for (int p = 0; p < grid_pruned.size(); p++) {
            double rad_sqr_distance;
            for (int j = 0; j < branches.size(); j++) {
                rad_sqr_distance = CGAL::squared_distance(grid_pruned[p], branches[j].back());
                if (rad_sqr_distance <= radius_list[branches[j].back()])
                    grid_closest_branches[p] = -1;
            }
        }
        Point<Inexact> testp(0,-12);
        for (int p = 0; p < grid_pruned.size(); p++) {
            double min_sqr_distance = INFINITY;
            double cur_sqr_distance;
            // Segments in a branch
            for (int j = 0; j < branches.size(); j++) {
                for (auto i = branches[j].begin(); i != branches[j].end(); i++) {
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
                    if (CGAL::squared_distance(grid_pruned[p], testp) < 0.01) {
                        std::cout << "going to: " << noBranchSegments[j] << std::endl;
                        std::cout << "from branch: " << branches[grid_closest_branches[p]][0] << std::endl;
                    }
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
        double total_length = 0.0;
        for (size_t i = 0; i < branches[j].size() - 1; i++) {
            total_length += std::sqrt(CGAL::squared_distance(branches[j][i], branches[j][i + 1]));
        }

        double retract_length = total_length * (retraction_percentage);
        
        double retracted_length = 0.0;
        size_t last_index = branches[j].size() - 1; 
        for (size_t i = 0; i < branches[j].size() - 1; i++) {
            double segment_length = std::sqrt(CGAL::squared_distance(branches[j][i], branches[j][i + 1]));

            if (retracted_length + segment_length >= retract_length) {
                // Find the new point to insert based on the remaining length to retract
                double remaining_length = retract_length - retracted_length;
                double ratio = remaining_length / segment_length;
                
                Point<Inexact> new_point = interpolate(branches[j][i], branches[j][i + 1], ratio);
                remove_vertex(branches[j][i + 1]);
                remove_vertex(branches[j][i]);
                add_edge(new_point, branches[j][i + 1]);
                branches[j][i] = new_point; // Replace the point at i-1 with the new point
                last_index = i; 
                retracted_length += segment_length;
                 
                break;
            }
        }

        // remove the retracted part of the branch
        branches[j].erase(branches[j].begin(),branches[j].begin() + last_index);
    }
}

//compute the last new leaf
Point<Inexact> MedialAxis::interpolate(const Point<Inexact>& start, const Point<Inexact>& end, double ratio) {
    auto x = start.x() + ratio * (end.x() - start.x());
    auto y = start.y() + ratio * (end.y() - start.y());
    return Point<Inexact>(x, y);
}


void MedialAxis::apply_modified_negative_offset(double constant_offset, double min_radius) {
    // iterate through each vertex-radius pair in the radius_list.
    for (auto& vertex_radius_pair : radius_list) {
        // reduce the radius by the constant offset.
        std :: cout << "Applying modified negative offset..." << std::endl;
        double new_radius = vertex_radius_pair.second - constant_offset;
        
        // ensure the new radius does not fall under the min threshold.
        // if it does, set it to the mini theshold.
        if (new_radius < min_radius) {
            new_radius = min_radius;
        }
        // print old vs new radious
        std::cout << "Old radius: " << vertex_radius_pair.second << " New radius: " << new_radius << std::endl;
        vertex_radius_pair.second = new_radius;
    }
}

    MedialAxis::MedialAxis(const Polygon<Inexact>& shape) {
        assert(shape.is_counterclockwise_oriented());
        polygon = shape;
        // create interior straight skeleton
        // TODO: maybe make a separate function to just compute this skeleton,
        // and call it in the constructor
        // TODO: maybe even put this under the Polygon class?
        iss = CGAL::create_interior_straight_skeleton_2(shape);
        if (!iss) {
            std::cout << "Failed creating interior straight skeleton." << std::endl;
            exit(1);
        }
        std::cout << "Successfully computed interior straight skeleton." << std::endl;

        for (auto halfedge = iss->halfedges_begin(); halfedge != iss->halfedges_end(); halfedge++) {
            if (halfedge->is_bisector()) {
                Point<Inexact> t = halfedge->vertex()->point();
                Point<Inexact> s = halfedge->opposite()->vertex()->point();
                add_edge(s, t);
            }
        }
    }

} // namespace cartocrow::medial_axis
