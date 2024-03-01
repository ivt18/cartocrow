#include "medial_axis.h"
#include "../core/core.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <cmath>

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

    void MedialAxis::calculate_weight_function() {
        std::cout << "Trying to print pair" << std::endl;   
        for (const auto& vertex: graph) {
            const Point<Inexact>& current_point = vertex.first;
            double radius = INFINITY;
            Segment<Inexact> segment;
            Segment<Inexact> segment2;
            std::cout << "vertex " << current_point << std::endl;
            for (auto edgeIt = polygon.edges_begin(); edgeIt != polygon.edges_end(); ++edgeIt) {
                Segment<Inexact> edge = *edgeIt;
                Inexact::FT dist = CGAL::squared_distance(current_point, edge);
                if (dist < radius) radius = dist;
            }
            radius_list.insert(std::make_pair(current_point, radius));
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
    }

    AdjacencyList<Inexact> MedialAxis::temporary_remove_branch(int index) {
        AdjacencyList<Inexact> output_graph = graph;
        Branch<Inexact> branch = branches[index];
        if (branch.size() == 0) {
            // If the branch only contains the junction point or is empty, nothing to remove.
            return output_graph;
        }

        // Iterate through all points in the branch except the last one (the junction point).
        for (size_t i = 0; i < branch.size(); ++i) {
            const Point<Inexact>& s=  branch[i];
            // Remove s from the adjacency lists of its neighbors.
            for (const auto& neighbor : output_graph[s]) {
                auto& neighborsList = output_graph[neighbor];
                neighborsList.remove(s);
            }

            // Finally, remove the point from the output_graph.
            output_graph.erase(s);
        }
        return output_graph;
    }

    double MedialAxis::get_branch_weight(int index) {
        // Compute total area
        double totalArea = 0.0;
        for (const auto& v : graph) {
            const Point<Inexact>& point = v.first;
            /* double radius = get_radius(point); */
            double radius = 2.0;
            totalArea += M_PI * radius * radius;
        }

        // Compute area without brach
        double nonBranchArea = 0.0;
        for (const auto& v : temporary_remove_branch(index)) {
            const Point<Inexact>& point = v.first;
            /* double radius = get_radius(point); */
            double radius = 2.0;
            nonBranchArea += M_PI * radius * radius;
        }

        // return weight
        return (1.0 - nonBranchArea / totalArea);
    }

    void MedialAxis::prune_points(double t) {
        compute_branches();
        while (true) {
            double minWeight = INFINITY;
            int minIndex = -1;

            // Find the branch with the minimum weight
            for (size_t i = 0; i < branches.size(); ++i) {
                double weight = get_branch_weight(i);
                if (weight < minWeight) {
                    minWeight = weight;
                    minIndex = i;
                }
            }
            if (minIndex == -1 || minWeight >= t) {
                std::cout << "Pruning finished" << std::endl;
                return;
            }
            remove_branch(minIndex);
        }
    }

    void MedialAxis::compute_grid(unsigned int cells_x, unsigned int cells_y, unsigned int points_per_cell_edge) {
        Point<Inexact> top_left, bottom_right;
        CGAL::Bbox_2 bbox = polygon.bbox();
        top_left = Point<Inexact>(bbox.xmin(), bbox.ymax());
        bottom_right = Point<Inexact>(bbox.xmax(), bbox.ymin());
        double cell_width = (bottom_right.x() - top_left.x()) / cells_x;
        double cell_height = (top_left.y() - bottom_right.y()) / cells_y;
        for (unsigned int curr_y = 0; curr_y <= cells_y; curr_y++) {
            std::vector<Point<Inexact>> row(cells_x * points_per_cell_edge + 1);
            double delta_y = cell_height * curr_y;
            for (unsigned int curr_x = 0; curr_x < cells_x; curr_x++) {
                double delta_x = cell_width * curr_x;
                row[curr_x] = Point<Inexact>(top_left.x() + delta_x, top_left.y() - delta_y);
            }
            row.back() = Point<Inexact>(bottom_right.x(), top_left.y() - delta_y); // add the right-most point in the row
            grid.push_back(row);
        }
    }

    MedialAxis::MedialAxis(const Polygon<Inexact>& shape) {
        assert(shape.is_counter_clockwise_oriented());
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
