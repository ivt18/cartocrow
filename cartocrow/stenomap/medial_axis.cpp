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

    MedialAxis::MedialAxis(const Polygon<Inexact>& shape) {
        assert(shape.is_counter_clockwise_oriented());
        
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

        for (auto halfedge = iss->halfedges_begin(); halfedge != iss->halfedges_end(); ++halfedge) {
            if (halfedge->is_bisector()) {
                Point<Inexact> t = halfedge->vertex()->point();
                Point<Inexact> s = halfedge->opposite()->vertex()->point();
                add_edge(s, t);
            }
        }
    }

} // namespace cartocrow::medial_axis
