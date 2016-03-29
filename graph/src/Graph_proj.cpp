//============================================================================
// Name        : graph_proj.cpp
// Author      : Prabhat Verma
// Version     : 0.1
// Copyright   : Fuck this shit!
// Description : Get a chain of words such that next word starts with an
//               alphabet that the current word ends in
//============================================================================

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cassert>
#include <stack>
#include <algorithm>
#include <set>

namespace MyGraph {

struct Node {
	typedef std::pair<long, Node*> Edge;
	std::vector<Edge> adjacencyList_;
	std::string name_;

	Node(const std::string& s) : name_(s) {
	}
};

struct Degree {
	size_t in_, out_;
	Degree() : in_(0), out_(0) {
	}
};

class Graph {
public:
	typedef std::map<std::string, Node*> Vertices;
	typedef std::map<std::string, Degree> Degrees;
	typedef std::map<std::string, bool> VisitedMap;

	void addVertex(const std::string& iName);
	void addEdge(const std::string& iFrom, const std::string& iTo, const long iCost);
	void display();
	const Vertices& getVertices() const;
	const Degrees& getDegrees() const;

private:
	Vertices vertices_;   // A map of vertices
	Degrees  degrees_;    // Maintain a map of in/out degrees of each vertex
};

const Graph::Degrees& Graph::getDegrees() const {
	return degrees_;
}

const Graph::Vertices& Graph::getVertices() const {
	return vertices_;
}

void Graph::display() {
	for(const auto& vertex : vertices_) {
		std::cout << "{{" << vertex.first << "\t"
				  << "In{" << degrees_[vertex.first].in_ << "} "
				  << "Out{" << degrees_[vertex.first].out_ << "} "
				  << "}} --> " << '\t';
		for(const auto& node : vertex.second->adjacencyList_) {
			std::cout << "{" << node.first << " ," << node.second->name_ << "}" << '\t';
		}
		std::cout << std::endl;
	}
}

void Graph::addVertex(const std::string& iName) {
	if(vertices_.find(iName) == vertices_.end()) {
		// Add a new vertex
		vertices_[iName] = (new Node{iName});
		// Add a new Degree element
		degrees_[iName] = Degree();
	}
//	else {
//		std::cout << "Vertex already exists!";
//	}
}

void Graph::addEdge(const std::string& iFrom, const std::string& iTo, const long iCost) {
	const auto& fromIt(vertices_.find(iFrom));
	const auto& toIt(vertices_.find(iTo));
	assert(fromIt != vertices_.end() || toIt != vertices_.end());
	vertices_[fromIt->first]->adjacencyList_.push_back(std::make_pair(iCost, toIt->second));
	// Update the in and out degree of respective vertices
	++degrees_[iFrom].out_;
	++degrees_[iTo].in_;
}


Graph transpose(const Graph& iGraph) {
	Graph transposedGraph;
	for(const auto& vertex : iGraph.getVertices()) {
		transposedGraph.addVertex(vertex.first);
		for(const auto& edge : vertex.second->adjacencyList_) {
			transposedGraph.addVertex(edge.second->name_);
			transposedGraph.addEdge(edge.second->name_, vertex.first, edge.first);
		}
	}
	return transposedGraph;
}

void DFSUtil(const Graph& iGraph, const std::string& iVertexName, Graph::VisitedMap& ioVisitedMap, std::stack<std::string>& oVisitedOrder) {
	// Mark node iVertexName as visited
	ioVisitedMap[iVertexName] = true;
	//std::cout << "Visited vertex --> " << iVertexName << std::endl;
	// For all neighbors of iVertexName
	//   If the neighbor is not already visited
	//     Call DFSUtil on that neighbor
	const Graph::Vertices& vertices = iGraph.getVertices();
	const std::vector<Node::Edge>& neighbors = vertices.find(iVertexName)->second->adjacencyList_;
	for(const auto& neighbor : neighbors) {
		if(!ioVisitedMap[neighbor.second->name_]) {
			DFSUtil(iGraph, neighbor.second->name_, ioVisitedMap, oVisitedOrder);
		}
	}

	oVisitedOrder.push(iVertexName);
}

void depthFirstSearch(const Graph& iGraph, const std::string& iVertexName, std::stack<std::string>& oVisitedOrder) {
	const Graph::Vertices& vertices = iGraph.getVertices();
	// Initialize the visited map of vertices
	Graph::VisitedMap visitedMap;
	for(const auto& vertex : vertices) {
		visitedMap.insert(std::make_pair(vertex.first, false));
	}
	DFSUtil(iGraph, iVertexName, visitedMap, oVisitedOrder);
}

size_t getConnectedComponents(const Graph& iGraph, std::vector<std::vector<std::string> >& ovvStrings) {
	size_t numConnectedComponents{0};
	const Graph::Vertices& vertices = iGraph.getVertices();
	// Initialize the visited map of vertices
	Graph::VisitedMap visitedMap;
	for(const auto& vertex : vertices) {
		visitedMap.insert(std::make_pair(vertex.first, false));
	}

	std::stack<std::string> visitedOrder;
	for(const auto& vertex : vertices) {
		if(!visitedMap[vertex.first]) {
			DFSUtil(iGraph, vertex.first, visitedMap, visitedOrder);
		}
	}

	Graph transposedGraph(transpose(iGraph));

	// Mark all vertices as not visited
	for(auto& element : visitedMap) {
		element.second = false;
	}

	while(!visitedOrder.empty()) {
		std::string vertex(visitedOrder.top());
		visitedOrder.pop();
		if(!visitedMap[vertex]) {
			std::stack<std::string> visitedOrderTransposed;
			DFSUtil(transposedGraph, vertex, visitedMap, visitedOrderTransposed);
			++numConnectedComponents;
			std::vector<std::string> vStrings;
			while(!visitedOrderTransposed.empty()) {
				vStrings.push_back(visitedOrderTransposed.top());
				visitedOrderTransposed.pop();
			}
			ovvStrings.push_back(vStrings);
		}
	}

	return numConnectedComponents;
}

void displayConnectedComponents(const std::vector<std::vector<std::string>>& ivvStrings) {
	size_t component(1);
	for(const auto& vStrings : ivvStrings) {
		std::cout << "component --> " << component << " \n";
		for(const auto& string : vStrings) {
			std::cout << string << "\t";
		}
		++component;
		std::cout << std::endl;
	}
}

auto isGraphEulerian(const Graph& iGraph) -> bool {
	long indegreeOutlier{0}, outdegreeOutlier{0}, degreeSameCount{0};
	std::string indegreeOutlierName, outdegreeOutlierName;
	const auto& degrees = iGraph.getDegrees();
	for(const auto& vertex : iGraph.getVertices()) {
		Degree degree = degrees.find(vertex.first)->second;
		if((degree.out_ - degree.in_) == 1) {
			++outdegreeOutlier;
			outdegreeOutlierName = vertex.first;
		} else if ((degree.in_ - degree.out_) == 1) {
			++indegreeOutlier;
			indegreeOutlierName = vertex.first;
		} else if (degree.in_ == degree.out_) {
			++degreeSameCount;
		} else {
			return false;
		}
	}

	// Do all vertices of non-zero degree form a single connected component?
	std::vector<std::vector<std::string> > vvStrings;
	if(1 != getConnectedComponents(iGraph, vvStrings)) {
		return false;
	}

	return true;
}

bool canAllNodesWithNonZeroDegreeBeVisited(const Graph iGraph) {
	const Graph::Vertices& vertices = iGraph.getVertices();
	const Graph::Degrees& degrees = iGraph.getDegrees();
	std::set<std::string> nonZeroDegreeVertices;
	for(const auto& degree : degrees) {
		if(0 != degree.second.in_ || 0 != degree.second.out_) {
			nonZeroDegreeVertices.insert(degree.first);
		}
	}

	std::set<std::string> verticesVisitedbyDFS;
	for(const auto& str : nonZeroDegreeVertices) {
		verticesVisitedbyDFS.clear();
		std::stack<std::string> visitedOrder;
		// Initialize the visited map of vertices
		Graph::VisitedMap visitedMap;
		for(const auto& vertex : vertices) {
			visitedMap.insert(std::make_pair(vertex.first, false));
		}
		DFSUtil(iGraph, str, visitedMap, visitedOrder);
		while(!visitedOrder.empty()) {
			verticesVisitedbyDFS.insert(visitedOrder.top());
			visitedOrder.pop();
		}
		std::vector<std::string> intersect;
		intersect.reserve(iGraph.getVertices().size());
		if(nonZeroDegreeVertices.size() == verticesVisitedbyDFS.size()) {
			// Sizes are equal. Check if all vertices are covered
			std::set_intersection(nonZeroDegreeVertices.begin(),
								  nonZeroDegreeVertices.end(),
								  verticesVisitedbyDFS.begin(),
								  verticesVisitedbyDFS.end(),
								  std::back_inserter(intersect)
							     );
			if(intersect.size() == nonZeroDegreeVertices.size()) {
				return true;
			}
		}
	}

	return false;
}

} // namespace Graph

int main(int argc, char** argv) {
	using Graph = MyGraph::Graph;

    size_t T(0), N(0);
	std::cin >> T >> N;
	for(size_t cnt(T); cnt > 0; --cnt) {
		Graph g;
		for(size_t strCnt(N); strCnt > 0; --strCnt) {
			std::string inputStr;
			std::cin >> inputStr;
			std::string startVertex(inputStr.substr(0,1));
			std::string endVertex(inputStr.size() > 1
							      ? inputStr.substr(inputStr.size()-1,1)
							      : inputStr.substr(0,1));
//			std::cout << "startVertex --> " << startVertex << "\t"
//					  << "endVertex --> " << endVertex << std::endl;
			g.addVertex(startVertex);
            g.addVertex(endVertex);
            g.addEdge(startVertex, endVertex, 1);
		}

//		std::cout << "Original --> \n";
//		g.display();
//		depthFirstSearch(g, "A");
//		Graph transposed(transpose(g));
//		std::cout << "Transposed --> \n";
//		transposed.display();

//		std::vector<std::vector<std::string> > vvStrings;
//		MyGraph::getConnectedComponents(g, vvStrings);
//		MyGraph::displayConnectedComponents(vvStrings);

		std::cout << (MyGraph::canAllNodesWithNonZeroDegreeBeVisited(g)
						? "Head to tail ordering is possible."
						: "Head to tail ordering is not possible.")
				  << std::endl;
	}
	return 0;
}
