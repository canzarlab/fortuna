#ifndef GRAPH_H
#define GRAPH_H

#include<iostream>
#include<map>
#include<vector>
#include<list>

class Graph
{
	public:

	class Node
	{
		public:
		
		std::string id;
		unsigned int depth;
	  std::map<std::string, Node*> children;		
		
		Node(std::string id) : id(id) { }
		Node(std::string id, unsigned int depth) : id(id), depth(depth) { }
	};
	
	bool insert_supstr(std::vector<std::string>& nodes);	
	bool insert_substr(std::vector<std::string>& nodes);	
	
	void clear();
	
	Graph() { root = new Graph::Node(" ", 0); }
	~Graph();
	
	private:
	
	Graph::Node* root;
};

#endif
