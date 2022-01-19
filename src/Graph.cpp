#ifndef GRAPH_CPP
#define GRAPH_CPP

#include "Graph.h"

bool Graph::insert_supstr(std::vector<std::string>& nodes)
{
	std::vector<Graph::Node*> L {root};
  bool flag = 0;

	for(auto& it : nodes)
		for(unsigned int i = 0; i < L.size(); ++i)
		{
			Graph::Node* node;
			
			if (L[i]->children.count(it))
				node = L[i]->children[it];
			else
			{
				node = new Graph::Node(it);
				L[i]->children[it] = node;
				flag = 1;
			}	
			
			if (i) 
				L[i] = node;	
			else 
				L.push_back(node);
		}
	
	return flag;
}

bool Graph::insert_substr(std::vector<std::string>& nodes)
{
	std::vector<Graph::Node*> L {root};

	for(auto& it : nodes)
	{
		unsigned int size = L.size();
		for(unsigned int i = 0; i < size; ++i)
		{
			Graph::Node* node;
			
			if (L[i]->children.count(it))
			{
				node = L[i]->children[it];
				if (node->children.count(" ")) return 0;
			}
			else
			{
				node = new Graph::Node(it, L[i]->depth + 1);
				L[i]->children[it] = node;
				if (node->depth == nodes.size())
					node->children[" "] = 0;
			}	
			
			if (i) 
				L[i] = node;	
			else 
				L.push_back(node);
		}
	}
	
	return 1;
}

void Graph::clear()
{
	std::vector<Graph::Node*> L {root};
	
	while(L.size())
	{
		Graph::Node* n = L.back();
		L.pop_back();
		
		if (!n) continue;
		
		for(auto& it : n->children)
			L.push_back(it.second);
			
		delete n;
	}
	
	root = new Graph::Node(" ", 0);
}

Graph::~Graph()
{
	std::vector<Graph::Node*> L {root};
	
	while(L.size())
	{
		Graph::Node* n = L.back();
		L.pop_back();
		
		if (!n) continue;
		
		for(auto& it : n->children)
			L.push_back(it.second);
			
		delete n;
	}
}
#endif
