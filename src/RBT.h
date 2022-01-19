#ifndef RBT_H
#define RBT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <vector>
#include <set>
#include <map>
#include <tuple>    
#include <algorithm>

using namespace std;

// Data container for the GTF subexon.
class Node
{
	public:

    // NodeId.
	string id;

	// Starting and ending coordinate.
	int s, e;

	// Location data.
	string chr;
	set<string> gene, tran;	

    // RBT stuff.
    Node *l, *r, *p;
    bool color;

	Node () : l(0), r(0), p(0), color(0)
	{ }

	Node(string id, int s, int e, Node* p) :
		id(id), s(s), e(e), chr(""), l(0), r(0), p(p), color(0)
	{ }

	Node(const Node& N) :
		id(N.id), s(N.s), e(N.e), chr(N.chr), gene(N.gene), tran(N.tran), l(0), r(0), p(N.p), color(0)
	{ }	

	// The amount of bases in the subexon.
	int size()
	{
		return e - s + 1;
	}
	
	// Some operators for easier comparison.
	bool operator<(const Node& rhs) const
	{
		return s < rhs.s;
	}
	
	bool operator<=(const Node& rhs) const
	{
		return s <= rhs.s;
	}
	
	bool operator==(const Node& rhs) const
	{
		return s == rhs.s;
	}

    void operator=(Node& rhs)
    {
        id = rhs.id;
        s = rhs.s, e = rhs.e;
    }

    bool contains(int x)
    {
        return s <= x && x <= e;
    }
};

class Transcript
{
	public:

	string id;
	vector<Node*> nodes;

	Transcript() : id("") { }

	void insert(Node* n)
	{
		if (find(nodes.begin(), nodes.end(), n) == nodes.end())
			nodes.push_back(n);
	}

	void sort()
	{
		std::sort(nodes.begin(), nodes.end(), 
			[](const Node* l, const Node* r)
		    { 
		    	return *l < *r; 
			}
		);
	}

	// Transcriptomic to genomic coordinates.
	// x - starting position within the transcript.
	// rev - is transcript on reverse strand?
	// trim - is transcript trimmed?
	int t2g(int x, bool rev, bool trim, int rl)
	{
		// Sorting nodes by their starting coordinates.
		vector<Node*>& v = nodes;
		
        // Trim left/right.
        trim = trim && v.size() > 1;
        bool tl = trim && v.front()->size() >= rl; 
        bool tr = trim && v.back()->size() >= rl;

		// Offseting start for reverse strand.
        if (rev) 
        {
            int c = 0;
            for (int i = 0; i < (int)v.size(); ++i)
                c += v[i]->size();
            x = c - (tl + tr + 1) * (rl - 1) - (x - 1);		
        }

		// Applying trimming criteria if applicable.
		if (tl) x += v.front()->size() - rl + 1;
			
		// Iterates through subexons in a transcript, finding to which one
		// does starting position x belong to.
		for (int i = 0; i < (int)v.size(); ++i)
			if (x <= v[i]->size()) 
				return max(-1, v[i]->s + x - 1);
            else
			    x -= v[i]->size();

		// x not found, default to -1.
		return -1;
	}

    int g2t(int x)
    {
        vector<Node*>& v = nodes;

        int c = 0;
        for (int i = 0; i < (int)v.size(); ++i)
            if (v[i]->s <= x && x <= v[i]->e)
                return c + x - v[i]->s + 1;           
            else
                c += v[i]->size();

        return -1;
    }
    
};

class RBT
{
    public:

    RBT() : nil(new Node("nil", 0, 0, 0)), root(nil) { }

    void insert(string id, int s, int e); 
	void insert(Node n);

    Node* search(Node* ptr, int x); 
    Node* search(int x) { return search(root, x); } 

    void remove(Node* ptr); 

    Node* min(Node* ptr); 
    Node* max(Node* ptr); 

    Node* successor(Node* ptr); 
    Node* predecessor(Node* ptr); 
    
    void clear(Node* ptr); 
    void clear() { clear(root); } 

    void inorder(Node* ptr); 
    void inorder() { inorder(root); } 

	Node* getRoot() { return root; }
	Node* getNil()  { return nil;  }

    ~RBT() { clear(root); delete nil; } 

    private:
    
    void rotateL(Node* x); 
    void rotateR(Node* y); 
    
    void insertRB(Node* x); 
    void removeRB(Node* x); 

    Node *nil, *root;
};

// GTF parser/container.
class GTF
{
	public:

	GTF(string filename)
	{
		if (filename == "") return;
		in.open(filename);
		while(getLine()) parseLine();
		line.clear();	
		in.close();			
	}

	Node* getNode(string chromosome, int x)
	{
        RBT& t = T[chromosome];
        Node* n = t.search(x);
        return n == t.getNil() ? 0 : n;
	}

	Node* getNode(string id)
	{
		auto it = C.find(id);
		if (it == C.end()) return 0;
        RBT& t = T[it->second.first];
        Node* n = t.search(it->second.second);
        return n == t.getNil() ? 0 : n;
	}

	string getChr(Node* n)
	{
		return C[n->id].first;
	}

	Transcript& getTrans(string id)
	{
		return S[id];
	}

	void genIntrons()
	{
		int k = 0;
		for (auto& it : T)
		{
			RBT& t = it.second;
			Node *n = t.min(t.getRoot()), *m;
			while (n && n != t.getNil())
			{
				m = t.successor(n);
				if (m->s - n->e > 1)
				{
					string c = n->chr;
					Node I;
					I.id  = "I" + to_string(++k);
					I.s   = n->e + 1;
					I.e   = m->s - 1;
					I.chr = it.first;
					t.insert(I);
					C[I.chr] = make_pair(I.chr, I.s);
					n = t.search(I.e + 1);
				}
				else 
					n = m;
			}
		}	
	}

	void genTranscripts()
	{
		for (auto& it : T)
		{
			RBT& t = it.second;
			for (Node* n = t.min(t.getRoot()); n && n != t.getNil(); n = t.successor(n))
			{
				if (n->id[0] == 'I') continue;
				for (auto& tr : n->tran)
				{
					Transcript& obj = S[tr]; 
					obj.insert(n);
					obj.id = tr;
				}
			}
		}			
	}

    map<string, RBT> T;
	map<string, Transcript> S;
	map<string, pair<string, int>> C;

	private:

	bool getLine()
	{
		for (string l; getline(in, l); )
		{
			l = l.substr(0, l.find("#"));			
			if (!l.size()) continue;
			line.clear();
			line.str(l);
			return true;
		}
		return false;
	}

	bool nextToken()
	{
		if (line >> token)
		{
			if (token.back() == ';')
				token = token.substr(0, token.size() - 1);
			if (token.back() == '\"')
				token = token.substr(1, token.size() - 2);						
			return true;
		}
		return false;
	}

	void parseLine()
	{
		string id = "", c = "", t = "", g = "";
		int s = -1, e = -1;

		for(int i = 0; nextToken(); ++i)
		{
			if (i == 0) 
				c = token; // chromosome
			else if (i == 2 && token != "subexon")
				return;	
			else if (i == 3)
				s = stoi(token);
			else if (i == 4)
				e = stoi(token);
			else if (token == "NodeId")
				nextToken(), id = token;
			else if (token == "transcript_id")
				nextToken(), t = token;
			else if (token == "gene_id")
				nextToken(), g = token;
		}

		int o = overlap(c, s, e);

        if (o == -1)
		{
			Node N;
			N.id  = id;
			N.s   = s;
			N.e   = e;
			N.chr = c;
			N.gene.insert(g);
			N.tran.insert(t);
            T[c].insert(N);
			C[id] = make_pair(c, s);
		}
		else 
		{
			Node* N = T[c].search(o);
			if (N->id == id) N->tran.insert(t);
			else cout << "[ warn] GTF node " << id << " intersects with " << N->id << endl; 
		}
	}

    int overlap(string c, int s, int e) // Think of something to do this in a better way.
    {
        for (RBT& t = T[c]; s <= e; ++s)
            if (t.search(s) != t.getNil())
                return s;          
        return -1;
    }

	string token;
	istringstream line;
	ifstream in;
};

#endif
