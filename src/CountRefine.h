#ifndef COUNTREFINE_H
#define COUNTREFINE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <vector>
#include <set>
#include <map>
#include <tuple>    
#include <algorithm>
#include <thread>

#include <htslib/sam.h>

#include "RBT.h"
#include "BamParser.h"

using namespace std;

class CountRefineOpt
{
	public:

	string bam;
	string gtf;

	string incnt, outcnt;
	string inalt, outalt;

	string ref;

	int rl;
};

class MT_CountRefine
{
	public:
	
	MT_CountRefine(CountRefineOpt& opt) : opt(opt), unq(0), gtf(opt.gtf)//, CNT(0)
	{
		cerr << "[cntrf] generating intron data ... " << flush;
		gtf.genIntrons();
		cerr << "done" << endl;
		cerr << "[cntrf] processing BAM file ... " << flush;
		ParseBam();
		cerr << "done" << endl;
		cerr << "[cntrf] generating transcript data ... " << flush;
		gtf.genTranscripts();
		ComputeOriginalNodes();
		cerr << "done" << endl;
		cerr << "[cntrf] processing alignments ... " << flush;
		ParseAlignments();
		cerr << "done" << endl;
		cerr << "[cntrf] aligned " << D.num_mapped << " additional reads" << flush; 
		if (opt.outalt != "")
			cerr << ", " << D.AltHeader.size() << " novel events" << endl;
		else
			cerr << endl;
		cerr << "[cntrf] augmenting/writing outputs ... " << flush;
		if (opt.outcnt != "") WriteCnt();
		if (opt.outalt != "") WriteAlt();
		if (opt.ref    != "") WriteRef();
		cerr << "done" << endl;	
	}

	private:

	void ComputeOriginalNodes()
	{
		for (auto& it : gtf.T)
		{
			RBT& t = it.second;
			for (Node* n = t.min(t.getRoot()); n && n != t.getNil(); n = t.successor(n))
			{
				if (n->id[0] == 'I') continue;
				string bid = BaseId(n->id);
				if (bid != n->id)
				{
					OriginalNodes[bid].push_back(n);
				}
			}
		}
	}
	
	void WriteCnt()
	{
		if (opt.incnt != "")
		{
			Tokenizer in(opt.incnt);
			while (in.nextLine())
			{
				string frag = in.getToken();
				Transcript T = Frag2Trans(frag);

				if (frag != Trans2Frag(gtf, T))
				{
					continue;
				}

				while(true)
				{
					string t = in.getToken();
					if (t == "") break;
					int s = stoi(t);

                    Transcript F = MT_CountRefine::TrimTrans(T, s, opt.rl);
					++D.FragCount[Trans2Frag(gtf, F)].count;		
				}
			}
		}	

		ofstream out(opt.outcnt);
		for (auto& it : D.FragCount)
			out << it.first << '\t' << it.second.count << '\n';
	}

	void WriteAlt()
	{
		map<tuple<string, int, int, string>, tuple<string, string, int>> M;

		if (opt.inalt != "")
		{
			Tokenizer in(opt.inalt);
			while (in.nextLine())
			{
				//string line = in.getLine();
                //int p = line.find_last_of('\t');
				string chr = in.getToken();
				int s = stoi(in.getToken());
                int e = stoi(in.getToken());
				string t = in.getToken();
                string g = in.getToken();
				string id = in.getToken();
				int c = stoi(in.getToken());
                //int c = stoi(line.substr(p + 1, line.size() - p - 1));
				//M[make_tuple(chr, s, e)] = c;
				M[make_tuple(chr, s, e, id)] = make_tuple(t, g, c);
			}
		}

		for (auto& it : D.AltCount)
		{
			Node*  l = it.first.first.back();
			Node*  r = it.first.second.front();
			string g = D.AltGenes[it.first];
			char   k = 0;

			for (string& t : D.AltHeader[it.first])
			{
				if (t.size()) 
				{
					tuple<string, string, int>& m = M[make_tuple(l->chr, l->e, r->s, AltData::EID2Str(k))];
					get<2>(m) += it.second;
					if (get<0>(m) == "")
					{
						get<0>(m) = t;
						get<1>(m) = g == "" ? "-" : g;
					}
				}
				++k;
			}
		}

		ofstream out(opt.outalt);

		for (auto& it : M)
		{
			out << get<0>(it.first) << '\t'
                << get<1>(it.first) << '\t'
                << get<2>(it.first) << '\t'
                << get<0>(it.second) << '\t' 
                << get<1>(it.second) << '\t'
                << get<3>(it.first) << '\t'
                << get<2>(it.second) << '\n';
		}

		for (auto& it : D.IntCount)
		{
			const Node* n = it.first;
			out << n->chr << '\t' << n->s << '\t' << n->e << "\t-\t-\tIR\t" << it.second << endl; 
		}
	}

	void WriteRef()
	{
		ofstream out(opt.ref);
		for (auto& it : gtf.T)
		{
			RBT& t = it.second;
			for (Node* n = t.min(t.getRoot()); n && n != t.getNil(); n = t.successor(n))
				if (n->id != BaseId(n->id) || n->id[0] == 'I')
					out << n->id << '\t' << it.first << '\t' << n->s << '\t' << n->e << endl;
		}
	}

	Transcript Cigar2Frag(string chr, string cig, int s)
	{
		vector<int> v = ParseCigarC(cig);
		RBT& tree = gtf.T[chr];
		Transcript T;
		
		for (int i = 0; i < (int)v.size(); s += v[i++])   
        {
            if (i % 2) continue;
			Node *l = tree.search(s + 1);
			Node *r = tree.search(s + v[i]);
			while (1)
			{
				T.nodes.push_back(l);
				if (l->s == r->s) break;
				l = tree.successor(l);
			}
		}
		return T;
	}

	void ProcessAlt(Transcript& T)
	{
		vector<vector<Node*>> V({vector<Node*>({T.nodes[0]})});

		for (int i = 1; i < T.nodes.size(); ++i)
			if (T.nodes[i]->s - T.nodes[i - 1]->e > 1)
				V.push_back(vector<Node*>({T.nodes[i]}));
			else
				V.back().push_back(T.nodes[i]);

		for (int i = 0; i < V.size() - 1; ++i)
		{
			pair<vector<Node*>, vector<Node*>> P(V[i], V[i + 1]);
			if (!D.IsAltP(P))
			{
				string genes = MT_Aligner::GetGenes(P.first.back(), P.second.front()); 
				bool f = false;
				for (auto& a : AltSeq(gtf, V[i], V[i + 1]).GetEvents())
				{
					D.AddAltH(P, AltData::Str2EID(a.second), a.first, genes);
					f = true;
				}
				D.AddAltP(P, f);
			}
			if (D.GetAltP(P)) D.AddAltC(P);
		}

		if (V.size() == 1)
        {
        	vector<Node*>& W = V[0];
        	if (W.size() == 1 && W[0]->id[0] == 'I' && W[0]->id == BaseId(W[0]->id))
        		D.AddInt(W[0]);
        }
        else for (int i = 0; i < V.size(); ++i)
        {
            vector<Node*>& W = V[i];
            if (W.size() < 2) continue;
            for (int j = 0; j < W.size(); ++j)
            {				
                if (W[j]->id[0] != 'I') continue;
				if (D.IntExists(W[j]))
					D.AddInt(W[j]);
                else if (!i && !j && W[j + 1]->id[0] != 'I')
                    D.AddInt(W[j]);
                else if (i == V.size() - 1 && j == W.size() - 1 && W[j - 1]->id[0] != 'I')
                    D.AddInt(W[j]);
    	        else if (j && j < W.size() - 1 && W[j - 1]->id[0] != 'I' && W[j + 1]->id[0] != 'I')
                    D.AddInt(W[j]);
            }
        }
	}

	static Transcript TrimTrans(Transcript& T, int s, int rl)
	{
		Transcript N;
		for (int i = 0, c = 0; i < T.nodes.size(); ++i)
		{
			c += T.nodes[i]->size();
			if (s + 1 <= c) N.nodes.push_back(T.nodes[i]); // s <= c
			if (s + rl <= c) break; // <=
		}
		return N;
	}

	Transcript Frag2Trans(string frag)
	{
		vector<string> v = SplitString(frag, '|');
		Transcript T;
		for (int i = 0; i < v.size(); ++i)
		{
			if (OriginalNodes.find(v[i]) == OriginalNodes.end())
			{
		        T.insert(gtf.getNode(v[i]));
			}
			else	
			{
				vector<Node*>& o = OriginalNodes[v[i]];
				for (int j = 0; j < o.size(); ++j)
				{
					T.insert(o[j]);
				}
			}
		}
		T.sort();
		return T;
	}

	string Trans2Frag(GTF& G, Transcript& T)
	{
		string s = "";
		for (int i = 0; i < T.nodes.size(); ++i)
		{
			Node *n = T.nodes[i];
			string bid = BaseId(n->id);

			if (bid == n->id || bid[0] == 'I')
				s += n->id + "|";
			else 
			{
				vector<Node*>& o = OriginalNodes[bid];
				bool f = i + o.size() <= T.nodes.size();
				for (int k = 0; f && k < o.size(); ++k)
				{
					f = T.nodes[i + k]->id == o[k]->id;
				}
				if (f)
				{
					s += bid + "|"; 
                    i += o.size() - 1;
				}
				else
					s += n->id + "|";				
			}
		}
		return s;
	}

	/*
	static string Trans2Frag(GTF& G, Transcript& T)
	{
		RBT& t = G.T[T.nodes[0]->chr];
		string s = "";
		for (int i = 0; i < T.nodes.size(); ++i)
		{
			Node *n = T.nodes[i];
			string bid = BaseId(n->id);

			if (bid == n->id || (i && bid == BaseId(t.predecessor(n)->id)) || bid[0] == 'I')
				s += n->id + "|";
			else 
			{
				bool f = true;
				int k = i + 1;
				while (k < T.nodes.size())
				{
					Node *m = T.nodes[k];
					if (BaseId(m->id) != bid) break;
					f = t.predecessor(m) == T.nodes[k - 1];
					if (!f) break;
					++k;
				}
				if (f && (BaseId(t.successor(T.nodes[k - 1])->id) != bid || k == T.nodes.size()))
					s += BaseId(n->id) + "|", i = k - 1;
				else
					s += n->id + "|";				
			}
		}
		return s;
	}
	*/

	void ParseAlignments()
	{
		for (int i = 0; i < RDS.size(); ++i)
		{
			Transcript T = Cigar2Frag(get<1>(RDS[i]), get<3>(RDS[i]), get<2>(RDS[i]));
			
			if (opt.outcnt != "")
			{	 
				D.AddFrag(MT_CountRefine::Trans2Frag(gtf, T), 0);
			}

			if (opt.outalt != "")
			{
				ProcessAlt(T);
			}
		}
	}

	void ParseBam()
	{
		IN = hts_open(opt.bam.c_str(), "r");
		HD = sam_hdr_read(IN);
		A = bam_init1();
		while(sam_read1(IN, HD, A) > 0) ParseLine();
		bam_destroy1(A);
		sam_close(IN);
	}

	void ParseLine()
	{
		if (A->core.flag != 0 && A->core.flag != 16) return;
		string id = bam_get_qname(A);
		string chr = HD->target_name[A->core.tid];
		int pos = A->core.pos;	
		string cig = ParseCigarB(bam_get_cigar(A), A->core.n_cigar);
		if (cig == "") return;
		vector<int> C = ParseCigarC(cig);
		if (!Valid(C, chr, pos)) return;
        if (C.size() > 1) RefineNodes(C, chr, pos);
		RDS.emplace_back(id, chr, pos, cig);
	}

    vector<int> ParseCigarC(string c)
    {
        vector<int> v; string s = "";
        bool f = false;

        for (int i = 0; i < (int)c.size(); ++i)
            if (isdigit(c[i]))
                s += c[i];
            else
            {
                if (c[i] != 'M' && c[i] != 'N') 
                    return vector<int>();
                f = f || c[i] == 'N';
                v.push_back(stoi(s));
                s = "";                
            }

        return f ? v : vector<int>();
    }

	string ParseCigarB(uint32_t* c, uint32_t n)
    {
        string s = ""; bool f = false;

        for (uint32_t i = 0; i < n; ++i)
		{
			uint32_t op = bam_cigar_op(c[i]); 
			if (op != 0 && op != 3) return "";
			f = f || op == 3;
			s += to_string(bam_cigar_oplen(c[i])) + bam_cigar_opchr(c[i]);
		}

        return f ? s : "";
    }

	bool Valid(vector<int>& v, string c, int s)
	{
		RBT& t = gtf.T[c];
		for (int i = 0; i < (int)v.size(); s += v[i++])   
			if (i % 2) 
				continue;
			else if (t.search(s + 1) == t.getNil() || t.search(s + v[i]) == t.getNil())
				return false;
		return true;
	}

	void NewNode(string id, string c, int s, int e, set<string>& t, set<string>& g)
	{
		Node N;
		N.id = id;
		N.s = s;
		N.e = e;
		N.chr = c;

		if (id[0] != 'I') 
		{
			N.gene = set<string>(g);
			N.tran = set<string>(t);
		}

		gtf.T[c].insert(N);
		gtf.C[id] = make_pair(c, s); 		
	}

    void RefineNodes(vector<int>& v, string c, int s)
    { 
        RBT& t = gtf.T[c];

        for (int i = 0; i < (int)v.size(); s += v[i++])   
        {
            if (i % 2) continue;

            Node* n = t.search(s + 1);
            Node* m = t.search(s + v[i]);

            if (n->s == s + 1 && m->e == s + v[i])
                continue;

            string nid = BaseId(n->id);
            int ns = n->s; int ne = n->e;
			string chr = n->chr;
			string mid = BaseId(m->id);
            int ms = m->s; int me = m->e;

            if (ns == ms)
            {
				set<string> nt(n->tran); 
				set<string> gt(n->gene); 

                if (!i && s + v[i] < ne)
                {
                    t.remove(n);
					NewNode(NewId(nid), c, ns, s + v[i], nt, gt);
					NewNode(NewId(nid), c, s + v[i] + 1, ne, nt, gt);
                }
                else if (i == (int)v.size() && ns < s + 1)
                {
                    t.remove(n);
					NewNode(NewId(nid), c, ns, s, nt, gt);
					NewNode(NewId(nid), c, s + 1, ne, nt, gt);
                }
                else if (i && i < (int)v.size())
                {
                    t.remove(n);
                    if (ns < s + 1)
						NewNode(NewId(nid), c, ns, s, nt, gt);
					NewNode(NewId(nid), c, s + 1, s + v[i], nt, gt);
                    if (s + v[i] < ne)
						NewNode(NewId(nid), c, s + v[i] + 1, ne, nt, gt);
                }
            }                
            else
            {
                if (i && ns < s + 1)
                {
					set<string> nt(n->tran);
                    set<string> gt(n->gene);
                    t.remove(n);
                    NewNode(NewId(nid), c, ns, s, nt, gt);
					NewNode(NewId(nid), c, s + 1, ne, nt, gt);
		    		m = t.search(ms);
                }    
                if (i < (int)v.size() - 1 && s + v[i] + 1 < me)
                {
					set<string> mt(m->tran);
                    set<string> gt(n->gene); 
                    t.remove(m);
					NewNode(NewId(mid), c, ms, s + v[i], mt, gt);
					NewNode(NewId(mid), c, s + v[i] + 1, me, mt, gt);
                }
            }   
        }
    }

    static string BaseId(string id)
    {
        return id.substr(0, id.find("."));
    }

    string NewId(string id)
    {
        return id + "." + to_string(++unq);
    }

	private:
	
	map<string, vector<Node*>> OriginalNodes;
	vector<tuple<string, string, int, string>> RDS;
	CountRefineOpt& opt;
	GTF gtf;

	AltData D;

	samFile*   IN;
	bam_hdr_t* HD;
	bam1_t*    A;	

	int unq;
};

void CountRefine(CountRefineOpt& opt);

#endif
