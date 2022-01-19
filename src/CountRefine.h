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

	void WriteCnt()
	{
		if (opt.incnt != "")
		{
			Tokenizer in(opt.incnt);
			while (in.nextLine())
			{
				string f = in.getToken(); 
				int    c = stoi(in.getToken());
				D.FragCount[f] += c;
			}
		}	

		ofstream out(opt.outcnt);
		for (auto& it : D.FragCount)
			out << it.first << '\t' << it.second << endl;
	}

	void WriteAlt()
	{
		if (opt.inalt != "")
		{
			Tokenizer in(opt.inalt);
			while (in.nextLine())
			{
				RBT& tree = gtf.T[in.getToken()];
				Node *l   = tree.search(stoi(in.getToken()));
				Node *r   = tree.search(stoi(in.getToken()));
				string a  = in.getToken();
				string b  = in.getToken();
				int    c  = stoi(in.getToken());
				pair<Node*, Node*> P(make_pair(l, r));
				vector<string>& v = D.AltHeader[P];
				if (!v.size()) v = vector<string>(6, "");
				v[AltData::Str2EID(b)] = a;
				if (!D.AltCount[P]) D.AltCount[P] += c;
				D.AltGenes[P] = in.getToken(); 
			}
		}

		ofstream out(opt.outalt);
		for (auto& it : D.AltCount)
		{
			Node*  l = it.first.first;
			Node*  r = it.first.second;
			string g = D.AltGenes[it.first];
			char   k = 0;

			for (string& t : D.AltHeader[it.first])
			{
				if (t.size()) out << l->chr<< '\t' << l->e << '\t' << r->s << '\t'
					              << t << '\t' << AltData::EID2Str(k) << '\t' 
                                  << it.second << '\t'
                                  << g << endl;				
				++k; 
			}
		}

		for (auto& it : D.IntCount)
		{
			const Node* n = it.first;
			out << n->chr << '\t' << n->s << '\t' << n->e << "\t-\tIR\t" << it.second << endl; 
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
			pair<Node*, Node*> P(V[i].back(), V[i + 1].front());
			if (!D.IsAltP(P))
			{
				string genes = MT_Aligner::GetGenes(P.first, P.second); 
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

	void ParseAlignments()
	{
		for (int i = 0; i < RDS.size(); ++i)
		{
			Transcript T = Cigar2Frag(get<1>(RDS[i]), get<3>(RDS[i]), get<2>(RDS[i]));	
			D.AddFrag(MT_CountRefine::Trans2Frag(gtf, T));
			ProcessAlt(T);
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

	void NewNode(string id, string c, int s, int e, set<string>& t)
	{
		Node N;
		N.id = id;
		N.s = s;
		N.e = e;
		N.chr = c;
		if (id[0] != 'I') N.tran = set<string>(t);
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

                if (!i && s + v[i] < ne)
                {
                    t.remove(n);
					NewNode(NewId(nid), c, ns, s + v[i], nt);
					NewNode(NewId(nid), c, s + v[i] + 1, ne, nt);
                }
                else if (i == (int)v.size() && ns < s + 1)
                {
                    t.remove(n);
					NewNode(NewId(nid), c, ns, s, nt);
					NewNode(NewId(nid), c, s + 1, ne, nt);
                }
                else if (i && i < (int)v.size())
                {
                    t.remove(n);
                    if (ns < s + 1)
						NewNode(NewId(nid), c, ns, s, nt);
					NewNode(NewId(nid), c, s + 1, s + v[i], nt);
                    if (s + v[i] < ne)
						NewNode(NewId(nid), c, s + v[i] + 1, ne, nt);
                }
            }                
            else
            {
                if (i && ns < s + 1)
                {
					set<string> nt(n->tran);
                    t.remove(n);
                    NewNode(NewId(nid), c, ns, s, nt);
					NewNode(NewId(nid), c, s + 1, ne, nt);
		    		m = t.search(ms);
                }    
                if (i < (int)v.size() - 1 && s + v[i] + 1 < me)
                {
					set<string> mt(m->tran);
                    t.remove(m);
					NewNode(NewId(mid), c, ms, s + v[i], mt);
					NewNode(NewId(mid), c, s + v[i] + 1, me, mt);
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