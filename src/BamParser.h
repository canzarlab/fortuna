/*

	time ~/MT/MT/build/mt --quant -cnt cnt.txt -rl 100 -thr 7 -gtf GRCh38.nrg.gtf -ind 11 -infq r.fq -miss 8 -infa ../hg38.fa -outfq r.reduced.fq

*/

#ifndef BAMPARSER_H
#define BAMPARSER_H

#include <stdio.h>
#include <stdlib.h>
#include <htslib/sam.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <map>
#include <tuple>
#include <mutex>
#include <vector>
#include <thread>

#include "RBT.h"
#include "SeqUtils.h"
#include "ZParser.h"

using namespace std;

typedef unsigned short thread_t;

class BamParserOpt
{
	public:

	string infq, infa;
	string outfq, cnt, alt;
	string gtf;
	string ind;
	string bam;

	htsFile* bamfp;

	int (*cb)(int, char**, void*);

	int thr;
	int rl;
	int n;

	BamParserOpt() : 
		infq(""), infa(""),
		outfq(""), cnt(""), alt(""),
		gtf(""), ind(""),
		thr(1), rl(0), n(0), bam("")
	{ }
};

class Tokenizer
{
	public:

	Tokenizer(string filename)
	{
		in.open(filename);		
	}

	string getToken()
	{
        return (ss >> token) ? token : "";
	}

	string getLine()
	{
		return line;
	}

    bool nextLine()
    {
        while(getline(in, line))
		{		
			if (!line.size()) continue;
			ss.clear();
			ss.str(line);
			return true;
		}
		return false;
    }

    private:

	string token, line;
	istringstream ss;
	ifstream in;
};

class AltSeq
{
	public:

	AltSeq(GTF& G, vector<Node*>& L, vector<Node*>& R) : 
		G(G), L(L), R(R) 
	{ }

	vector<pair<string, string>> GetEvents()
	{
		vector<pair<string, string>> V;
		if (NovelIntron()) ProcessEvents(V);
		return V;
	}

	private:

	bool NovelIntron()
	{
		Node *l = L.back(), *r = R.front();
		if (l->id[0] == 'I' || r->id[0] == 'I') return true;
		vector<string> T(SetIntersection(l->tran, r->tran));
		for (int i = 0; i < T.size(); ++i)
		{
			vector<Node*>& V = G.getTrans(T[i]).nodes;
			if (Search(V, 0, V.size(), l) == Search(V, 0, V.size(), r) - 1)
				return false;
		}
		return true;
	}

	bool IsIntronic()
	{
		bool l = false, r = false;
		for (int i = 0; !l && i < L.size(); ++i)
			if  (L[i]->id[0] != 'I') l = true;
		if (!l) return true;
		for (int i = 0; !r && i < R.size(); ++i)
			if  (R[i]->id[0] != 'I') r = true;
		return !r;
	}

	vector<string> GetTrans()
	{
		set<string> A, B;
		for (int i = 0; i < L.size(); ++i)
			A.insert(L[i]->tran.begin(), L[i]->tran.end());
		for (int i = 0; i < R.size(); ++i)
			B.insert(R[i]->tran.begin(), R[i]->tran.end());
		return SetIntersection(A, B);
	}		

	int LFlank(vector<Node*>& V)
	{
		int k = -1;		
		for (int i = L.size() - 1; k == -1 && i >= 0; --i)
			k = Search(V, 0, V.size(), L[i]);		
		if (k > -1) 
			while(k < V.size() - 1 && V[k + 1]->s - V[k]->e == 1 && V[k + 1]->id[0] != 'I') ++k; 
		return k;
	}	

	int RFlank(vector<Node*>& V)
	{
		int k = -1;		
		for (int i = 0; k == -1 && i < R.size(); ++i)
			k = Search(V, 0, V.size(), R[i]);
		if (k > -1)
			while(k && V[k]->s - V[k - 1]->e == 1 && V[k - 1]->id[0] != 'I') --k; 
		return k;
	}	

	int Search(vector<Node*> V, int s, int e, Node* x)
	{
		if (!(e - s)) return -1;
		int m = s + (e - s) / 2;
		if (V[m]->s == x->s) return m;
		if (V[m]->s > x->s) return Search(V, s, m, x);
		return Search(V, m + 1, e, x);
	}

	vector<string> SetIntersection(set<string>& A, set<string>& B)
	{
		vector<string> R;
		set_intersection(A.begin(), A.end(), B.begin(), B.end(), back_inserter(R));
		return R;
	}

	void ProcessEvents(vector<pair<string, string>>& V)
	{
		vector<string> TS(GetTrans());
		if (!TS.size()) return;

		Node *l = L.back(), *r = R.front();
		string es = "", ie = "", ad = "", aa = "", ap = "", un = "";

		if (!IsIntronic()) for (auto& t : TS)
		{
			vector<Node*>& W = G.getTrans(t).nodes;
			int a = LFlank(W), b = RFlank(W);

			if (a == -1 || b == -1) continue;

			if (b < a && a - b > 1)
				ie += t + ",";
			else if (W[a] != l && W[b] == r && b - a == 1) 
				ad += t + ",";
			else if (W[a] == l && W[b] != r && b - a == 1)
				aa += t + ",";
			else if (W[a] != l && W[b] != r && b - a == 1)
				ap += t + ",";
			else if (b - a > 1 && W[a + 1]->s - l->e > 1 && r->s - W[b - 1]->e > 1)
				es += t + ",";
		}
		else
			un = "-";

		if (es.size())
			V.emplace_back(es.substr(0, es.size() - 1), "ES");	
		if (ad.size())
			V.emplace_back(ad.substr(0, ad.size() - 1), "AD");
		if (aa.size())
			V.emplace_back(aa.substr(0, aa.size() - 1), "AA");
		if (ap.size())
			V.emplace_back(ap.substr(0, ap.size() - 1), "AP");
		if (ie.size())
			V.emplace_back(ie.substr(0, ie.size() - 1), "IE");
		if (un.size())
			V.emplace_back(un, "XX");
	}

	GTF& G;
	vector<Node*> &L, &R;
};

class FragRecord
{
	public:

	vector<int> positions;
	int count;
};

class AltData
{
	public:

	map<string, FragRecord> FragCount;
	map<string, int> ReadCount;

	map<pair<vector<Node*>, vector<Node*>>, vector<string>> AltHeader;
	map<pair<vector<Node*>, vector<Node*>>, int> AltCount;
	map<pair<vector<Node*>, vector<Node*>>, bool> AltProcessed;
	map<pair<vector<Node*>, vector<Node*>>, string> AltGenes;
	
	map<Node*, int> IntCount;

	int num_unmapped;
	int num_mapped;
	int num_novel;

	AltData() : num_mapped(0), num_unmapped(0), num_novel(0) { }

	void AddFrag(string f, int s)
	{
		lock_guard<mutex> lock(writer_lock);
		FragRecord& R = FragCount[f];
		R.positions.push_back(s);
		++R.count;
		++num_mapped;
	}

	void AddRead(string r, int e)
	{
		lock_guard<mutex> lock(writer_lock);
		ReadCount[r] = e;
		++num_unmapped;
	}
	
	bool HeaderExists(pair<vector<Node*>, vector<Node*>>& I)
	{
		lock_guard<mutex> lock(writer_lock);
		return AltHeader.find(I) != AltHeader.end();
	}

	bool IsAltP(pair<vector<Node*>, vector<Node*>>& I)
    {
        lock_guard<mutex> lock(writer_lock);
        return AltProcessed.find(I) != AltProcessed.end();
    }

	bool GetAltP(pair<vector<Node*>, vector<Node*>>& I)
	{
		lock_guard<mutex> lock(writer_lock);
        return AltProcessed[I];
	}

	void AddAltP(pair<vector<Node*>, vector<Node*>>& I, bool b)
    {
        lock_guard<mutex> lock(writer_lock);
        AltProcessed[I] = true;
    }

	void AddAltC(pair<vector<Node*>, vector<Node*>>& I)
	{
		lock_guard<mutex> lock(writer_lock);
		if (AltHeader.find(I) != AltHeader.end())
			++AltCount[I], ++num_novel;
	}

	void AddAltH(pair<vector<Node*>, vector<Node*>>& I, char eid, string t, string g)
	{
		lock_guard<mutex> lock(writer_lock);
		vector<string>& v = AltHeader[I];
		if (!v.size()) v = vector<string>(6, "");
		v[eid] = t;
		AltGenes[I] = g;
	}

	bool IntExists(Node* I)
	{
		lock_guard<mutex> lock(writer_lock);
		return IntCount.find(I) != IntCount.end();
	}

	void AddInt(Node* I)
	{
		lock_guard<mutex> lock(writer_lock);
		++IntCount[I];
	}

	static string EID2Str(char c)
	{
		if (c == 0) return "ES";
		if (c == 1) return "IE";
		if (c == 2) return "AD";
		if (c == 3) return "AA";
		if (c == 4) return "AP";
		if (c == 5) return "XX";
		return "";
	}

	static char Str2EID(string s)
	{
		if (s == "ES") return 0;
		if (s == "IE") return 1;
		if (s == "AD") return 2;
		if (s == "AA") return 3;
		if (s == "AP") return 4;
		if (s == "XX") return 5;
		return -1;
	}

	private:
	
	mutex writer_lock;
};

class MT_Aligner
{
	public:

	MT_Aligner(BamParserOpt& opt) : opt(opt), gtf(opt.gtf), CNT(0)
	{
		if (opt.infa  != "")
		{
			cerr << "[quant] loading FASTA file ... " << flush;
			LoadFA(); 
			cerr << "done" << endl;
		}
		if (opt.gtf   != "") 
		{
			cerr << "[quant] loading GTF file ... " << flush;
			LoadGTF();
			cerr << "done" << endl;
		}
		LoadKallisto();
		//if (buffer.size()) FlushBuffer();
		cerr << "[quant] aligned " << D.num_mapped << " / " 
             << D.num_unmapped + D.num_mapped << " reads" << flush; 
		if (opt.alt != "")
			cerr << ", " << D.AltHeader.size() << " novel events" << endl;
		else
			cerr << endl;
		cerr << "[quant] writing outputs ... " << flush;
		if (opt.outfq != "") WriteFQ();
		if (opt.cnt   != "") WriteCnt();
		if (opt.alt   != "") WriteAlt();
		cerr << "done" << endl;
	}

	bool KallistoCallback(bam_hdr_t* H, const std::vector<bam1_t>& bv, intptr_t id)
	{
		/*if (threadId == thread::id())
		{
			lock_guard<mutex> lock(writer_lock);
			threadId = this_thread::get_id();	
		}*/
		vector<const bam1_t*> B;
		for (auto& b : bv) ProcessAlignment(H, &b, B);		
		/*if (this_thread::get_id() == threadId && buffer.size() > 8192 * 16) */
	    FlushBuffer(H, B);
		return true;
	}

	BamParserOpt& opt;

	private:

	//thread::id threadId;
	
	//bam_hdr_t* header;

	void FlushBuffer(bam_hdr_t* H, vector<const bam1_t*>& B)
	{
		lock_guard<mutex> lock(writer_lock);
		for (const bam1_t* b : B)
			int r = sam_write1(opt.bamfp, H, b);
	}

	mutex writer_lock;
	int CNT;

	map<string, string> fa;			
	AltData D;	
	GTF gtf;

	void LoadFA()
	{
		Tokenizer T(opt.infa);
		for (string s = "", t = T.getToken(); T.nextLine(); t = T.getToken())
			if (t[0] == '>')
				s = t.substr(1, t.find(" "));
			else if (s != "")
				fa[s] += t;
	}

	void LoadGTF()
	{
		if (opt.alt != "") gtf.genTranscripts();
	}

	void LoadKallisto()
	{
		vector<string> args =
		{
			"kallisto",
			"quant",
			"-t", to_string(opt.thr),
			"-l", to_string(opt.rl),
			"-s", "1",
			"-o", ".",
			"-i", opt.ind,
			"--single",
			"--pseudobam",
			opt.infq
		};
		int argc = args.size();
		char* argv[argc];

		for (unsigned int i = 0; i < argc; ++i)
			argv[i] = const_cast<char*>(args[i].c_str());

		(opt.cb)(argc, argv, (void*)this);
	}

	void WriteFQ()
	{
		if (GetFileExtension(opt.infq) == ".gz")
		{
			ZTokenizer FQ(opt.infq);
			if (GetFileExtension(opt.outfq) != ".gz")
				opt.outfq += ".gz";
			ZWriter out(opt.outfq);
			for (int k = 0; k < D.ReadCount.size(); )
		    {
		        if (!FQ.nextLine()) break;
		        string read = FQ.getToken().substr(1);
		        auto p = D.ReadCount.find(read);
		        bool b = p != D.ReadCount.end() && (p->second > opt.n || p->second == -1);
				string s = "";
		        if (b) s += FQ.getLine() + " " + to_string(p->second) + "\n";
		        if (!FQ.nextLine()) break;
		        if (b) s += FQ.getLine() + "\n";
		        if (!FQ.nextLine()) break;
		        if (b) s += FQ.getLine() + "\n";
		        if (!FQ.nextLine()) break;
		        if (b) s += FQ.getLine() + "\n";
				if (b) out.addStr(s);
				if (b) ++k;
		    }
			out.addStr("", 1);
		}
		else
		{
			Tokenizer FQ(opt.infq);
			ofstream out(opt.outfq);
		    for (int k = 0; k < D.ReadCount.size(); )
		    {
		        if (!FQ.nextLine()) break;
		        string read = FQ.getToken().substr(1);
		        auto p = D.ReadCount.find(read);
		        bool b = p != D.ReadCount.end() && (p->second > opt.n || p->second == -1);
		        if (b) out << FQ.getLine() << ' ' << p->second << "\n";
		        if (!FQ.nextLine()) break;
		        if (b) out << FQ.getLine() << "\n";
		        if (!FQ.nextLine()) break;
		        if (b) out << FQ.getLine() << "\n";
		        if (!FQ.nextLine()) break;
		        if (b) out << FQ.getLine() << "\n";
				if (k % 10000 == 0) out << flush;
				if (b) ++k;
		    }
			out << flush;
		}
	}

	void WriteCnt()
	{
		ofstream out(opt.cnt);

		if (opt.outfq == "")
		{
			for (auto& it : D.FragCount)
				out << it.first << '\t' << it.second.count << endl;
		}
		else
		{
			for (auto& it : D.FragCount)
			{
				out << it.first;
				for (int i = 0; i < it.second.positions.size(); ++i)
				{
					out << '\t' << it.second.positions[i];
				}
				out << '\n';
			}			
		}
	}	

	void WriteAlt()
	{
		map<string, int> M;
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
					string s = l->chr + "\t" + to_string(l->e) + "\t" + to_string(r->s) + "\t"
                             + t + "\t" + g + "\t" + AltData::EID2Str(k) + "\t";
                    M[s] += it.second;  
				}
				++k;
			}
		}

		ofstream out(opt.alt);
		for (auto& it : M)
		{
			out << it.first << "\t" << it.second << "\n";
		}
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
	}

	void ProcessAlignment(bam_hdr_t* H, const bam1_t* A, vector<const bam1_t*>& B)
	{
		string read = bam_get_qname(A);

		if (A->core.flag == 4) 
		{
			if (opt.outfq != "")
				D.AddRead(read, -1);
			else
			{
				lock_guard<mutex> lock(writer_lock);
				++D.num_unmapped;
			}
			return;
		}
		if (A->core.flag != 0 && A->core.flag != 16)
			return;

		string frag = H->target_name[A->core.tid];
		Transcript T = MT_Aligner::Frag2Trans(gtf, frag);
		int rl = A->core.l_qseq;

		if (!T.nodes.size()) 
		{
			if (opt.outfq != "")
				D.AddRead(read, -1);
			else
			{
				lock_guard<mutex> lock(writer_lock);
				++D.num_unmapped;
			}
			return;
		}

		int d = -1;
		
		if (opt.n)
		{
			string a = Hts2Seq(bam_get_seq(A), rl);
			string b = MT_Aligner::Trans2Seq(fa, T, opt.rl);
			if ((int)b.size() <= A->core.pos + opt.rl - 1)
				d = opt.rl + 1;
			else
				d = SeqDistance(a, b.substr(A->core.pos, opt.rl));
		}	
	
		if (d > opt.n)
		{
			if (opt.outfq != "")
				D.AddRead(read, -1);
			else
			{
				lock_guard<mutex> lock(writer_lock);
				++D.num_unmapped;
			}
			return;
		}

		if (opt.cnt != "" || opt.alt != "")
		{
			T = MT_Aligner::TrimTrans(T, A->core.pos, opt.rl, rl);
		}

		if (opt.cnt != "")
		{
			D.AddFrag(MT_Aligner::Trans2Frag(T), T.pos);
		}
		else
		{
			lock_guard<mutex> lock(writer_lock);
			++D.num_mapped;
		}
		
		if (opt.alt != "") 
			ProcessAlt(T);

		if (opt.bam != "")
		{
			B.push_back(A);
		}
	}

	public:

	static Transcript Frag2Trans(GTF& gtf, string frag)
	{
		vector<string> v = SplitString(frag, '|');
		Transcript T;
		for (int i = 0; i < v.size(); ++i)
		{
			Node* n = gtf.getNode(v[i]); 
			if (!n) return Transcript();
            T.insert(n);
		}
		T.sort();
		return T;
	}	

	static string Trans2Seq(map<string, string>& FA, Transcript& T, int rl)
	{
		string& c = FA[T.nodes.front()->chr];
		int sn = c.size(); string s = "";

		for (int i = 0; i < T.nodes.size(); ++i)
		{
			Node *n = T.nodes[i];	
			int st = n->s - 1;
			int en = n->e;
	
			if (rl && T.nodes.size() > 1 && !i && n->size() > rl)
				st = en - rl + 1;

			if (rl && T.nodes.size() > 1 && i == T.nodes.size() - 1 && n->size() > rl)
				en = st + rl - 1;

			for (int j = st; j < en && j < sn; ++j)
				s += c[j]; 
		}
		return s;
	}

	static Transcript TrimTrans(Transcript& T, int s, int rl, int len)
	{
		Transcript N;
		if (T.nodes.size() > 1 && T.nodes[0]->size() > rl - 1) 
			s += T.nodes[0]->size() - rl + 1;
		N.pos = s;
		for (int i = 0, c = 0; i < T.nodes.size(); ++i)
		{
			c += T.nodes[i]->size();
			if (s + 1 <= c) 
				N.nodes.push_back(T.nodes[i]); // s <= c
			else	
				N.pos -= T.nodes[i]->size();
			if (s + len <= c) break; // <=
		}
		return N;
	}

	static string Trans2Frag(Transcript& T)
	{
		string s = "";
		for (int i = 0; i < T.nodes.size(); ++i)
			s += T.nodes[i]->id + "|";
		return s;
	}

	static string GetGenes(Node* l, Node* r)
	{
		vector<string> v(SetIntersection(l->gene, r->gene));
		string g = "";		
		for (int i = 0; i < v.size(); ++i)
			g += v[i] + (i == v.size() - 1 ? "" : ",");
		return g;
	}

	private:

	static vector<string> SetIntersection(set<string>& A, set<string>& B)
	{
		vector<string> R;
		set_intersection(A.begin(), A.end(), B.begin(), B.end(), back_inserter(R));
		return R;
	}

	static string GetFileExtension(string& name)
	{
		size_t p = name.find_last_of(".");
		if (p == string::npos) return "";
		return name.substr(p + 1);
	}
};

void BamParse(BamParserOpt& opt);

#endif
