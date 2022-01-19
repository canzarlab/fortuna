#ifndef BAMPARSER_CPP
#define BAMPARSER_CPP

#include "BamParser.h"




/*
vector<Node*> GetFragment(GTF& gtf, string f, int rl)
{
	vector<string> v = split(f, '|');
	vector<Node*> o;

	for (unsigned int i = 0; i < v.size(); ++i)
		o.push_back(gtf.getNode(v[i]));

	return o;
}

vector<Node*> TrimFragment(vector<Node*>& v, int s, int rl)
{
	// ADD CHROMOSOME AND TRANSCRIPTS/GENES TO NODES
	return vector<Node*>();
}

string GetSequence(map<string, string>& fa, vector<Node*>& v, string chr, int rl)
{
	string f = fa[chr];

	if (v.size() == 1)
		return f.substr(v[0]->s - 1, v[0]->size());
	
	string s;

	if (v[0]->size() < rl)
		s = f.substr(v[0]->s - 1, v[0]->size());	
	else
		s = f.substr(v[0]->e - rl + 1, rl - 1);	

	for (unsigned int i = 1; i < v.size() - 1; ++i)
		s += f.substr(v[i]->e - rl + 1, rl - 1);	

	if (v.back()->size() < rl)
		s += f.substr(v.back()->s - 1, v.back()->size());	
	else
		s += f.substr(v.back()->s - 1, rl - 1);

	return s;
}

string Frag2Str(vector<Node*>& v)
{
	string s = "";
	for (unsigned int i = 0; i < v.size(); ++i)
		s += v[i]->id + "|";
	return s;
}

void AS_ParseFragment(GTF& gtf, vector<Node*>& v, ofstream& alt)
{
}
*/

void BamParse(BamParserOpt& opt)
{
	MT_Aligner A(opt);

/*
    map<string, string> fa;
	if (opt.infa != "")
	{
		Tokenizer T(opt.infa);

		for (string s = ""; T.nextLine();)
			if (T.getToken()[0] == '>')
				s = T.getToken(), s = s.substr(1, s.find(" "));
			else if (s != "")
				fa[s] += T.getToken();				
	}

	GTF gtf(opt.gtf); 
	gtf.genTranscripts();

	samFile*   IN = hts_open(opt.bam.c_str(), "r");
	bam_hdr_t* HD = sam_hdr_read(IN);
	bam1_t*    A  = bam_init1();

	map<string, int> V;
	map<string, int> C;

	ofstream alt;
	if (opt.alt != "") alt.open(opt.alt);

	while(sam_read1(IN, HD, A) > 0)
		if (A->core.flag & 4)
    		V[bam_get_qname(A)] = -1; 
		else if (A->core.flag == 0 || A->core.flag == 16) 
        {
			vector<Node*> v = GetFragment(gtf, HD->target_name[A->core.tid], opt.rl);
			int d = -1;

			if (opt.n)
			{
				string a = conv(bam_get_seq(A), opt.rl);
	        	string b = GetSequence(fa, v, gtf.getChr(v[0]), opt.rl);

				if ((int)b.size() < A->core.pos + opt.rl)
					d = opt.n + 1;
				else
					d = dist(a, b.substr(A->core.pos, opt.rl));
			}

			v = TrimFragment(v, A->core.pos, opt.rl);
			string f = Frag2Str(v);

			if (opt.infq != "" && d > opt.n)
		        V[bam_get_qname(A)] = d;              

			if (opt.cnt != "" && d > opt.n)
				C[f] += 1;

			if (opt.alt != "")
				AS_ParseFragment(gtf, v, alt);
        }

	alt.close();
	bam_destroy1(A);
	sam_close(IN);

    if (V.size())
    {
		Tokenizer FQ(opt.infq);
        ofstream out(opt.outfq);
        int c = 0;

        while (1)
        {
            if (!FQ.nextLine()) break;
            string read = FQ.getToken().substr(1);
            auto p = V.find(read);
            bool b = p != V.end() && (p->second > opt.n || p->second == -1);
            if (b) out << FQ.getLine() << ' ' << p->second << "\n";
            if (!FQ.nextLine()) break;
            if (b) out << FQ.getLine() << "\n";
            if (!FQ.nextLine()) break;
            if (b) out << FQ.getLine() << "\n";
            if (!FQ.nextLine()) break;
            if (b) out << FQ.getLine() << endl;
            if (b) ++c;
        }
        cout << c << " reads." << endl;
    }

	if (C.size())
	{
		ofstream out(opt.cnt);
		int c = 0;

		for (auto it : C)
		{
			out << it.first << '\t' << it.second << endl;
			c += it.second;
		}

		cout << c << " counts." << endl;
	}*/
}


#endif
