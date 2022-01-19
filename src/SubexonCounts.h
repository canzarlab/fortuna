#ifndef SUBEXONCOUNTS_H
#define SUBEXONCOUNTS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <map>
#include <tuple>
#include <mutex>
#include <vector>

#include <htslib/sam.h>

#include "BamParser.h"

using namespace std;

class SubexonCountOpt
{
	public:

	string incnt, outcnt;

	SubexonCountOpt() : 
		incnt(""), outcnt("")
	{ }
};

void SubexonCount(SubexonCountOpt& opt)
{
	cerr << "[trans] loading CNT file ... " << flush;
	Tokenizer in(opt.incnt);
	map<string, int> M;
	while (in.nextLine())
	{
		string f = in.getToken(); 
		int    c = stoi(in.getToken());
		vector<string> v = SplitString(f, '|');
		for (int i = 0; i < v.size(); ++i)
			if (v[i].size()) M[v[i]] += c;
	}
	cerr << "done" << endl;

	cerr << "[trans] writing CNT file ... " << flush;
	ofstream out(opt.outcnt);
	for (auto& it : M)
		out << it.first << '\t' << it.second << endl;
	cerr << "done" << endl;
}

#endif
