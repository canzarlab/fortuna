#ifndef SEQUTILS_CPP
#define SEQUTILS_CPP

#include "SeqUtils.h"

char BaseComplement(char c)
{
	switch (c)
	{
		case 'A' : return 'T';
		case 'C' : return 'G';
		case 'G' : return 'C';
		case 'T' : return 'A';
		default : return 'N';
	};
}

/*int SeqDistance(string a, string b)
{
	int n = a.size(), c1 = 0, c2 = 0;
    for (int i = 0; i < n; ++i)
	{
        c1 += a[i] != b[i];
		c2 += a[i] != BaseComplement(b[n - i]);
	}
    return min(c1, c2);
}*/

int SeqDistance(string a, string b)
{
	int c = 0;
    for (int i = 0; i < a.size(); ++i)
        c += toupper(a[i]) != toupper(b[i]);
    return c;
}

string Hts2Seq(unsigned char* A, int n)
{
	string s(n, ' ');
	for (int i = 0; i < n; ++i)
		switch (bam_seqi(A, i))
		{
			case 1:  s[i] = 'A'; break;
			case 2:  s[i] = 'C'; break;
			case 4:  s[i] = 'G'; break;
			case 8:  s[i] = 'T'; break;
			default: s[i] = 'N';
		}
	return s;
}

vector<string> SplitString(string s, char c)
{
  	vector<string> V({""});
	for(unsigned int i = 0; i < s.size(); i++)
		if (s[i] != c)
			V.back() += s[i];
		else if (i < s.size() - 1)
			V.push_back("");			
	return V;
}

template<typename T> vector<T> SetIntersection(set<T>& A, set<T>& B)
{
	vector<T> R;
    set_intersection(A.begin(), A.end(), B.begin(), B.end(), back_inserter(R));
	return R;
}

template<typename T> vector<T> SetUnion(set<T>& A, set<T>& B)
{
	vector<T> R;
    set_union(A.begin(), A.end(), B.begin(), B.end(), back_inserter(R));
	return R;
}

#endif
