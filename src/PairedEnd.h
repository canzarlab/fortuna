//  ./fortuna --pairs -lf l.txt -rf r.txt -out o.txt -lp "/1" -rp "/2"

#ifndef PAIREDEND_H
#define PAIREDEND_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "BamParser.h"

using namespace std;

class PairedEndOpt
{
	public:

    string leftFile, rightFile;
	string leftPostfix, rightPostfix;
    string outFile;

	PairedEndOpt() : 
		leftFile(""), rightFile(""),
		leftPostfix(""), rightPostfix(""),
		outFile("")
	{ }
};

map<string, string> ParseFile(string file, string readPostfix)
{
	Tokenizer T(file);
	map<string, string> reads;

	int postfixSize = readPostfix.size();

	while (T.nextLine())
	{
		string fragmentId = T.getToken();

		for (string readId = T.getToken(); readId != ""; readId = T.getToken())
		{
			int start = readId.size() - postfixSize;

			if (start > 0)
			{
				if (readId.substr(start, postfixSize) == readPostfix)
				{
					readId = readId.substr(0, start);
				}
			}

			reads[readId] = fragmentId;
		}
	}

	return reads;
}

void PairedEndCount(PairedEndOpt& opt)
{
	map<string, string> lReads = ParseFile(opt.leftFile, opt.leftPostfix);
    map<string, string> rReads = ParseFile(opt.rightFile, opt.rightPostfix);
	map<string, int> pairs;

	for (auto& lRead : lReads)
	{
		string readId = lRead.first;
        string lFragId = lRead.second;

		if (rReads.find(readId) != rReads.end())
		{
			++pairs[lFragId + "^" + rReads[readId]];
		}
		else
		{
			++pairs[lFragId + "^"];
		}
	}

	for (auto& rRead : rReads)
	{
		string readId = rRead.first;
        string rFragId = rRead.second;

		if (lReads.find(readId) == lReads.end())
		{
			++pairs["^" + rFragId];
		}
	}

	ofstream out(opt.outFile);

	for (auto& pair : pairs)
	{
		out << pair.first << '\t' << pair.second << "\n";
	}
	out << flush;

	out.close();
}

#endif
