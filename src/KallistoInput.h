#ifndef KALLISTOINPUT_H
#define KALLISTOINPUT_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <ctime>

#include "Graph.h"

typedef std::pair<unsigned int, unsigned int> UIntPair;

class GffIn
{
  public:

	class Options
  {
		public:
		std::string  gffFile;
		std::string  fstFile;
		std::string  tmpPath;
		unsigned int readLength;
		unsigned int minIntrLength;
		unsigned int maxIntrLength;
		unsigned int minExonLength;
		unsigned int maxExonLength;
		unsigned int maxCombos;
		unsigned int maxOverhang;		
		unsigned int exonSkip;
		std::string  largeGeneOpt;
		bool         forceSingleExons;
		bool         samePreMRNA;
		bool         silent;
		
		Options()
		{
			gffFile          = "";
			fstFile          = "";
		  tmpPath          = "";
			readLength       = 75;
		  minIntrLength    = 0;
		  maxIntrLength    = 0;
			minExonLength    = 0;
			maxExonLength    = 0;
			maxCombos        = 0;
			maxOverhang      = 0;
			exonSkip         = 0;		
			largeGeneOpt     = "00"; 
			forceSingleExons = 0;
			samePreMRNA      = 0;
			silent           = 0;
		}
	};

	GffIn::Options opt;

	GffIn(GffIn::Options& opt) : opt(opt) { }
	~GffIn() { clear(); }

	void parse();
	void write();
	void clear(); 	
	
	void outFasta(std::string file);
	void outGff(std::string file);

	private:

	static const unsigned int MAX_LGO_LEVEL = 2;

	class Block
	{
		public:

		std::string geneId;
		std::string nodeId;
		std::string chromId;
		std::string strandId;
		std::string spliceEnd;
		std::string spliceType;
		
		unsigned int sId;
		unsigned int eId;

		unsigned int length() { return eId - sId + 1; }
	} NULL_BLOCK;
	
	std::map<std::string, GffIn::Block*>                                   M; // Actual subexons

	std::map<std::string, std::vector<std::string>>                        G; // Genes
	std::map<std::string, std::map<std::string, std::vector<std::string>>> T; // Genes-Transcripts
	std::map<std::string, std::map<std::string, std::string>>              E; // Exons
	std::map<std::string, std::vector<UIntPair>>                           P; // PreMRNA
	std::map<std::string, std::vector<std::string>>                        Z; // Subexons-Transcripts
	
	class Data
	{
		public:

		std::string    					    O; // Output
		std::vector<std::string>    V; // Input
		std::vector<std::string>    C; // Current combination
		std::map<std::string, bool> S; // Singletons
		Graph                       G; // Combinations graph
		
		unsigned int combsNum;
		unsigned int largeOpt;
		std::string  transcr;  // Current transcript

		Data(std::vector<std::string>& V) : V(V)
		{
			combsNum = 0;
			largeOpt = 0;
		}
		
		void large(unsigned int optLevel)
		{
			largeOpt = optLevel;
			combsNum = 0;
			V.clear();
			C.clear();
			G.clear();
			S.clear();
			O = "";
		}
	};
	
	bool DFS(GffIn::Data& D, unsigned int id, unsigned int head, unsigned int sum);
	
	static std::string spliceSum(std::string a, std::string b);
	static std::string sideSum(std::string l, std::string r);
};

#endif
