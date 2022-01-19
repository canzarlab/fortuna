#ifndef MAIN_CPP
#define MAIN_CPP

#include "config.h"
#include "ArgParser.cpp"
#include "KallistoInput.h"
#include "Graph.h"
#include "BamParser.h"
#include "CountRefine.h"
#include "SubexonCounts.h"

#include <kallisto/src/main.cpp>

#include <iostream>
#include <iomanip>

#include <stdio.h>
#include <zlib.h>

const int v_hi = 1;
const int v_mi = 1;
const int v_lo = 0;

const string tool = "fortuna";

string GetVersion()
{
	return to_string(v_hi) + "." + to_string(v_mi) + "." + to_string(v_lo);
}

int main(int argc, char* const argv[])
{
	ArgParser A(argc, argv);
	
	if (argc == 1)
	{
		cerr << "\n\t" << tool << "(R) " << GetVersion() << " tools:\n\n"
             << "\t --index  " << "\tgenerates a reusable alignment index\n"
		     << "\t --quant  " << "\tquantification step powered by kallisto\n"
		     << "\t --refine " << "\trefines the coutns using realigned unmapped reads\n"
             << "\t --trans "  << "\ttransforms a count file into subexon count file\n"
		     << "\t --version" << "\toutputs the current tool version\n"
		     << "\t --cite   " << "\toutputs how should the tool be cited\n"
             << "\n\tRun the tool without additional parameters to display help.\n"
             << endl;
	}
	else if (A.exists("index"))
	{
		if (argc == 2)
		{
			cerr << "\n\tindex options:\n"
				 << "\n\t -ingtf  <FILE>\t" << "input refined and grouped GTF file (*)"
				 << "\n\t -infa   <FILE>\t" << "input chromosome FASTA file (*)"
				 << "\n\t -outfa  <FILE>\t" << "outputs intermediate index fasta (*)"
				 << "\n\t -ind    <FILE>\t" << "outputs the final index file for quantification (*)"
			     << "\n\t -rl     <INT> \t"  << "target read length (*)"
                 << "\n\t -tmp    <FOLD>\t" << "temporary folder (defaults to \"./\")"
                 << "\n\t -lgo    <STR> \t"  << "large gene optimization strategy (default \"11\")"
                 << "\n\t -Mc     <INT> \t"  << "maximum number of fragments per gene (default 25000)"
                 << "\n\t -mil    <INT> \t"  << "minimum generated intron length (default 0)"
                 << "\n\t -Mil    <INT> \t"  << "maximum generated intron length (default 0)"
                 << "\n\t -mel    <INT> \t"  << "minimum generated exon length (default 0)"
                 << "\n\t -Mel    <INT> \t"  << "maximum generated exon length (default 0)" 
                 << "\n\t -exs    <INT> \t"  << "maximum number of skipped exons (default 0)"
                 << "\n\t -Moh    <INT> \t"  << "maximum overhang on fragment ends (default 0)"
                 << "\n\t --pmrna       \t"  << "generate frgaments from known premrna (default off)"
                 << "\n\t --fse         \t"  << "generate all single subexons (default off)"
                 << "\n\t -outgtf <FILE>\t"  << "outputs a gtf file describing fragments"
				 << "\n\n\tOption -lgo takes a string of up to two characters from the set "
                 << "\n\t{0, 2, 1} where 0 is the least restrictive option, 2 is the middle "
                 << "\n\tground and 1 is the most restrictive option. If the parameter "
                 << "\n\tstring consists of 2 characters, index will use the first strategy "
                 << "\n\tand fall back gradually to the second one if -Mc is exceeded. \n"
                 << endl;

			return 0;
		}


		GffIn::Options opt; 

		opt.silent = A.exists("silent");
	
		if (A.exists("ingtf"))
			opt.gffFile = A("ingtf");
		else
		{
			std::cerr << "[error] GTF file not specified." << std::endl;
			return 0;
		}		

		if (A.exists("infa"))
			opt.fstFile = A("infa");
		else
		{
			std::cerr << "[error] Input FASTA file not specified." << std::endl;
			return 0;
		}		

		if (A.exists("tmp"))
			opt.tmpPath = A("tmp");
		
		if (A.exists("rl"))
			opt.readLength = std::atoi(A("rl").c_str());
		else 
		{
			std::cerr << "[error] Read length not specified." << std::endl;
			return 0;
		}
		
		if (!opt.silent)	
			std::cout << "[param] readLength: " << opt.readLength << std::endl;
	
		if (A.exists("mil"))
			opt.minIntrLength = std::atoi(A("mil").c_str());
	
		if (!opt.silent)	
			std::cout << "[param] minIntrLength: " << opt.minIntrLength << std::endl;
	
		if (A.exists("Mil"))
			opt.maxIntrLength = std::atoi(A("Mil").c_str());
	
		if (!opt.silent)	
			std::cout << "[param] maxIntrLength: " << opt.maxIntrLength << std::endl;
		
		if (A.exists("mel"))
			opt.minExonLength = std::atoi(A("mel").c_str());
			
		if (!opt.silent)
			std::cout << "[param] minExonLength: " << opt.minExonLength << std::endl;
		
		if (A.exists("Mel"))
			opt.maxExonLength = std::atoi(A("Mel").c_str());

		if (!opt.silent)
			std::cout << "[param] maxExonLength: " << opt.maxExonLength << std::endl;

		if (A.exists("exs"))
    		opt.exonSkip = std::atoi(A("exs").c_str());

		if (!opt.silent)
			std::cout << "[param] exonSkip: " << opt.exonSkip << std::endl;

		if (A.exists("Mc"))
			opt.maxCombos = std::atoi(A("Mc").c_str());
	
		if (!opt.silent)
			std::cout << "[param] maxCombos: " << opt.maxCombos << std::endl;

		if (opt.maxCombos && A.exists("lgo"))
			opt.largeGeneOpt = A("lgo").c_str();
	
		if (opt.largeGeneOpt.size() == 1)
			opt.largeGeneOpt += opt.largeGeneOpt;
	
		if (!opt.silent)	
			std::cout << "[param] largeGeneOpt: " << opt.largeGeneOpt << std::endl;

		if (A.exists("Moh"))
			opt.maxOverhang = std::atoi(A("Moh").c_str());
	
		if (!opt.silent)	
			std::cout << "[param] maxOverhang: " << opt.maxOverhang << std::endl;

		opt.forceSingleExons = A.exists("fse");

		if (!opt.silent)
			std::cout << "[param] forceSingleExons: " << opt.forceSingleExons << std::endl;

		opt.samePreMRNA = A.exists("pmrna");

		if (!opt.silent)
			std::cout << "[param] samePreMRNA: " << opt.samePreMRNA << std::endl;
		
		GffIn in(opt); 
		in.parse(); 
		in.write();
		
		if(A.exists("outfa"))
			in.outFasta(A("outfa")); 
		else
		{
			std::cerr << "[error] Output FASTA file not specified." << std::endl;
			return 0;
		} 
	
		if(A.exists("outgtf"))
			in.outGff(A("outgtf"));

		string ind;
		if(A.exists("ind"))
			ind = A("ind");
		else
		{
			std::cerr << "[error] Output kallisto index not specified." << std::endl;
			return 0;
		}

		vector<string> args =
		{
			"kallisto",
			"index",
			"-i", ind,
			A("outfa")
		};
		int argc = args.size();
		char* argv[argc];

		for (unsigned int i = 0; i < argc; ++i)
			argv[i] = const_cast<char*>(args[i].c_str());
	
		main_kallisto(argc, argv, 0);

		if (!opt.silent) std::cout << std::endl;
	}
	else if (A.exists("quant"))
	{
		if (argc == 2)
		{
			cerr << "\n\tquant options:\n"
				 << "\n\t -rl     <INT> \t" << "target read length (*)"
                 << "\n\t -gtf    <FILE>\t" << "input refined and grouped GTF file (*)"
                 << "\n\t -infq   <FILE>\t" << "input FASTQ file with reads (*)" 
                 << "\n\t -ind    <FILE>\t" << "input index file (*)"
                 << "\n\t -fa     <FILE>\t" << "input transcriptome FASTA file"
				 << "\n\t -outfq  <FILE>\t" << "output FASTQ file with misaligned reads"
				 << "\n\t -cnt    <FILE>\t" << "output fragment count file"
				 << "\n\t -alt    <FILE>\t" << "output alternative splicing file"
				 << "\n\t -bam    <FILE>\t" << "output pseudoalignment BAM file"
			     << "\n\t -miss   <INT> \t" << "maximum number of per alignment mismatches"
                 << "\n\t -thr    <INT> \t" << "number of worker threads"
                 << "\n\n\tAlignment is done by kallisto. A read is considered to be misaligned "
                 << "\n\tif its unaligned by kallisto or if it has more than \"-miss\" mismatches.\n"
                 << endl;
			return 0;
		}		

		BamParserOpt opt;

		if (!A.exists("outfq") && !A.exists("cnt") && !A.exists("alt") && !A.exists("bam"))
		{
			std::cerr << "[error] No outputs specified." << std::endl;
			return 0;
		}

		if (A.exists("outfq"))
			opt.outfq = A("outfq");

		if (A.exists("infq"))
			opt.infq = A("infq");
		else
		{
			std::cerr << "[error] Input FASTQ file not specified." << std::endl;
			return 0;
		}

		if (A.exists("cnt"))
			opt.cnt = A("cnt");

		if (A.exists("gtf"))
			opt.gtf = A("gtf");
		else
		{
			std::cerr << "[error] GTF file not specified." << std::endl;
			return 0;
		}

		if (A.exists("alt"))
		{
			opt.alt = A("alt");

			if (A.exists("fa"))
				opt.infa = A("fa");
			else
			{
				std::cerr << "[error] FASTA file not specified." << std::endl;
				return 0;
			}

			if (A.exists("miss"))
				opt.n = stoi(A("miss"));
			else
			{
				std::cerr << "[error] Maximum number of mismatches not specified." << std::endl;
				return 0;
			}
		}

		if (A.exists("miss"))
		{
			opt.n = std::stoi(A("miss"));

			if (A.exists("fa"))
				opt.infa = A("fa");
			else
			{
				std::cerr << "[error] FASTA file not specified." << std::endl;
				return 0;
			}

			if (A.exists("gtf"))
				opt.gtf = A("gtf");
			else
			{
				std::cerr << "[error] GTF file not specified." << std::endl;
				return 0;
			}
		}

		if (A.exists("rl"))
			opt.rl = std::stoi(A("rl"));
		else
		{
			std::cerr << "[error] Read length not specified." << std::endl;
			return 0;
		} 

		if (A.exists("ind"))
			opt.ind = A("ind");
		else
		{
			std::cerr << "[error] Index file not specified." << std::endl;
			return 0;
		}

		if (A.exists("thr"))
			opt.thr = stoi(A("thr"));

		if (A.exists("bam"))
			opt.bam = A("bam");

		opt.cb = &main_kallisto;

		BamParse(opt);
	}
	else if (A.exists("refine"))
	{
		if (argc == 2)
		{
			cerr << "\n\trefine options:\n"
				 << "\n\t -rl     <INT> \t" << "target read length (*)"
                 << "\n\t -gtf    <FILE>\t" << "input refined and grouped GTF file (*)"
				 << "\n\t -bam    <FILE>\t" << "input genome BAM file with realigned reads (*)"
				 << "\n\t -outcnt <FILE>\t" << "output fragment count file"
				 << "\n\t -outalt <FILE>\t" << "output alternative splicing file"
                 << "\n\t -ref    <FILE>\t" << "output refinement reference file"
				 << "\n\t -incnt  <FILE>\t" << "input fragment count file (default none)"
				 << "\n\t -inalt  <FILE>\t" << "input alternative splicing file (default none)"
                 << "\n\n\tOptional inputs specified by \"-incnt\" and \"-inalt\" are going "
                 << "\n\tto be augmented with the new information from the BAM file and stored "
                 << "\n\tto \"-outcnt\" and \"-outalt\" respectively. Reference file contains "
                 << "\n\tthe information on newly created/refined subexons and introns.\n"
                 << endl;
			return 0;
		}
		
		CountRefineOpt opt;

		if (A.exists("bam"))
			opt.bam = A("bam");
		else
		{
			std::cerr << "[error] BAM file not specified." << std::endl;
			return 0;
		}
		
		if (A.exists("gtf"))
			opt.gtf = A("gtf");
		else
		{
			std::cerr << "[error] GTF file not specified." << std::endl;
			return 0;
		}

		if (A.exists("rl"))
			opt.rl = std::stoi(A("rl"));
		else
		{
			std::cerr << "[error] Read Length not specified." << std::endl;
			return 0;
		}

		if (!A.exists("outcnt") && !A.exists("outalt"))
		{
			std::cerr << "[error] Outputs not specified." << std::endl;
			return 0;
		}

		if (A.exists("outcnt"))
			opt.outcnt = A("outcnt");

		if (A.exists("outalt"))
			opt.outalt = A("outalt");

		if (A.exists("ref"))
			opt.ref = A("ref");

		if (A.exists("incnt"))
			opt.incnt = A("incnt");
	
		if (A.exists("inalt"))
			opt.inalt = A("inalt");

		CountRefine(opt);
	}	
	else if (A.exists("trans"))
	{
		if (argc == 2)
		{
			cerr << "\n\ttrans options:\n"
				 << "\n\t -incnt  <FILE>\t" << "input chromosome FASTA file (*)"
				 << "\n\t -outcnt <FILE>\t" << "outputs intermediate index fasta (*)"
				 << "\n\n\tTransforms regular count file into subexon counts."
                 << endl;
			return 0;
		}

		SubexonCountOpt opt;

		if (A.exists("incnt"))
			opt.incnt = A("incnt");
		else
		{
			std::cerr << "[error] input CNT file not specified." << std::endl;
			return 0;
		}

		if (A.exists("outcnt"))
			opt.outcnt = A("outcnt");
		else
		{
			std::cerr << "[error] output CNT file not specified." << std::endl;
			return 0;
		}

		SubexonCount(opt);
	}
	else if (A.exists("version"))
	{
		cerr << "\n\t" << tool << "(R) " << GetVersion() << "\n\n";
	}
	else if (A.exists("cite"))
	{
		cerr << "\n\t" << tool << "(R), Smart fridge crew (C)2020\n\n";
	}
	else if (A.exists("test"))
	{
		/*string infile = "";

		if (A.exists("infile"))
			infile = A("infile");
		else
		{
			std::cerr << "[error] input file not specified." << std::endl;
			return 0;
		}		

		ZParser_Test(infile);*/
	}

	return 0;
}

#endif
