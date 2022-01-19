#ifndef KALLISTOINPUT_CPP
#define KALLISTOINPUT_CPP

#include "KallistoInput.h"

/*
 ./RNAseq --MT1 -ingff "../data/toy/1.gtf" -infst "../data/toy/1.fa" -tmp "." -outgff "../outs/1.gtf" -outfst "../outs/1.fa" -rl 25 -mil 0 -Mil 0 -mel 0 -Mel 0 -Mc 10000 -Moh 0 -lgo "22" --fse --pmrna

*/

void GffIn::parse()
{
	clock_t time = clock(); 
	
	if (!opt.silent) 
		std::cout << "[build] Parsing GFF... " << std::flush;
		
	std::ifstream in(opt.gffFile);	
	for(std::string line; std::getline(in, line);)
	{
		if (!line.size()) continue;

	  unsigned char k = 0;
		std::string r[9];
	  
	  for (unsigned int i = 0; i < line.length(); ++i)
			if (line[i] == '\t')
				++k;
	  	else
	  		r[k] += line[i];

	  if (r[2] != "subexon") continue;
	  
		std::string a[16], b[16];
    k = 0;
    
    for (unsigned int i = 0; i < r[8].length(); ++i)
    	if (r[8][i] == ' ' || i + 1 == r[8].length())
				++k;
	  	else if (r[8][i] != '\"' && r[8][i] != ';')
	  		b[k] += r[8][i];

		for (unsigned int i = 0; i < k; i += 2)
			if (b[i] == "NodeId")
				a[11] = b[i + 1];
			else if (b[i] == "gene_id")
				a[13] = b[i + 1];
      			else if (b[i] == "transcript_id")
				a[15] = b[i + 1];
			else if (b[i] == "SpliceEnd")
				a[7] = b[i + 1];

		if (!M.count(a[11])) 
		{
			GffIn::Block* block = new GffIn::Block();
			M[a[11]] = block;

			block->chromId = r[0];

			/*if (r[0] == "MT") // TMP block->chromId = r[0];
				block->chromId = "chrM";
			else
				block->chromId = "chr" + r[0];*/

			block->geneId   = a[13];
			block->sId      = std::stoi(r[3]); 
			block->eId      = std::stoi(r[4]);		
		  block->strandId = r[6];
		  block->nodeId   = a[11];
		  
		  if (r[6] == "-" && a[7] == "L")
		  	block->spliceEnd = "R";
		  else if (r[6] == "-" && a[7] == "R")
		  	block->spliceEnd = "L";
		  else
		  	block->spliceEnd = a[7];

		  G[a[13]].push_back(a[11]); 
     }
	 else
	 {
		if (r[6] == "-" && a[7] == "L")
		  	M[a[11]]->spliceEnd = sideSum(M[a[11]]->spliceEnd, "R");
		  else if (r[6] == "-" && a[7] == "R")
		  	M[a[11]]->spliceEnd = sideSum(M[a[11]]->spliceEnd, "L");
		  else
		  	M[a[11]]->spliceEnd = sideSum(M[a[11]]->spliceEnd, a[7]);		
	 }		
		if (a[13] == M[a[11]]->geneId && std::find(T[a[13]][a[15]].begin(), T[a[13]][a[15]].end(), a[11]) == T[a[13]][a[15]].end())
			T[a[13]][a[15]].push_back(a[11]);
		Z[a[11]].push_back(a[15]);
	}
	in.close();
	
	if (!opt.silent)
		std::cout << "creating data structures... " << std::flush;
	
	for (auto& it : Z)
		std::sort(it.second.begin(), it.second.end());

	for(auto& it : G)
	{
		auto il = it.second.begin();
		auto ir = it.second.end();
		std::sort(il, ir, [this] (const std::string& a, const std::string& b) 
		{ 
			GffIn::Block* l = M[a];
			GffIn::Block* r = M[b];
			return l->sId < r->sId; 
		});
	}
	
	for(auto& git : T)
	{
		for(auto& it : git.second)
		{
			std::vector<std::string>& V = it.second;
		
			auto il = V.begin();
			auto ir = V.end();
			std::sort(il, ir, [this] (const std::string& a, const std::string& b) 
			{ 
				GffIn::Block* l = M[a];
				GffIn::Block* r = M[b];
				return l->sId < r->sId; 
			});
			
			unsigned int e = V.size();
			
			if (opt.minExonLength || opt.maxExonLength)
			{
				unsigned int s = 0;
				for(unsigned int i = 0; i < e; ++i)
					if (i == e - 1 || M[V[i + 1]]->sId - M[V[i]]->eId > 1)
					{
						E[V[s]][V[i]] = it.first;
						s = i + 1;
					}
			}
			
			for(unsigned int i = 0; i < e; ++i)
			{
				GffIn::Block* block = M[V[i]];
				if (block->spliceType == "X")
					continue;
				else if (e == 1)
				{
					block->spliceType = GffIn::spliceSum(block->spliceType, "S"); 
					block->spliceType = GffIn::spliceSum(block->spliceType, "E"); 
				}
				else if (!i)
				{
					if (M[V[1]]->sId - block->eId > 1)
						block->spliceType = GffIn::spliceSum(block->spliceType, "S-");
					else
						block->spliceType = GffIn::spliceSum(block->spliceType, "S");
				}
				else if (i == e - 1)
				{
					if (block->sId - M[V[e - 2]]->eId > 1)
						block->spliceType = GffIn::spliceSum(block->spliceType, "E-");
					else
						block->spliceType = GffIn::spliceSum(block->spliceType, "E");
				}
				else
					block->spliceType = "X";
			}
					
			if (opt.samePreMRNA)
			{
				unsigned int l = M[V[0]]->sId;
				unsigned int r = M[V[e - 1]]->eId;
				P[git.first].push_back(std::make_pair(l, r));
			}
		}
	}
	
	if (opt.samePreMRNA)
		for (auto& it : P)
		{
			auto il = it.second.begin();
			auto ir = it.second.end();
			std::sort(il, ir, [this] (const UIntPair& a, const UIntPair& b) 
			{ 
				return a.first < b.first; 
			});
		}

	if (!opt.silent)
		std::cout << "parsed " << M.size() << " subexons from " << G.size() << " genes in " << (clock() - time) / (double)CLOCKS_PER_SEC << "s." << std::endl; 
}

std::string GffIn::sideSum(std::string l, std::string r)
{
	if (l == "-") return r;
	if (r == "-" || l == r) return l;
	return "B";
}

void GffIn::write()
{
	clock_t time = clock();
  unsigned int s = 0;
  
	if (!opt.silent)
		std::cout << "[build] Creating combinations... " << std::flush;

	std::map<std::string, std::map<std::string, std::vector<std::string>>> TT(T);

  std::ofstream out(opt.tmpPath + "tmp");
  for(auto& it : G)
	{
		unsigned int lgo_s = opt.largeGeneOpt[0] - '0';
		unsigned int lgo_e = opt.largeGeneOpt[1] - '0';
	
		GffIn::Data D(it.second);
		if (!lgo_s)
		{ 
			DFS(D, 0, 0, 0);
			if (lgo_e) lgo_s = GffIn::MAX_LGO_LEVEL; 
		}
		else
			D.combsNum = opt.maxCombos; 

		for (unsigned int k = lgo_s; lgo_e && k >= lgo_e; --k)
			if (D.combsNum >= opt.maxCombos)
			{
				if (k == 2) for(auto& i : TT[it.first])
				{
					
					unsigned int se = M[i.second[0]]->sId;
					unsigned int ee = M[i.second[i.second.size() - 1]]->sId;
	
					for (unsigned int j = 0; j < it.second.size() && M[it.second[j]]->sId < ee; ++j)
						if (M[it.second[j]]->sId < se) continue;
						else if (std::find(i.second.begin(), i.second.end(), it.second[j]) == i.second.end())
							i.second.push_back(it.second[j]);
					
					std::sort(i.second.begin(), i.second.end(), [this] (const std::string& a, const std::string& b) 
					{ 
						return M[a]->sId < M[b]->sId; 
					});
				}

				D.large(k);
				if (k == 2) for(auto i : TT[it.first])
				{
					D.transcr = i.first;
					D.V = i.second;
					DFS(D, 0, 0, 0);
					D.C.clear();
				}
				else for(auto i : T[it.first])
                                {
                                        D.transcr = i.first;
                                        D.V = i.second;
                                        DFS(D, 0, 0, 0);
                                        D.C.clear();
                                }
			}
		s += !opt.maxCombos || D.combsNum < opt.maxCombos;
		out << "G:" << it.first << std::endl << D.O;
	}
	out << "G:THE_END" << std::endl;
	out.close();
	
	if (!opt.silent)
		std::cout << "processed " << s * 100.0 / G.size() << "\% genes in " << (clock() - time) / (double)CLOCKS_PER_SEC << "s." << std::endl;
}

bool GffIn::DFS(GffIn::Data& D, unsigned int id, unsigned int head, unsigned int sum)
{
	if (sum + 2 >= opt.readLength) return 0;
	
	std::string   cs = "";
	GffIn::Block* cb = 0;
	
	if (head)
	{
		cs = D.V[id - 1];
		cb = M[cs];
	}
	
	bool fOne = 0;
	
	for(unsigned int i = id + 1; i <= D.V.size(); i++)
		if (!head)
		{
			D.C.push_back(D.V[i - 1]);
			DFS(D, i, i, 0);
			D.C.pop_back();
		}
		else do 
		{
			if (opt.maxCombos && D.combsNum == opt.maxCombos)
				break;
			
			if (D.largeOpt && D.largeOpt == 1 && i > id + 1) 
				break;
			
			std::string   ns = D.V[i - 1];
			GffIn::Block* nb = M[ns];

			if (opt.exonSkip) // limit the exon skipping
			{
				unsigned int l = 0;
				for(unsigned int k = id - 1; k < i - 1; ++k)
					if (M[D.V[k + 1]]->sId - M[D.V[k]]->eId > 1) ++l;
				if (l > opt.exonSkip) break;
			}

		 	unsigned int l = nb->sId - cb->eId - 1;
			if (l) 
		 	{ 
				if (opt.maxIntrLength && l > opt.maxIntrLength) break; 
				if (opt.minIntrLength && l < opt.minIntrLength) break;
			
				if (cb->spliceEnd     == "L") break;
	      if (cb->spliceEnd     == "-") break;
	      if (nb->spliceEnd     == "R") break;
	      if (nb->spliceEnd     == "-") break;
	        
	      if (cb->spliceType[0] == 'B') break;
	      if (cb->spliceType[0] == 'E') break;
	      if (nb->spliceType[0] == 'B') break;
	      if (nb->spliceType[0] == 'S') break;

	      if (!D.largeOpt && D.C.size() > 2 && E[cs][cs] != "") 
	      {
		    	if (opt.minExonLength && cb->length() < opt.minExonLength)
		    		break;
		    	if (opt.maxExonLength && cb->length() > opt.maxExonLength)
		    		break;
	      }
	    }
			else 
			{
				std::vector<std::string> CE;
				std::set_intersection(Z[cs].begin(), Z[cs].end(), Z[ns].begin(), Z[ns].end(), std::back_inserter(CE));
				if (!CE.size()) break;
			}

			D.C.push_back(ns);
			if (id != head)
				fOne = DFS(D, i, head, sum + cb->length()) || fOne;
			else
				fOne = DFS(D, i, head, 0) || fOne; // TMP
			D.C.pop_back();
		} while(0);

	if (!head) return 0;

	do
	{
		bool fSingle = !opt.forceSingleExons || id != head;
		if (!id) break;
		if (fSingle && opt.maxCombos && D.combsNum == opt.maxCombos) break; 
		if (sum + 2 > opt.readLength) break; 
		if (id == head && cb->length() < opt.readLength) break;
		if (id == head && D.S.count(D.C[0])) break; 
		if (id != head && fOne) break; // TMP
		if (id != head && sum + M[D.C[0]]->length() + cb->length() < opt.readLength) break; 

		if (!D.largeOpt && opt.maxOverhang && opt.readLength < opt.maxOverhang + sum && D.C.size() > 2) 
		{
			unsigned int l = D.C.size() - 1;
			if (M[D.C[1]]->sId - M[D.C[0]]->eId     > 1) break;
			if (M[D.C[l]]->sId - M[D.C[l - 1]]->eId > 1) break;
		}
		
		if(!D.largeOpt && opt.samePreMRNA)
		{
			bool fSame = 0;
			for(auto& it : P[cb->geneId])
			{
				fSame = it.first <= M[D.C[0]]->sId && cb->eId <= it.second;
				if (fSame) break;
			}
			if (!fSame) break;
		} 

		if (D.largeOpt == 2)
		{
			std::vector<std::string> CE;
			bool flag = 0;

			for (unsigned int i = 0; i < D.C.size(); ++i)
			{
				std::string fe = D.C[i]; CE.push_back(fe);
				std::string se = (i == D.C.size() - 1) ? "" : D.C[i + 1]; 
				if (se == "" || !(M[se]->sId - M[fe]->eId - 1))			
				{					
					for (unsigned int j = 0; !flag && j < CE.size(); ++j)
						flag = std::find(Z[CE[j]].begin(), Z[CE[j]].end(), D.transcr) != Z[CE[j]].end();											
					if (!flag) break;
					CE.clear(); 
					flag = (i == D.C.size() - 1) ? flag : 0;
				}
			}

			if (!flag) break;
		}

		if (id != head && !D.G.insert_supstr(D.C)) break; 
	
		if (id == head) 
			D.S[D.C[0]] = 1;
			
		if (fSingle) 
			++D.combsNum;
		
		for(auto& it : D.C)
			D.O += it + ' ';
		D.O += '\n';

		fOne = 1;
	} while(0);

	return fOne;
}

void GffIn::outFasta(std::string file) 
{
	clock_t time = 0;

	if (!opt.silent) 
	{
		std::cout << "[build] Reading FASTA input... " << std::flush;
		time = clock();
	}

	std::ifstream in(opt.fstFile);
	std::map<std::string, std::string> C;
	std::string chromId = "";
	
	for(std::string line; std::getline(in, line);)
	{
		if (!line.size()) continue;
		if (line[0] == '>')
			chromId = line.substr(1, line.find(' ', 0) - 1);
		else if (chromId != "")
			C[chromId] += line;
	}
  in.close();

	std::ofstream out(file);
	std::string geneId = "";

	if (!opt.silent) 
	{
		std::cout << "writing FASTA output... " << std::flush;
	}

	in.open(opt.tmpPath + "tmp");
	for(std::string line; std::getline(in, line);)
		if (!line.size()) 
			continue;
		else if (line.substr(0, 2) == "G:")
			geneId = line.substr(2, line.size());
		else if (geneId != "")	
		{
			unsigned int  num = 0;
			unsigned int  sum = 0;
			GffIn::Block* fst = 0;
			GffIn::Block* lst = 0;
			std::string s = "";
			std::string t = "";

			out << ">";	
			for(unsigned int i = 0; i < line.size(); ++i)
				if (line[i] == ' ')
				{
					GffIn::Block* b = M[s]; 
					t += C[b->chromId].substr(b->sId - 1, b->length());	
					out << b->nodeId << "|";
					if (!num)
						fst = b;
					else if (i == line.size() - 1)
						lst = b;
					else 
						sum += b->length();
					s = "";
					num++;
				}	
				else
					s += line[i];

			if (num >= 2 && t.size() >= opt.readLength + 2)
			{	
				unsigned int l = opt.readLength - 1;
				unsigned int s = 0;
				unsigned int e = t.size(); 
			
				if (fst->length() > l)
				{
					s  = fst->length() - l;
					e -= s;
				}
				
				if (lst->length() > l)
					e += l - lst->length();
			
				t = t.substr(s, e);
			}
		
			for (auto& c : t) c = toupper((unsigned char)c);	
			out << std::endl << t << std::endl;
		}
		
	in.close();
	out.close();
	
	if (!opt.silent)
		std::cout << "done in " << (clock() - time) / (double)CLOCKS_PER_SEC << "s." << std::endl;
}

void GffIn::outGff(std::string file) 
{
	clock_t time = 0;

	if (!opt.silent)
	{
		std::cout << "[build] Writing GFF output... " << std::flush;
		time = clock();
	}

	std::string geneId = "";
	std::ofstream out(file);
	std::ifstream in(opt.tmpPath + "tmp");
	unsigned int t = 0;

	for(std::string line; std::getline(in, line);)
		if (!line.size()) continue;
		else if (line.substr(0, 2) == "G:") 
			geneId = line.substr(2, line.size()); 
		else if (geneId != "") 
		{
			std::string s = "";
			++t;
			
			for(unsigned int i = 0; i < line.size(); ++i)
				if (line[i] == ' ')
				{
					Block* b = M[s];
					out << b->chromId << '\t';
					out << ".\texon\t";
					out << b->sId << '\t';
					out << b->eId << '\t';
					out << ".\t";
					out << b->strandId << '\t';
					out << ".\t";
					out << "gene_id \"" << b->geneId << "\"; ";
					out << "node_id \"" << b->nodeId << "\"; ";
					out << "transcript_id \"" << t << "\";";
					out << std::endl;
					s = "";
				}	
				else
					s += line[i];
		}
	
	in.close();
	out.close();
	
	if (!opt.silent)
		std::cout << "done in " << (clock() - time) / (double)CLOCKS_PER_SEC << "s." << std::endl;
}

std::string GffIn::spliceSum(std::string c, std::string n)
{
	if (c == "")  return n;
	if (n == "")  return c;
  if (c == n)   return n;
  
  if (c == "X" || n == "X")   return "X";	
	
	if (c == "S-" && n == "E-") return "X";
  if (c == "E-" && n == "S-") return "X";

	if (c == "S-" || c == "E-") return c;
	if (n == "S-" || n == "E-") return n;
	
	if (c == "S" && n == "E")   return "B";
  if (c == "E" && n == "S")   return "B";
	
	return c;
}

void GffIn::clear()
{
	for(auto& it : M) delete it.second;
	std::string name = opt.tmpPath + "tmp"; 
	std::remove(name.c_str());
	M.clear();
	G.clear();
	T.clear();
	E.clear();
	P.clear();
}

#endif
