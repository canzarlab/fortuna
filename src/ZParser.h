#ifndef ZPARSER_H
#define ZPARSER_H

#include <iostream>
#include <sstream>

//#include <bits/stdc++.h>
#include <stdio.h>
#include <zlib.h>

using namespace std;

#define CHUNK 256 * 1024

class ZTokenizer
{
	public:

	ZTokenizer(string infile) : 
		strm({0}), buff(""), line(""), token("")
	{
		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;
		strm.next_in = in;
		strm.avail_in = 0;
		inflateInit2(&strm, 15 | 32);
		file = fopen(infile.c_str(), "rb");
	}

	bool nextLine()
    {
        while(1)
		{
			if (feof(file) && !buff.size()) 
				break;

			size_t p = buff.find('\n');

			if (p == std::string::npos)
			{		
				strm.avail_in = fread(in, sizeof(char), sizeof(in), file);
				strm.next_in = in;
				do
				{
					strm.avail_out = CHUNK;
					strm.next_out  = out;
					if (!Valid(inflate(&strm, Z_NO_FLUSH)))
						return false;
					AddToBuffer();
				}
				while (!strm.avail_out);
			}
			else
			{
				line = buff.substr(0, p);
				buff = buff.substr(p + 1);
				if (!line.size()) continue;
				ss.clear();
				ss.str(line);
				return true;	
			}
		}

		return false;
    }

	string getToken()
	{
        return (ss >> token) ? token : "";
	}

	string getLine()
	{
		return line;
	}

	~ZTokenizer() 
	{
		inflateEnd(&strm);
		fclose(file); 
	}

	private:

	void AddToBuffer()
	{
		for (int i = 0; i < CHUNK && out[i] != '\0'; ++i)
			buff += out[i];
	}

	bool Valid(int e)
	{
		return (e == Z_OK || e == Z_STREAM_END); // || e == Z_BUF_ERROR);
	}

	string buff, line, token;
	istringstream ss;

	unsigned char in [CHUNK];
	unsigned char out[CHUNK];
	z_stream strm;
	FILE* file;
};

class ZWriter
{
	public:

	ZWriter(string infile) : strm({0})
	{ 
		strm.zalloc = Z_NULL;
    	strm.zfree  = Z_NULL;
    	strm.opaque = Z_NULL;
		file = fopen(infile.c_str(), "wb");
		deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 | 16, 8, Z_DEFAULT_STRATEGY);
	}

	bool addStr(string str, bool fin = 0) 
	{ 
		strm.avail_in = str.size();
	    strm.next_in  = (unsigned char*)str.data();

	    do 
		{
			strm.avail_out = CHUNK;
	        strm.next_out  = out;
			if (deflate(&strm, fin ? Z_FINISH : Z_NO_FLUSH) == Z_STREAM_ERROR)
				return false;
			int have = CHUNK - strm.avail_out;
	        if (fwrite(out, 1, have, file) != have || ferror(file)) 
	            return false;
		} while (strm.avail_out == 0);    

		return true;
	}

	~ZWriter() 
	{
		deflateEnd(&strm);
		fclose(file); 
	}

	private:

    unsigned char out[CHUNK];
	z_stream strm;
	FILE* file;
};

/*
void ZParser_Test(string infile)
{
	ZWriter W("2.txt.gz");
	for (int i = 0; i < 1000000; ++i)
	{
		string s = "";
        s += "@READ400000000.123456789.1.2.3.4.5.6.7.8.9 LUKA\n";
		s += "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n";
        s += "+\n";
        s += "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
		W.addStr(s);
	}

	ZTokenizer T("2.txt.gz");
	while (T.nextLine())
		cerr << T.getLine() << endl;	
}
*/

#endif
