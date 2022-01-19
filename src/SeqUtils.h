#ifndef SEQUTILS_H
#define SEQUTILS_H

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <htslib/sam.h>

using namespace std;

char BaseComplement(char c);

int SeqDistance(string a, string b);

string Hts2Seq(unsigned char* A, int n);

vector<string> SplitString(string s, char c);

template<class T> vector<T> SetIntersection(set<T>& A, set<T>& B);

template<class T> vector<T> SetUnion(set<T>& A, set<T>& B);

#endif
