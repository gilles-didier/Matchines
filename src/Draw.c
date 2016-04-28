/*
    MATCHING MACHINES PACKAGE provides various softwares for studying the asymptotic behavior of pattern matching algorithm and for designing efficient algorithms
    Copyright (C) 2016  Gilles DIDIER and Laurent TICHIT

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/




#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "Utils.h"
#include "Text.h"
#include "Bernoulli.h"
#include "MatchingMachine.h"
#include "PatternMatchingAlgorithms.h"
#include "AsymptoticSpeed.h"
#include "Explore.h"

#define INC_BUFFER_CHAR 300
#define SIZE_BUFFER_CHAR 300
//#define ALPHABET "abcdefghijklmnopqrstuvwxyz"
#define ALPHABET "ab"
#define HELPMESSAGE "\n\nNAME\n\tdraw - return full lattice and matching machines of various pattern matching algorithms in dot '.gv' format\n\t\nSYNOPSIS\n\tdraw [OPTIONS] [PATTERN]\n\nDESCRIPTION\n\treturn full lattice and matching machines of various pattern matching algorithms corresponding to a given pattern in dot '.gv' format\n\t\n\t-a [ALPHABET]\n\t\tset the alphabet\n\t-b [FILE]\n\t\tset the file containing the Bernoulli model\n\t-n [ORDER]\n\t\tset the max order of the heuristic\n\t-h\n\t\tdisplay help\n\n"
#define RANDOM_SIZE 10000


//./draw -a ab -b bernoulli0.csv abaa

int main(int argc, char **argv) {
	char *inputFileBernoulli = NULL, output[SIZE_BUFFER_CHAR], *alphabet = NULL, *motif, option[256];
	TypeSymbol *pattern;
	TypeCardinal cardinal;
	int i, iarg, nbFunct = 10, nbNew = 3, length;
	TypeBernoulli *bernoulli;
	FILE * fb;
	
	TypeMatchingMachine *(*straFunct[])(TypeCardinal, TypeSymbol*, int) = {getNaiveMatchingMachine, getMorrisPrattMatchingMachine, getKnuthMorrisPrattMatchingMachine, getQuickSearchMatchingMachine, getHorspoolMatchingMachine, getFJSMatchingMachine, getTVSBMatchingMachine, getEBOMMatchingMachine, getHASHqMatchingMachine, getReverseMatchingMachine};
	char (*nameFunct[]) = {"Naive", "Morris-Pratt", "Knuth-Morris-Pratt", "Quicksearch", "Horspool", "FJSS", "TVSBS", "EBOM", "Hashq", "Reverse"};
	
	for(i=0; i<256; i++)
		option[i] = 0;
	for(iarg=1; iarg<argc && *(argv[iarg]) == '-'; iarg++) {
		int j;
		for(j=1; argv[iarg][j] != '\0'; j++)
			option[(unsigned char) argv[iarg][j]] = 1;
		if(option['a']) {
			option['a'] = 0;
			if((iarg+1)<argc)
				alphabet = argv[++iarg];
			else
				exitProg(ErrorArgument, "need an alphabet after option a");
		}
		if(option['b']) {
			option['b'] = 0;
			if((iarg+1)<argc)
				inputFileBernoulli = argv[++iarg];
			else
				exitProg(ErrorArgument, "need a file name after option b");
		}
		if(option['n']) {
			option['n'] = 0;
			if((iarg+1)<argc && sscanf(argv[iarg+1], "%d", &nbNew) == 1)
				iarg++;
			else
				exitProg(ErrorArgument, "need a length after -l");
		}
		if(option['h']) {
			printf("%s\n", HELPMESSAGE);
			exitProg(ExitOk, NULL);
		}
	}
	if(iarg<argc)
		motif = argv[iarg++];
	else {
		fprintf(stderr, "Missing pattern argument.\n");
		exit(1);
	}		
	if(alphabet == NULL)
		alphabet = getAlphabet(motif, NULL);
	length = strlen(motif);
	cardinal = strlen(alphabet);
	pattern = toSymbolSequence(motif, length, alphabet);
	if(inputFileBernoulli == NULL || !(fb = fopen(inputFileBernoulli, "r"))) {
		bernoulli = getUniformBernoulli(cardinal);
		fprintf(stderr, "No Bernoulli file provided or error. Use uniform probabilities.\n");;
	} else {
		bernoulli = freadBernoulli(fb, alphabet);
		fclose(fb);
	}
	{
		FILE *fo;
		TypeLattice *lattice;
		lattice = getLattice(cardinal, pattern, length);
		sprintf(output, "Lattice_%s.gv", motif);
		if((fo = fopen(output, "w"))) {
			fprintLatticeDot(fo, lattice, alphabet);
			fclose(fo);
		}
		sprintf(output, "Lattice_%s.tex", motif);
		if((fo = fopen(output, "w"))) {
			fprintLatticeGasTex(fo, lattice, alphabet);
			fclose(fo);
		}
		freeLattice(lattice);
	}
	for(i=0; i<nbFunct; i++) {
		FILE *fo;
		TypeMatchingMachine *stra, *expa;
		stra = (*straFunct[i])(cardinal, pattern, length);
		sprintf(output, "%s_%s_std.gv", nameFunct[i], motif);
		if((fo = fopen(output, "w"))) {
			fprintMatchingMachineDot(fo, stra, alphabet);
			fclose(fo);
		}
		expa = expandMatchingMachine(stra);
		sprintf(output, "%s_%s_exp.gv", nameFunct[i], motif);
		if((fo = fopen(output, "w"))) {
			fprintMatchingMachineDot(fo, expa, alphabet);
			fclose(fo);
		}
		freeMatchingMachine(stra);
		freeMatchingMachine(expa);
	}	
	for(i=1; i<=nbNew; i++) {
		FILE *fo;
		TypeMatchingMachine *stra;
//		stra = getNewKMatchingMachine(cardinal, pattern, length, bernoulli, i);
		stra = getPolyKMatchingMachine(cardinal, pattern, length, bernoulli, i);
		sprintf(output, "%d%s_%s_std.gv", i, HEURISTIC, motif);
		if((fo = fopen(output, "w"))) {
			fprintMatchingMachineDot(fo, stra, alphabet);
			fclose(fo);
		}
		freeMatchingMachine(stra);
	}
	if(strlen(motif)<=4) {
		FILE *fo;
		TypeMatchingMachine *stra;
		stra = explore(cardinal, pattern, length, bernoulli);
		sprintf(output, "%s_%s_std.gv", FASTEST, motif);
		if((fo = fopen(output, "w"))) {
			fprintMatchingMachineDot(fo, stra, alphabet);
			fclose(fo);
		}
		freeMatchingMachine(stra);
	}
	free((void*)pattern);
	freeBernoulli(bernoulli);
	return EXIT_SUCCESS;
}
