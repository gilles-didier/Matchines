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
#include <math.h>

#include "Utils.h"
#include "MatchingMachine.h"
#include "PatternMatchingAlgorithms.h"
#include "AsymptoticSpeed.h"
#include "Explore.h"
#include "Text.h"
#include "TableUtils.h"


//./empirical -a acgt -p acccggaat ../data/genome/ecoli.txt 

#define INC_BUFFER_CHAR 300

#define SIZE_BUFFER_CHAR 300
//#define ALPHABET "abcdefghijklmnopqrstuvwxyz"
#define ALPHABET "ab"
#define HELPMESSAGE "\n\nNAME\n\ttheoretical - compute the asymptotic speed of various pattern matching algorithms under a Bernoulli model\n\t\nSYNOPSIS\n\ttheoretical [OPTIONS]\n\nDESCRIPTION\n\treturn a table containing the asymptotic speeds under a Bernoulli model of a given pattern or of all the patterns of a given size.\n\toptions -a and -b have to be filled\n\t\n\t-a [ALPHABET]\n\t\tset the alphabet\n\t-p [PATTERN]\n\t\tset the pattern\n\t-l [LENGTH]\n\t\tset the length of the patterns to be tested (not used if a pattern is given via option -p)\n\t-s [NUMBER]\n\t\tset the number of pattern samples to be tested (not used if a pattern is given via option -p)\n\t\tIf this option is not used, all the pattern of the selected length are tested.\n\t-x [FILE]\n\t\tset the file name of the output (default is 'table_result.xxx')\n\t-n [ORDER]\n\t\tset the max order of the heuristic\n\t-u\n\t\tevaluate heuristic and the fastest strategy computed from the uniform Bernoulli model\n\t-f [format]\n\t\tset the format of the output (l: LaTex tabular format, c/default: CSV format)\n\t-r [TYPE]\n\t\ttype of random used for sampling patterns (active with option 's'):\n\t\t\tm\tuse Bernoulli model from option b\n\t\t\tu\tuse uniform model\n\t-b [FILE]\n\t\tset the file containing the Bernoulli model\n\t-h\n\t\tdisplay help\n\n"
#define RANDOM_SIZE 10000
#define OUTPUT_DEFAULT "table_result"



void computeLine(double *val, TypeMatchingMachine *(*straFunct[])(TypeCardinal, TypeSymbol*, int), int nbFunct, int nbNew, TypeBernoulli *bernoulli, TypeBernoulli *bernoulli2, TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	int i, ind = 0;
	for(i=0; i<nbFunct; i++) {
		TypeMatchingMachine *stra, *expa;
		stra = (*straFunct[i])(cardinal, pattern, length);
		expa = expandMatchingMachine(stra);
		val[ind++] = getAsymptoticSpeed(expa, bernoulli);
		freeMatchingMachine(stra);
		freeMatchingMachine(expa);
	}
	if(bernoulli2 != NULL) {
		for(i=1; i<=nbNew; i++) {
			TypeMatchingMachine *stra;
	//		stra = getNewKMatchingMachine(cardinal, pattern, length, bernoulli2, i);
			stra = getPolyKMatchingMachine(cardinal, pattern, length, bernoulli2, i);
			val[ind++] = getAsymptoticSpeed(stra, bernoulli);
			freeMatchingMachine(stra);
		}
		if(length<=4) {
			TypeMatchingMachine *stra;
			stra = explore(cardinal, pattern, length, bernoulli2);
			val[ind++] = getAsymptoticSpeed(stra, bernoulli);
			freeMatchingMachine(stra);
		}
	}
	for(i=1; i<=nbNew; i++) {
		TypeMatchingMachine *stra;
//		stra = getNewKMatchingMachine(cardinal, pattern, length, bernoulli, i);
		stra = getPolyKMatchingMachine(cardinal, pattern, length, bernoulli, i);
		val[ind++] = getAsymptoticSpeed(stra, bernoulli);
		freeMatchingMachine(stra);
	}
	if(length<=4) {
		TypeMatchingMachine *stra;
		stra = explore(cardinal, pattern, length, bernoulli);
		val[ind++] = getAsymptoticSpeed(stra, bernoulli);
		freeMatchingMachine(stra);
	}
}

int main(int argc, char **argv) {
	char *inputFileBernoulli = NULL, *output = NULL, *alphabet = NULL, *motif = NULL, option[256], format='c', rtype = 'u', TMP[300];
	int i, iarg, length = 3, nbFunct = 9, nbNew = 3, p, unif = 0, nbSample = 0, sgn = 2;
	TypeBernoulli *bernoulli = NULL, *bernoulli2 = NULL, *bernoulliUnif;
	FILE *fo;
	TypeSymbol *pattern, cardinal;
	TypeMatchingMachine *(*straFunct[])(TypeCardinal, TypeSymbol*, int) = {getNaiveMatchingMachine, getMorrisPrattMatchingMachine, getKnuthMorrisPrattMatchingMachine, getQuickSearchMatchingMachine, getHorspoolMatchingMachine, getFJSMatchingMachine, getTVSBMatchingMachine, getEBOMMatchingMachine, getHASHqMatchingMachine, getReverseMatchingMachine};
	char (*nameFunct[]) = {"Naive", "Morris-Pratt", "Knuth-Morris-Pratt", "Quicksearch", "Horspool", "FJS", "TVSBS", "EBOM", "Hashq", "Reverse"};
	
	for(i=0; i<256; i++)
		option[i] = 0;
	for(iarg=1; iarg<argc && *(argv[iarg]) == '-'; iarg++) {
		int j;
		for(j=1; argv[iarg][j] != '\0'; j++)
			option[(unsigned char) argv[iarg][j]] = 1;
		if(option['a']) {
			option['a'] = 0;
			if((iarg+1)<argc) {
				alphabet = argv[++iarg];
			} else
				exitProg(ErrorArgument, "need the alphabet string after option a");
		}
		if(option['f']) {
			option['f'] = 0;
			if((iarg+1)<argc && sscanf(argv[iarg+1], "%c", &format) == 1)
				iarg++;
			else
				exitProg(ErrorArgument, "l or c is expected after -f");
		}
		if(option['e']) {
			option['e'] = 0;
			if((iarg+1)<argc && sscanf(argv[iarg+1], "%d", &sgn) == 1)
				iarg++;
			else
				exitProg(ErrorArgument, "need a number after -d");
		}
		if(option['s']) {
			option['s'] = 0;
			if((iarg+1)<argc && sscanf(argv[iarg+1], "%d", &nbSample) == 1)
				iarg++;
			else
				exitProg(ErrorArgument, "need a length after -l");
		}
		if(option['n']) {
			option['n'] = 0;
			if((iarg+1)<argc && sscanf(argv[iarg+1], "%d", &nbNew) == 1)
				iarg++;
			else
				exitProg(ErrorArgument, "need a length after -l");
		}
		if(option['u']) {
			option['u'] = 0;
			unif = 1;
		}
		if(option['p']) {
			option['p'] = 0;
			if((iarg+1)<argc) {
				motif = argv[++iarg];
			} else
				exitProg(ErrorArgument, "need a pattern after option p");
		}
		if(option['l']) {
			option['l'] = 0;
			if((iarg+1)<argc && sscanf(argv[iarg+1], "%d", &length) == 1)
				iarg++;
			else
				exitProg(ErrorArgument, "need a length after -l");
		}
		if(option['x']) {
			option['x'] = 0;
			if((iarg+1)<argc) {
				output = argv[++iarg];
			} else
				exitProg(ErrorArgument, "a file name is expected after -x");
		}			
		if(option['b']) {
			option['b'] = 0;
			if((iarg+1)<argc) {
				inputFileBernoulli = argv[++iarg];
			} else
				exitProg(ErrorArgument, "a file name is expected after -b");
		}
		if(option['h']) {
			printf("%s\n", HELPMESSAGE);
			exitProg(ExitOk, NULL);
		}
	}
	if(alphabet == NULL) {
		alphabet = ALPHABET;
		fprintf(stderr, "No alphabet provided. Use {a,b}.\n");;
	}
	cardinal = strlen(alphabet);
	if(inputFileBernoulli == NULL) {
		fprintf(stderr, "No Bernoulli file provided. Use uniform probabilities.\n");;
		bernoulli = getUniformBernoulli(cardinal);
	} else {
		FILE *fb;
		if(!(fb = fopen(inputFileBernoulli, "r"))) {
			fprintf(stderr, "Can't open file %s. Use uniform probabilities.\n",inputFileBernoulli);
			bernoulli = getUniformBernoulli(cardinal);
		} else {
			bernoulli = freadBernoulli(fb, alphabet);
			fclose(fb);
		}
	}
	bernoulliUnif = getUniformBernoulli(cardinal);
	if(unif)
		bernoulli2 = bernoulliUnif;
	if(output == NULL) {
		switch(format) {
			case 'l':
				sprintf(TMP, "%s.tex", OUTPUT_DEFAULT);
				break;
			case 'c':
			default:
				sprintf(TMP, "%s.csv", OUTPUT_DEFAULT);
		}
		output = TMP;
	}
	if(motif != NULL)
		length = strlen(motif);
	if((fo = fopen(output, "w"))) {
		char **col;
		int ncol, c;
		if(unif) {
			if(length<=4)
				ncol = nbFunct+2*(nbNew+1);
			else
				ncol = nbFunct+2*nbNew;
		} else {
			if(length<=4)
				ncol = nbFunct+nbNew+1;
			else
				ncol = nbFunct+nbNew;
		}
		col = (char**) malloc(ncol*sizeof(char*));
		c=0;
		for(i=0; i<nbFunct; i++)
			col[c++] = nameFunct[i];
		if(unif) {
			for(i=0; i<nbNew; i++) {
				sprintf(TMP, "%s %d%s", UNIFORM, i+1, HEURISTIC);
				col[c] = (char*) malloc((strlen(TMP)+1)*sizeof(char));
				strcpy(col[c], TMP);
				c++;
			}
			if(length <= 4) {
				sprintf(TMP, "%s %s", UNIFORM, FASTEST);
				col[c] = (char*) malloc((strlen(TMP)+1)*sizeof(char));
				strcpy(col[c], TMP);
				c++;
			}
		}
		for(i=0; i<nbNew; i++) {
			sprintf(TMP, "%d%s", i+1, HEURISTIC);
			col[c] = (char*) malloc((strlen(TMP)+1)*sizeof(char));
			strcpy(col[c], TMP);
			c++;
		}
		if(length <= 4)
			col[c++] = FASTEST;
		if(motif == NULL) {
			double **val;
			int nrow;
			char **row;
			if(nbSample <= 0) {
				nrow = (int) floor(pow((double) cardinal, length));
				val = (double**) malloc(nrow*sizeof(double));
				row = (char**) malloc(nrow*sizeof(char*));
				pattern = (TypeSymbol*) malloc(length*sizeof(TypeSymbol));
				for(p=0; p<length; p++)
					pattern[p] = 0;
				i = 0;
				do {
					row[i] = (char*) malloc((length+1)*sizeof(char));
					row[i][length] = '\0';
					for(p=0; p<length; p++)
						row[i][p] = alphabet[pattern[p]];
					val[i] = (double*) malloc(ncol*sizeof(double));
					computeLine(val[i], straFunct, nbFunct, nbNew, bernoulli, bernoulli2, cardinal, pattern, length);
					for(p=length-1; p>=0 && pattern[p]==cardinal-1; p--)
						pattern[p] = 0;
					if(p>=0)
						pattern[p]++;
					i++;
				} while(p >= 0);
				fprintTable(fo, col, row, val, nrow, ncol, 1, 1, format, sgn);
			} else {
				nrow = nbSample;
				val = (double**) malloc(nrow*sizeof(double));
				row = (char**) malloc(nrow*sizeof(char*));
				pattern = (TypeSymbol*) malloc(length*sizeof(TypeSymbol));
				for(i=0; i<nrow; i++) {
					switch(rtype) {
						case 'm':
							fillRandomText(pattern, length, bernoulli);
							break;
						case 'u':
						default:
							fillRandomText(pattern, length, bernoulliUnif);
					}
					row[i] = (char*) malloc((length+1)*sizeof(char));
					row[i][length] = '\0';
					for(p=0; p<length; p++)
						row[i][p] = alphabet[pattern[p]];
					val[i] = (double*) malloc(ncol*sizeof(double));
					computeLine(val[i], straFunct, nbFunct, nbNew, bernoulli, bernoulli2, cardinal, pattern, length);
				}
				fprintTable(fo, col, row, val, nrow, ncol, 1, 1, format, sgn);
			}
			for(i=0; i<nrow; i++) {
				free((void*)row[i]);
				free((void*)val[i]);
			}
			free((void*)row);
			free((void*)val);
		} else {
			double **val;
			int nrow;
			char **row;
			nrow = 1;
			val = (double**) malloc(nrow*sizeof(double));
			row = (char**) malloc(nrow*sizeof(char*));
			pattern = toSymbolSequence(motif, length, alphabet);
			row[0] = motif;
			val[0] = (double*) malloc(ncol*sizeof(double));
			computeLine(val[0], straFunct, nbFunct, nbNew, bernoulli, bernoulliUnif, cardinal, pattern, length);
			fprintTable(fo, col, row, val, nrow, ncol, 1, 1, format, sgn);
			free((void*)val[0]);
			free((void*)row);
			free((void*)val);
		}
		if(unif) {
			if(length <= 4) {
				for(i=0; i<=nbNew; i++)
					free((void*)col[nbFunct+i]);
				for(i=0; i<nbNew; i++)
					free((void*)col[nbFunct+nbNew+1+i]);
			} else {
				for(i=0; i<nbNew; i++)
					free((void*)col[nbFunct+i]);
				for(i=0; i<nbNew; i++)
					free((void*)col[nbFunct+nbNew+i]);
			}
		} else {
			for(i=0; i<nbNew; i++)
				free((void*)col[nbFunct+i]);
		}
		free((void*)col);
		free((void*)pattern);
		fclose(fo);
	}
	freeBernoulli(bernoulli);
	freeBernoulli(bernoulliUnif);
	return EXIT_SUCCESS;
}

