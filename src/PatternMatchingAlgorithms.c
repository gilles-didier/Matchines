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
#include <string.h>
#include <math.h>

#include "PatternMatchingAlgorithms.h"
#include "Utils.h"


static int *getBorder(TypeSymbol *pattern, int length);
static int *getStrictBorder(TypeSymbol *pattern, int length);
static int* getBoyerMooreBadChar(TypeCardinal cardinal, TypeSymbol *m, int l);
static void FJSpreQsBc(TypeCardinal cardinal, TypeSymbol *x, int m, int *qbc);
static void FJSpreKmp(TypeSymbol *x, int m, int *kmpNexy);
static void TVSBSpreBrBc(TypeCardinal cardinal, TypeSymbol *x, int m, int **brBc);
static void SApre(TypeCardinal cardinal, TypeSymbol *x, int m, unsigned int *S);
static void addPos(TypeState *st, int pos);
static int fillHASHqBadChar(TypeCardinal cardinal, TypeSymbol *pattern, int m, int *shift);


/*return st union {pos}*/
void addPos(TypeState *st, int pos) {
	*st = (*st) | (1<<pos);
}

TypeMatchingMachine *getNaiveMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	TypeMatchingMachine *matchine;
	int i;
	TypeCardinal a;
	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	matchine->cardinal = cardinal;
	matchine->size = length;
	matchine->trans = (TypeState**) malloc(matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(matchine->size*sizeof(int*));
	matchine->next = (int*) malloc(matchine->size*sizeof(int));
	matchine->sizeTerm = 1;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	for(i=0; i<length; i++) {
		matchine->trans[i] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[i] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[i] = i;
		for(a=0; a<pattern[i]; a++) {
//printf("a %u/%u\n", a, matchine->cardinal);
			matchine->trans[i][a] = 0;
			matchine->shift[i][a] = 1;
		}
		if(i<length-1) {
			matchine->trans[i][pattern[i]] = i+1;
			matchine->shift[i][pattern[i]] = 0;
		} else {
			matchine->trans[i][pattern[i]] = 0;
			matchine->shift[i][pattern[i]] = 1;
		}
		for(a=pattern[i]+1; a<matchine->cardinal; a++) {
			matchine->trans[i][a] = 0;
			matchine->shift[i][a] = 1;
		}
	}
	matchine->term[0] = length-1;
	matchine->xterm[0] = pattern[matchine->next[matchine->term[0]]];
	return matchine;
}


TypeMatchingMachine *getMatchAutomataMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	TypeMatchingMachine *matchine;
	int i;
	TypeCardinal a;
	TypeLattice *lattice = getLattice(cardinal, pattern, length);
	TypeState *tabSt;
	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	matchine->cardinal = cardinal;
	matchine->size = length;
	matchine->trans = (TypeState**) malloc(matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(matchine->size*sizeof(int*));
	matchine->next = (int*) malloc(matchine->size*sizeof(int));
	matchine->sizeTerm = 1;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	tabSt = (TypeState*) malloc((length)*sizeof(TypeState));
	tabSt[0] = 0;
	for(i=1; i<length; i++) {
		tabSt[i] = tabSt[i-1];
		addPos(&(tabSt[i]), i-1);
	}
	for(i=0; i<length; i++) {
		matchine->trans[i] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[i] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[i] = i;
		for(a=0; a<pattern[i]; a++) {
			matchine->trans[i][a] = MAX(0, i-lattice->shift[tabSt[i]][0][a]+1);
			matchine->shift[i][a] = lattice->shift[tabSt[i]][0][a];
		}
		if(i<length-1) {
			matchine->trans[i][pattern[i]] = i+1;
			matchine->shift[i][pattern[i]] = 0;
		} else {
			matchine->trans[i][pattern[i]] = MAX(0, i-lattice->shift[tabSt[i]][0][a]+1);
			matchine->shift[i][pattern[i]] = lattice->shift[tabSt[i]][0][a];
		}			
		for(a=pattern[i]+1; a<matchine->cardinal; a++) {
			matchine->trans[i][a] = MAX(0, i-lattice->shift[tabSt[i]][0][a]+1);
			matchine->shift[i][a] = lattice->shift[tabSt[i]][0][a];
		}
	}
	matchine->term[0] = length-1;
	matchine->xterm[0] = pattern[matchine->next[matchine->term[0]]];
	free((void*)tabSt);
	freeLattice(lattice);
	return matchine;
}	


int *getBorder(TypeSymbol *pattern, int length) {
	int i, j, *border;
	border = (int*) malloc((length+1)*sizeof(int));
	i = 0; 
	j = -1;
	border[0] = -1;
	while(i<length) {
		while(j>=0 && pattern[i] != pattern[j])
			j = border[j];
		border[++i] = ++j;
	}
	return border;
}


TypeMatchingMachine *getMorrisPrattMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	TypeMatchingMachine *matchine;
	int i, *border;
	TypeCardinal a;
	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	matchine->cardinal = cardinal;
	matchine->size = length;
	matchine->trans = (TypeState**) malloc(matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(matchine->size*sizeof(int*));
	matchine->next = (int*) malloc(matchine->size*sizeof(int));
	matchine->sizeTerm = 1;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	border = getBorder(pattern, length);
	for(i=0; i<length-1; i++) {
		matchine->trans[i] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[i] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[i] = i;
		for(a=0; a<pattern[i]; a++) {
			matchine->trans[i][a] = MAX(0, border[i]);
			matchine->shift[i][a] = i-border[i];
		}
		matchine->trans[i][pattern[i]] = i+1;
		matchine->shift[i][pattern[i]] = 0;
		for(a=pattern[i]+1; a<matchine->cardinal; a++) {
			matchine->trans[i][a] = MAX(0, border[i]);
			matchine->shift[i][a] = i-border[i];
		}
	}
	matchine->trans[length-1] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
	matchine->shift[length-1] = (int*) malloc(matchine->cardinal*sizeof(int));
	matchine->next[length-1] = length-1;
	for(a=0; a<pattern[length-1]; a++) {
		matchine->trans[length-1][a] = MAX(0, border[length-1]);
		matchine->shift[length-1][a] = length-1-border[i];
	}
	matchine->trans[length-1][pattern[length-1]] = border[length];
	matchine->shift[length-1][pattern[length-1]] = length-border[length];
	for(a=pattern[length-1]+1; a<matchine->cardinal; a++) {
		matchine->trans[length-1][a] = MAX(0, border[length-1]);
		matchine->shift[length-1][a] = length-1-border[i];
	}
	matchine->term[0] = length-1;
	matchine->xterm[0] = pattern[matchine->next[matchine->term[0]]];
	free((void*)border);
	return matchine;
}


int *getStrictBorder(TypeSymbol *pattern, int length) {
	int i, j, *border;
	border = (int*) malloc((length+1)*sizeof(int));
	i = 0; 
	j = -1;
	border[0] = -1;
	while(i<length) {
		while(j>=0 && pattern[i] != pattern[j])
			j = border[j];
		border[++i] = ++j;
		if(i<length && pattern[i] == pattern[j])
			border[i] = border[j];
	}
	return border;
}


TypeMatchingMachine *getKnuthMorrisPrattMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	TypeMatchingMachine *matchine;
	int i, *border;
	TypeCardinal a;
	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	matchine->cardinal = cardinal;
	matchine->size = length;
	matchine->trans = (TypeState**) malloc(matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(matchine->size*sizeof(int*));
	matchine->next = (int*) malloc(matchine->size*sizeof(int));
	matchine->sizeTerm = 1;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	border = getStrictBorder(pattern, length);	
	for(i=0; i<length-1; i++) {
		matchine->trans[i] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[i] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[i] = i;
		for(a=0; a<pattern[i]; a++) {
			matchine->trans[i][a] = MAX(0, border[i]);
			matchine->shift[i][a] = i-border[i];
		}
		matchine->trans[i][pattern[i]] = i+1;
		matchine->shift[i][pattern[i]] = 0;
		for(a=pattern[i]+1; a<matchine->cardinal; a++) {
			matchine->trans[i][a] = MAX(0, border[i]);
			matchine->shift[i][a] = i-border[i];
		}
	}
	matchine->trans[length-1] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
	matchine->shift[length-1] = (int*) malloc(matchine->cardinal*sizeof(int));
	matchine->next[length-1] = length-1;
	for(a=0; a<pattern[length-1]; a++) {
		matchine->trans[length-1][a] = MAX(0, border[length-1]);
		matchine->shift[length-1][a] = length-1-border[length-1];
	}
	matchine->trans[length-1][pattern[length-1]] = border[length];
	matchine->shift[length-1][pattern[length-1]] = length-border[length];
	for(a=pattern[length-1]+1; a<matchine->cardinal; a++) {
		matchine->trans[length-1][a] = MAX(0, border[length-1]);
		matchine->shift[length-1][a] = length-1-border[length-1];
	}
	matchine->term[0] = length-1;
	matchine->xterm[0] = pattern[matchine->next[matchine->term[0]]];
	free((void*)border);
	return matchine;
}


int* getBoyerMooreBadChar(TypeCardinal cardinal, TypeSymbol *m, int l) {
   int i;
   int *bmBc = (int*) malloc(cardinal*sizeof(int));
   for(i=0; i<cardinal; ++i)
      bmBc[i] = l;
   for(i=0; i<l-1; ++i)
      bmBc[m[i]] = l-i-1;
   return bmBc;
}


TypeMatchingMachine *getHorspoolMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	TypeMatchingMachine *matchine;
	int i, *boyerMooreBadChar, last;
	TypeCardinal a;
	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	boyerMooreBadChar = getBoyerMooreBadChar(cardinal, pattern, length);
	matchine->cardinal = cardinal;
	matchine->size = length;
	matchine->trans = (TypeState**) malloc(matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(matchine->size*sizeof(int*));
	matchine->next = (int*) malloc(matchine->size*sizeof(int));
	matchine->sizeTerm = 1;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	matchine->trans[0] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
	matchine->shift[0] = (int*) malloc(matchine->cardinal*sizeof(int));
	matchine->next[0] = length-1;
	for(a=0; a<pattern[length-1]; a++) {
		matchine->trans[0][a] = 0;
		matchine->shift[0][a] = boyerMooreBadChar[a];
	}
	if(length > 1) {
		matchine->trans[0][pattern[length-1]] = 1;
		matchine->shift[0][pattern[length-1]] = 0;
	} else {
		matchine->trans[0][pattern[length-1]] = 0;
		matchine->shift[0][pattern[length-1]] = 1;
	}
	for(a=pattern[length-1]+1; a<matchine->cardinal; a++) {
		matchine->trans[0][a] = 0;
		matchine->shift[0][a] = boyerMooreBadChar[a];
	}
	last = boyerMooreBadChar[pattern[length-1]];
	for(i=1; i<length; i++) {
		matchine->trans[i] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[i] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[i] = length-i-1;
		for(a=0; a<pattern[length-1-i]; a++) {
			matchine->trans[i][a] = 0;
			matchine->shift[i][a] = last;
		}
		if(i<length-1) {
			matchine->trans[i][pattern[length-1-i]] = i+1;
			matchine->shift[i][pattern[length-1-i]] = 0;
		} else {
			matchine->trans[i][pattern[length-1-i]] = 0;
			matchine->shift[i][pattern[length-1-i]] = last;
		}
		for(a=pattern[length-1-i]+1; a<matchine->cardinal; a++) {
			matchine->trans[i][a] = 0;
			matchine->shift[i][a] = last;
		}
	}
	matchine->term[0] = length-1;
	matchine->xterm[0] = pattern[matchine->next[matchine->term[0]]];
	free((void*)boyerMooreBadChar);
	return matchine;
}



// HASHq Lecroq Wu Manber
#define WSIZE 255
#define RANK 3


int fillHASHqBadChar(TypeCardinal cardinal, TypeSymbol *pattern, int m, int *shift) {
	int i, r, sh1;
	TypeSymbol h;
	
	for(i=0; i<WSIZE; ++i)
		shift[i] = m-RANK+1;
	h = pattern[0];
	for(r=1; r<RANK; r++)
		h = ((h<<1) + pattern[r]);
	shift[h] = m-RANK;
	for(i=RANK; i<m-1; ++i) {
		h = pattern[i-RANK+1];
		for(r=1; r<RANK; r++)
			h = ((h<<1)+pattern[r+i-RANK+1]);
		shift[h%WSIZE] = m-1-i;
	}
	h = pattern[m-RANK];
	for(r=1; r<RANK; r++)
		h = ((h<<1) + pattern[r+m-RANK]);
	sh1 = shift[h%WSIZE];
	shift[h%WSIZE] = 0;
	if(sh1==0)
		sh1 = 1;
	return sh1;
}


TypeMatchingMachine *getHASHqMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	TypeMatchingMachine *matchine;
	int i, k, sh1, *HASHqBadChar;
	TypeCardinal a;
	int debut, fin;
	
	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	HASHqBadChar = (int *) malloc (WSIZE * sizeof(int)); 
	sh1 = fillHASHqBadChar(cardinal, pattern, length, HASHqBadChar);
	matchine->cardinal = cardinal;
	matchine->size = (pow(cardinal, RANK)-1)/(cardinal-1)+length;
	matchine->trans = (TypeState**) malloc(matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(matchine->size*sizeof(int*));
	matchine->next = (int*) malloc(matchine->size*sizeof(int));
	matchine->sizeTerm = 1;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	debut = 0;
	fin = 1;
	for (i=0; i<RANK-1; i++) {
		for (k=debut; k<fin ; k++) {
			matchine->trans[k] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
			matchine->shift[k] = (int*) malloc(matchine->cardinal*sizeof(int));
			matchine->next[k] = length-RANK+i;
			for(a=0; a<matchine->cardinal; a++) {
				matchine->trans[k][a] = fin + ((k-debut)<<1) + a;
				matchine->shift[k][a] = 0;
			}
		}
		debut = fin;
		fin = debut + pow(cardinal, i+1);
	}
	for (k=debut; k<fin ; k++) {
		matchine->trans[k] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[k] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[k] = length-1;
		for(a=0; a<matchine->cardinal; a++) {
			matchine->shift[k][a] = HASHqBadChar[(((k-debut)<<1) + a)%WSIZE];
			if (matchine->shift[k][a]!=0) {
				matchine->trans[k][a] = 0;
			} else {
				matchine->trans[k][a] = fin;
			}
		}
	}
	for (k=0; k<length; k++) {
		matchine->trans[fin+k] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[fin+k] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[fin+k] = k;
		for(a=0; a<pattern[k]; a++) {
			matchine->trans[fin+k][a] = 0;
			matchine->shift[fin+k][a] = sh1;
		}
		if (k< length-1) {
			matchine->trans[fin+k][pattern[k]] = fin+k+1;
			matchine->shift[fin+k][pattern[k]] = 0;
		} else {
			matchine->trans[fin+k][pattern[k]] = 0;
			matchine->shift[fin+k][pattern[k]] = sh1;
		}
		for(a=pattern[k]+1; a<matchine->cardinal; a++) {
			matchine->trans[fin+k][a] = 0;
			matchine->shift[fin+k][a] = sh1;
		}
	}
	matchine->term[0] = fin+length-1;
	matchine->xterm[0] = pattern[matchine->next[matchine->term[0]]];
	free((void*)HASHqBadChar);
	return matchine;
}





TypeMatchingMachine *getQuickSearchMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	TypeMatchingMachine *matchine;
	int i, *quickSearchBadChar;
	TypeCardinal a;
	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	quickSearchBadChar = (int*) malloc(cardinal * sizeof(int));
	FJSpreQsBc(cardinal, pattern, length, quickSearchBadChar);

	matchine->cardinal = cardinal;
	matchine->size = length+1;
	matchine->trans = (TypeState**) malloc(matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(matchine->size*sizeof(int*));
	matchine->next = (int*) malloc(matchine->size*sizeof(int));
	matchine->sizeTerm = 1;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	for(i=0; i<length; i++) {
		matchine->trans[i] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[i] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[i] = i;
		for(a=0; a<pattern[i]; a++) {
			matchine->trans[i][a] = length;
			matchine->shift[i][a] = 0;
		}
		if(i<length-1) {
			matchine->trans[i][pattern[i]] = i+1;
			matchine->shift[i][pattern[i]] = 0;
		} else {
			matchine->trans[i][pattern[i]] = length;
			matchine->shift[i][pattern[i]] = 0;
		}
		for(a=pattern[i]+1; a<matchine->cardinal; a++) {
			matchine->trans[i][a] = length;
			matchine->shift[i][a] = 0;
		}
	}
	matchine->trans[length] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
	matchine->shift[length] = (int*) malloc(matchine->cardinal*sizeof(int));
	matchine->next[length] = length;
	for(a=0; a<matchine->cardinal; a++) {
		matchine->trans[length][a] = 0;
		matchine->shift[length][a] = quickSearchBadChar[a];
	}
	matchine->term[0] = length-1;
	matchine->xterm[0] = pattern[matchine->next[matchine->term[0]]];
	free((void*)quickSearchBadChar);
	return matchine;
}


void FJSpreQsBc(TypeCardinal cardinal, TypeSymbol *x, int m, int *qbc) {
   int i;
   for(i=0; i<cardinal; i++)
		qbc[i]=m+1;
   for(i=0; i<m; i++)
		qbc[x[i]]=m-i;
}


void FJSpreKmp(TypeSymbol *x, int m, int *kmpNexy) {
   int i, j;
   i = 0;
   j = kmpNexy[0] = -1;
   while(i < m) {
      while(j>-1 && x[i] != x[j])
         j = kmpNexy[j];
      i++;
      j++;
      if(i<m && x[i] == x[j])
         kmpNexy[i] = kmpNexy[j];
      else
         kmpNexy[i] = j;
   }
}


TypeMatchingMachine *getFJSMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	TypeMatchingMachine *matchine;
	int i, *qsbc, *kmp;
	TypeCardinal a;
	
	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	qsbc = (int*) malloc(cardinal*sizeof(int));
	kmp = (int*) malloc((length+1)*sizeof(int));
	/* Preprocessing */
	FJSpreQsBc(cardinal, pattern, length, qsbc);

	FJSpreKmp(pattern, length, kmp);
	matchine->cardinal = cardinal;
	matchine->size = length+2;
	matchine->trans = (TypeState**) malloc(matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(matchine->size*sizeof(int*));
	matchine->next = (int*) malloc(matchine->size*sizeof(int));
	matchine->sizeTerm = 1;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	matchine->trans[0] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
	matchine->shift[0] = (int*) malloc(matchine->cardinal*sizeof(int));
	matchine->next[0] = length-1;
	for(a=0; a<pattern[length-1]; a++) {
		matchine->trans[0][a] = 1;
		matchine->shift[0][a] = 1;
	}
	matchine->trans[0][pattern[length-1]] = 2;
	matchine->shift[0][pattern[length-1]] = 0;
	for(a=pattern[length-1]+1; a<matchine->cardinal; a++) {
		matchine->trans[0][a] = 1;
		matchine->shift[0][a] = 1;
	}
	
	matchine->trans[1] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
	matchine->shift[1] = (int*) malloc(matchine->cardinal*sizeof(int));
	matchine->next[1] = length-1;
	for(a=0; a<matchine->cardinal; a++) {
		matchine->trans[1][a] = 0;
		matchine->shift[1][a] = qsbc[a]-1;
	}
	for(i=0; i<length; i++) {
		matchine->trans[i+2] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[i+2] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[i+2] = i;
		for(a=0; a<pattern[i]; a++) {
			matchine->trans[i+2][a] = 0;
			matchine->shift[i+2][a] = i-kmp[i];
		}
		if(i<length-1) {
			matchine->trans[i+2][pattern[i]] = i+3;
			matchine->shift[i+2][pattern[i]] = 0;
		} else {
			matchine->trans[i+2][pattern[i]] = 0;
			matchine->shift[i+2][pattern[i]] = length-kmp[length];
		}			
		for(a=pattern[i]+1; a<matchine->cardinal; a++) {
			matchine->trans[i+2][a] = 0;
			matchine->shift[i+2][a] = i-kmp[i];
		}
	}
	matchine->term[0] = length+1;
	matchine->xterm[0] = pattern[matchine->next[matchine->term[0]]];
	free((void*)kmp);
	free((void*)qsbc);
	return matchine;
}


void TVSBSpreBrBc(TypeCardinal cardinal, TypeSymbol *x, int m, int **brBc) {
   int i;
   TypeCardinal a, b;
   for(a=0; a<cardinal; ++a)
      for(b=0; b<cardinal; ++b)
         brBc[a][b] = m+2;
   for(a=0; a<cardinal; ++a)
      brBc[a][x[0]] = m+1;
   for(i=0; i<m-1; ++i)
      brBc[x[i]][x[i+1]] = m - i;
   for(a=0; a<cardinal; ++a)
      brBc[x[m-1]][a] = 1;
}

TypeMatchingMachine *getTVSBMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	TypeMatchingMachine *matchine;
	int i, **BrBc;
	TypeCardinal a, b;
	
	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	BrBc = (int**) malloc(cardinal*sizeof(int*));
	for(i=0; i<cardinal; i++)
		BrBc[i] = (int*) malloc(cardinal*sizeof(int));
// Preprocessing 
	TVSBSpreBrBc(cardinal, pattern, length, BrBc);

	matchine->cardinal = cardinal;
	matchine->size = length+cardinal+1;
	matchine->trans = (TypeState**) malloc(matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(matchine->size*sizeof(int*));
	matchine->next = (int*) malloc(matchine->size*sizeof(int));
	matchine->sizeTerm = 1;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	matchine->trans[0] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
	matchine->shift[0] = (int*) malloc(matchine->cardinal*sizeof(int));
	matchine->next[0] = length-1;
	for(a=0; a<pattern[matchine->next[0]]; a++) {
		matchine->trans[0][a] = length;
		matchine->shift[0][a] = 0;
	}
	matchine->trans[0][pattern[matchine->next[0]]] = length-1;
	matchine->shift[0][pattern[matchine->next[0]]] = 0;
	for(a=pattern[matchine->next[0]]+1; a<matchine->cardinal; a++) {
		matchine->trans[0][a] = length;
		matchine->shift[0][a] = 0;
	}
	for(i=length-2; i>=1; i--) {
		matchine->trans[i] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[i] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[i] = i;
		for(a=0; a<pattern[matchine->next[i]]; a++) {
			matchine->trans[i][a] = length;
			matchine->shift[i][a] = 0;
		}
		if(i>1) {
			matchine->trans[i][pattern[matchine->next[i]]] = i-1;
			matchine->shift[i][pattern[matchine->next[i]]] = 0;
		} else {
			matchine->trans[i][pattern[matchine->next[i]]] = length;
			matchine->shift[i][pattern[matchine->next[i]]] = 0;
		}		
		for(a=pattern[matchine->next[i]]+1; a<matchine->cardinal; a++) {
			matchine->trans[i][a] = length;
			matchine->shift[i][a] = 0;
		}
		
	}
	matchine->trans[length-1] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
	matchine->shift[length-1] = (int*) malloc(matchine->cardinal*sizeof(int));
	matchine->next[length-1] = 0;
	for(a=0; a<pattern[matchine->next[length-1]]; a++) {
		matchine->trans[length-1][a] = length;
		matchine->shift[length-1][a] = 0;
	}
	matchine->trans[length-1][pattern[matchine->next[length-1]]] = length-2;
	matchine->shift[length-1][pattern[matchine->next[length-1]]] = 0;
	for(a=pattern[matchine->next[length-1]]+1; a<matchine->cardinal; a++) {
		matchine->trans[length-1][a] = length;
		matchine->shift[length-1][a] = 0;
	}
	matchine->trans[length] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
	matchine->shift[length] = (int*) malloc(matchine->cardinal*sizeof(int));
	matchine->next[length] = length;
	for(a=0; a<matchine->cardinal; a++) {
		matchine->trans[length][a] = length+a+1;
		matchine->shift[length][a] = 0;
	}
	for(b=0; b<matchine->cardinal; b++) {
		matchine->trans[length+b+1] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[length+b+1] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[length+b+1] = length+1;
		for(a=0; a<matchine->cardinal; a++) {
			matchine->trans[length+b+1][a] = 0;
			matchine->shift[length+b+1][a] = BrBc[b][a];
		}
	}
	matchine->term[0] = 1;
	matchine->xterm[0] = pattern[matchine->next[matchine->term[0]]];
	for(i=0; i<cardinal; i++)
		free((void*)BrBc[i]);
	free((void*)BrBc);
	return matchine;
}

void SApre(TypeCardinal cardinal, TypeSymbol *x, int m, unsigned int *S) { 
	unsigned int j; 
	int i;
	TypeCardinal a;
	for(a=0; a<cardinal; ++a)
		S[a] = 0; 
	for(i=0, j=1; i<m; ++i, j <<= 1)
		S[x[i]] |= j; 
} 



TypeMatchingMachine *getSAMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	TypeMatchingMachine *matchine;
	int i, size, *state, *stack, is;
	TypeCardinal a;
	unsigned int *S;
	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	S = (unsigned int*) malloc(cardinal*sizeof(int));
	/* Preprocessing */
	SApre(cardinal, pattern, length, S);
	matchine->cardinal = cardinal;
	size = 1<<length;
	matchine->size = 1;
	matchine->trans = (TypeState**) malloc(size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(size*sizeof(int*));
	matchine->next = (int*) malloc(size*sizeof(int));
	matchine->sizeTerm = 1;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	stack = (int*) malloc(size*sizeof(int));
	state = (int*) malloc(size*sizeof(int));
	for(i=0; i<size; i++)
		state[i] = SPECIAL;
	state[0] = 0;
	stack[0] = 0;
	is = 1;
	while(is>0) {
		unsigned int D, E;
		int st;
		D = stack[--is];
		st = state[D];
		matchine->trans[st] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[st] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[st] = 0;
		for(a=0; a<matchine->cardinal; a++) {
			E = ((D<<1) | 1) & S[a];
			if(state[E] == SPECIAL) {
				state[E] = matchine->size++;
				stack[is++] = E;
			}
			matchine->trans[st][a] = state[E];
			matchine->shift[st][a] = 1;
		}
	}
	matchine->trans = (TypeState**) realloc(matchine->trans, matchine->size*sizeof(TypeState*));
	matchine->term[0] = 0;	
	matchine->xterm[0] = pattern[matchine->next[matchine->term[0]]];
	free((void*)S);
	free((void*)state);
	free((void*)stack);
	return matchine;
}


TypeMatchingMachine *getEBOMMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	TypeMatchingMachine *matchine;
	int i, st, j, p, k, q=1, iMinus1;
	TypeCardinal a, c;
	int *S, **FT, **trans;

	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
   // Preprocessing
	S = (int*) malloc((length+1)*sizeof(int));
	trans = (int**) malloc((length+2)*sizeof(int*));
	FT = (int**) malloc(cardinal*sizeof(int*));
	for (i=0; i<=length+1; i++)
		trans[i] = (int*) malloc(sizeof(int)*cardinal);
	for (i=0; i<=length+1; i++)
		for (j=0; j<cardinal; j++)
			trans[i][j]=SPECIAL;
	S[length] = length + 1;
	for (i=length; i>0; i--) {
		iMinus1 = i-1;
		c = pattern[iMinus1];
		trans[i][c] = iMinus1;
		p = S[i];
		while (p <= length && (q = trans[p][c]) ==  SPECIAL) {
			trans[p][c] = iMinus1;
			p = S[p];
		}
		S[iMinus1] = (p == length + 1 ? length : q);
	}
	/* Construct the FirstTransition table */
	for (i=0; i<cardinal; i++) {
		q = trans[length][i];
		FT[i] =(int*) malloc(cardinal*sizeof(int));
		for (j=0; j<cardinal; j++)
			if (q != SPECIAL)
				FT[i][j] = trans[q][j];
			else
				FT[i][j] = SPECIAL;
	}
	matchine->cardinal = cardinal;
	matchine->size = cardinal+length*length+2;
	matchine->trans = (TypeState**) malloc(matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(matchine->size*sizeof(int*));
	matchine->next = (int*) malloc(matchine->size*sizeof(int));
	matchine->sizeTerm = 1;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	matchine->trans[0] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
	matchine->shift[0] = (int*) malloc(matchine->cardinal*sizeof(int));
	matchine->next[0] = length-1;
	for(a=0; a<cardinal; a++) {
		matchine->trans[0][a] = a+1;
		matchine->shift[0][a] = 0;
	}
	for(st=1; st<=cardinal; st++) {
		matchine->trans[st] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[st] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[st] = length-2;
		for(a=0; a<cardinal; a++) {
			if(FT[st-1][a] == SPECIAL) {
				matchine->trans[st][a] = 0;
				matchine->shift[st][a] = length-1;
			} else {
				matchine->trans[st][a] = (FT[st-1][a]*length+(length-3))+cardinal+1;
				matchine->shift[st][a] = 0;
			}
		}
	}
	
	for(q=0; q<length; q++) {
		for(k=1; k<length; k++) {
			st = (q*length)+k;
			matchine->trans[cardinal+1+st] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
			matchine->shift[cardinal+1+st] = (int*) malloc(matchine->cardinal*sizeof(int));
			matchine->next[cardinal+1+st] = k;
			for(a=0; a<cardinal; a++) {
				if(trans[q][a] == SPECIAL) {
					matchine->trans[cardinal+1+st][a] = 0;
					matchine->shift[cardinal+1+st][a] = k+1;
				} else {
					matchine->trans[cardinal+1+st][a] = trans[q][a]*length+(k-1)+(cardinal+1);
					matchine->shift[cardinal+1+st][a] = 0;
				}
			}
		}
		st = q*length;
		matchine->trans[cardinal+1+st] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[cardinal+1+st] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[cardinal+1+st] = 0;
		for(a=0; a<cardinal; a++) {
			if(trans[q][a] == SPECIAL) {
				matchine->trans[cardinal+1+st][a] = 0;
				matchine->shift[cardinal+1+st][a] = 1;
			} else {
				matchine->trans[cardinal+1+st][a] = cardinal+length*length+1;
				matchine->shift[cardinal+1+st][a] = 0;
			}
		}
	}
	matchine->trans[cardinal+length*length+1] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
	matchine->shift[cardinal+length*length+1] = (int*) malloc(matchine->cardinal*sizeof(int));
	matchine->next[cardinal+length*length+1] = 0;
	for(a=0; a<cardinal; a++) {
		matchine->trans[cardinal+length*length+1][a] = 0;
		matchine->shift[cardinal+length*length+1][a] = 1;
	}
	matchine->term[0] = 0;
	matchine->xterm[0] = pattern[matchine->next[matchine->term[0]]];
	free((void*)S);
	for (i=0; i<length+2; i++)
		free((void*)trans[i]);
	free((void*)trans);
	for (a=0; a<cardinal; a++)
		free((void*)FT[	a]);
	free((void*)FT);
	return matchine;
}

TypeMatchingMachine *getBOMMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	TypeMatchingMachine *matchine;
	int i, st, j, p, q=1;
	TypeCardinal a, c;
	int *S;
	int **trans;
	int iMinus1;
	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
   // Preprocessing
	S = (int*) malloc((length+1)*sizeof(int));
	trans = (int**) malloc((length+2)*sizeof(int*));
	for (i=0; i<=length+1; i++)
		trans[i] = (int*) malloc(sizeof(int)*cardinal);
	for (i=0; i<=length+1; i++)
		for (j=0; j<cardinal; j++)
			trans[i][j]=SPECIAL;
	S[length] = length + 1;
	for (i=length; i>0; i--) {
		iMinus1 = i-1;
		c = pattern[iMinus1];
		trans[i][c] = iMinus1;
		p = S[i];
		while (p <= length && (q = trans[p][c]) ==  SPECIAL) {
			trans[p][c] = iMinus1;
			p = S[p];
		}
		S[iMinus1] = (p == length + 1 ? length : q);
	}
	matchine->cardinal = cardinal;
	matchine->size = length;
	matchine->trans = (TypeState**) malloc(matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(matchine->size*sizeof(int*));
	matchine->next = (int*) malloc(matchine->size*sizeof(int));
	matchine->sizeTerm = 1;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	matchine->trans[0] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
	matchine->shift[0] = (int*) malloc(matchine->cardinal*sizeof(int));
	matchine->next[0] = length-1;
	for(a=0; a<cardinal; a++) {
		if(trans[length][a] == SPECIAL) {
			matchine->trans[0][a] = 0;
			matchine->shift[0][a] = length;
		} else {
			matchine->trans[0][a] = trans[length][a];
			matchine->shift[0][a] = 0;
		}
	}
	matchine->trans[1] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
	matchine->shift[1] = (int*) malloc(matchine->cardinal*sizeof(int));
	matchine->next[1] = 0;
	for(a=0; a<cardinal; a++) {
		matchine->trans[1][a] = 0;
		matchine->shift[1][a] = 1;
	}
	for(st=2; st<length; st++) {
		matchine->trans[st] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		matchine->shift[st] = (int*) malloc(matchine->cardinal*sizeof(int));
		matchine->next[st] = st-1;
		for(a=0; a<cardinal; a++) {
			if(trans[st][a] == SPECIAL) {
				matchine->trans[st][a] = 0;
				matchine->shift[st][a] = st;
			} else {
				matchine->trans[st][a] = trans[st][a];
				matchine->shift[st][a] = 0;
			}
		}
	}
	matchine->term[0] = 1;
	matchine->xterm[0] = pattern[matchine->next[matchine->term[0]]];
	free((void*)S);
	for (i=0; i<=length+1; i++)
		free((void*)trans[i]);
	free((void*)trans);
	return matchine;
}
