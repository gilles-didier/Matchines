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
#include <limits.h>
#include <math.h>

#include "Lattice.h"
#include "Utils.h"
#include "PowerSet.h"

#define ST_FREE -1
#define ST_STACKED -2
#define INC_STRATEGY_TMP_SIZE_DICT 1000

static int isIn(int pos, TypeState st);
static int getLastBorder(TypeSymbol *pattern, int length);
static int getTransIndex(int pos, TypeState st);
static void addPos(TypeState *st, int pos);
static void remPos(TypeState *st, int pos);


/*return st union {pos}*/
void addPos(TypeState *st, int pos) {
	*st = (*st) | (1<<pos);
}

/*return st minus {pos}*/
void remPos(TypeState *st, int pos) {
	*st = (*st) & ~(1<<pos);
}

/*return 1 if pos in st*/
int isIn(int pos, TypeState st) {
	return (st & (1<<pos)) != 0;
}

/*return the set/state of the set/table set*/
void setToState(int *set, int size, TypeState *st) {
	int i;
	*st = 0;
	for(i=0; i<size; i++)
		addPos(st, set[i]);
}

/*fill the set/table set with the set/state st*/
void stateToSet(TypeState st, int length, int *set, int *size) {
	int i;
	*size = 0;
	for(i=0; i<length; i++)
		if(isIn(i, st))
			set[(*size)++] = i;
}

/*return the cardinal of st*/
int stateCardinal(TypeState st, int length) {
	int i, card = 0;
	for(i=0; i<length; i++)
		if(isIn(i, st))
			card++;
	return card;
}

/*return the greater proper suffix of pattern which a prefix*/
int getLastBorder(TypeSymbol *pattern, int length) {
	int i, j , k, *border;
	border = (int*) malloc((length+1)*sizeof(int));
	i = 0; 
	j = -1;
	border[0] = -1;
	while(i<length) {
		while(j>=0 && pattern[i] != pattern[j])
			j = border[j];
		border[++i] = ++j;
	}
	k = border[length];
	free((void*)border);
	return k;
}

/*return the order of pos among the complementary of st*/
int getTransIndex(int pos, TypeState st) {
	int ind=0, i;
	for(i=0; i<pos; i++)
		if(!isIn(i, st))
			ind++;
	return ind;
}

/*return the full lattice of pattern*/
TypeLattice *getLattice(TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	TypeCardinal a;
	TypeLattice *lattice;
	TypeState size, stMatch;
	int **occPrec, sing, pos, level, i, j, k, *set, *comp, levelComp, last, shMatch;
	
	lattice = (TypeLattice*) malloc(sizeof(TypeLattice));
	lattice->length = length;
	lattice->cardinal = cardinal;	
	occPrec = (int**) malloc(lattice->length*sizeof(int*));
	for(pos=0; pos<lattice->length; pos++) {
		occPrec[pos] = (int*) malloc(lattice->cardinal*sizeof(int));
	}
	for(a=0; a<lattice->cardinal; a++)
		occPrec[0][a] = -1;
	for(pos=1; pos<lattice->length; pos++) {
		for(a=0; a<pattern[pos-1]; a++)
			occPrec[pos][a] = occPrec[pos-1][a];
		occPrec[pos][pattern[pos-1]] = pos-1;
		for(a=pattern[pos-1]+1; a<lattice->cardinal; a++)
			occPrec[pos][a] = occPrec[pos-1][a];
	}
	size = (1<<lattice->length) - 1;
	lattice->trans = (TypeState***) malloc(size*sizeof(TypeState**));
	lattice->shift = (int***) malloc(size*sizeof(int**));
	lattice->trans[0] = (TypeState**) malloc(lattice->length*sizeof(TypeState*));
	lattice->shift[0] = (int**) malloc(lattice->length*sizeof(int*));
	for(pos=0; pos<lattice->length; pos++) {
		TypeState dest;
		lattice->trans[0][pos] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState));
		lattice->shift[0][pos] = (int*) malloc(lattice->cardinal*sizeof(int));
		for(a=0; a<pattern[pos]; a++) {
			dest = 0;
			if(occPrec[pos][a] >= 0)
				addPos(&dest, occPrec[pos][a]);
			lattice->trans[0][pos][a] = dest;
			lattice->shift[0][pos][a] = pos-occPrec[pos][a];
		}
		dest = 0;
		addPos(&dest, pos);
		lattice->trans[0][pos][pattern[pos]] = dest;
		lattice->shift[0][pos][pattern[pos]] = 0;
		for(a=pattern[pos]+1; a<lattice->cardinal; a++) {
			dest = 0;
			if(occPrec[pos][a] >= 0)
				addPos(&dest, occPrec[pos][a]);
			lattice->trans[0][pos][a] = dest;
			lattice->shift[0][pos][a] = pos-occPrec[pos][a];
		}
	}
	last = getLastBorder(pattern, length);
	stMatch=0;
	for(i=0; i<last; i++)
		addPos(&stMatch, i);
	shMatch = lattice->length-last;
	for(sing=0; sing<lattice->length; sing++) {
		TypeState st=0;
		addPos(&st, sing);
		lattice->trans[st] = (TypeState**) malloc((lattice->length-1)*sizeof(TypeState*));
		lattice->shift[st] = (int**) malloc((lattice->length-1)*sizeof(int*));
		for(pos=0; pos<sing; pos++) {
			TypeState dest, sty;
			int shift;
			lattice->trans[st][pos] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState*));
			lattice->shift[st][pos] = (int*) malloc(lattice->cardinal*sizeof(int*));
			for(a=0; a<pattern[pos]; a++) {
				shift = lattice->shift[0][pos][a];
				sty = lattice->trans[0][pos][a];
				if(sty>0) {
					lattice->trans[st][pos][a] = lattice->trans[sty][sing-shift-1][pattern[sing]];
					lattice->shift[st][pos][a] = lattice->shift[sty][sing-shift-1][pattern[sing]]+shift;
				} else {
					lattice->trans[st][pos][a] = lattice->trans[sty][sing-shift][pattern[sing]];
					lattice->shift[st][pos][a] = lattice->shift[sty][sing-shift][pattern[sing]]+shift;
				}
			}
			dest = st;
			addPos(&dest, pos);
			lattice->trans[st][pos][pattern[pos]] = dest;
			lattice->shift[st][pos][pattern[pos]] = 0;
			for(a=pattern[pos]+1; a<lattice->cardinal; a++) {
				shift = lattice->shift[0][pos][a];
				sty = lattice->trans[0][pos][a];
				if(sty>0) {
					lattice->trans[st][pos][a] = lattice->trans[sty][sing-shift-1][pattern[sing]];
					lattice->shift[st][pos][a] = lattice->shift[sty][sing-shift-1][pattern[sing]]+shift;
				} else {
					lattice->trans[st][pos][a] = lattice->trans[sty][sing-shift][pattern[sing]];
					lattice->shift[st][pos][a] = lattice->shift[sty][sing-shift][pattern[sing]]+shift;
				}
			}
		}
		for(pos=sing+1; pos<lattice->length; pos++) {
			TypeState dest, sty;
			int shift;
			lattice->trans[st][pos-1] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState*));
			lattice->shift[st][pos-1] = (int*) malloc(lattice->cardinal*sizeof(int*));
			for(a=0; a<pattern[pos]; a++) {
				shift = sing-occPrec[sing][pattern[sing]];
				sty = 0;
				if(shift <= sing)
					addPos(&sty, sing-shift);
				if(sty>0) {
					lattice->trans[st][pos-1][a] = lattice->trans[sty][pos-shift-1][a];
					lattice->shift[st][pos-1][a] = lattice->shift[sty][pos-shift-1][a]+shift;
				} else {
					lattice->trans[st][pos-1][a] = lattice->trans[sty][pos-shift][a];
					lattice->shift[st][pos-1][a] = lattice->shift[sty][pos-shift][a]+shift;
				}
			}
			dest = st;
			addPos(&dest, pos);
			lattice->trans[st][pos-1][pattern[pos]] = dest;
			lattice->shift[st][pos-1][pattern[pos]] = 0;
			for(a=pattern[pos]+1; a<lattice->cardinal; a++) {
				shift = sing-occPrec[sing][pattern[sing]];
				sty = 0;
				if(shift <= sing)
					addPos(&sty, sing-shift);
				if(sty>0) {
					lattice->trans[st][pos-1][a] = lattice->trans[sty][pos-shift-1][a];
					lattice->shift[st][pos-1][a] = lattice->shift[sty][pos-shift-1][a]+shift;
				} else {
					lattice->trans[st][pos-1][a] = lattice->trans[sty][pos-shift][a];
					lattice->shift[st][pos-1][a] = lattice->shift[sty][pos-shift][a]+shift;
				}
			}
		}			
	}
	for(pos=0; pos<lattice->length; pos++)
		free((void*)occPrec[pos]);
	free((void*)occPrec);
	set = (int*) malloc(lattice->length*sizeof(int));
	comp = (int*) malloc(lattice->length*sizeof(int));
	for(level=2; level<lattice->length; level++) {
		levelComp = lattice->length-level;
		for(i=0; i<level; i++)
			set[i] = i;
		do {
			TypeState st;
			setToState(set, level, &st);
			stateToSet(~st, lattice->length, comp, &i); //we don't care about i : it is always equal to lattice->length-level, aka levelComp
			lattice->trans[st] = (TypeState**) malloc(levelComp*sizeof(TypeState*));
			lattice->shift[st] = (int**) malloc(levelComp*sizeof(int*));
			for(i=0; i<levelComp; i++) {
				TypeState stx, sty;
				int shift, ind;
				lattice->trans[st][i] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState));
				lattice->shift[st][i] = (int*) malloc(lattice->cardinal*sizeof(int));
				for(a=0; a<pattern[comp[i]]; a++) {
					stx = st;
					remPos(&stx, set[level-1]);
					if(set[level-1]<comp[i]) {
						shift = lattice->shift[stx][i+1][a];
						sty = lattice->trans[stx][i+1][a];
					} else {
						shift = lattice->shift[stx][i][a];
						sty = lattice->trans[stx][i][a];
					}
					if(set[level-1] >= shift) {
						ind = getTransIndex(set[level-1]-shift, sty);
						lattice->trans[st][i][a] = lattice->trans[sty][ind][pattern[set[level-1]]];
						lattice->shift[st][i][a] = lattice->shift[sty][ind][pattern[set[level-1]]]+shift;
					} else {
						lattice->trans[st][i][a] = sty;
						lattice->shift[st][i][a] = shift;
					}
				}
				if(level<(lattice->length-1)) {
					sty = st;
					addPos(&sty, comp[i]);
					lattice->trans[st][i][pattern[comp[i]]] = sty;
					lattice->shift[st][i][pattern[comp[i]]] = 0;
				} else {
					lattice->trans[st][i][pattern[comp[i]]] = stMatch;
					lattice->shift[st][i][pattern[comp[i]]] = shMatch;
				}
				for(a=pattern[comp[i]]+1; a<lattice->cardinal; a++) {
					stx = st;
					remPos(&stx, set[level-1]);
					if(set[level-1]<comp[i]) {
						shift = lattice->shift[stx][i+1][a];
						sty = lattice->trans[stx][i+1][a];
					} else {
						shift = lattice->shift[stx][i][a];
						sty = lattice->trans[stx][i][a];
					}
					if(set[level-1] >= shift) {
						ind = getTransIndex(set[level-1]-shift, sty);
						lattice->trans[st][i][a] = lattice->trans[sty][ind][pattern[set[level-1]]];
						lattice->shift[st][i][a] = lattice->shift[sty][ind][pattern[set[level-1]]]+shift;
					} else {
						lattice->trans[st][i][a] = sty;
						lattice->shift[st][i][a] = shift;
					}
				}
			}
			for(j=level-1; j>=0 && set[j]>=(lattice->length-level+j); j--);
			if(j>=0) {
				set[j]++;
				for(k=j+1; k<level; k++)
					set[k] = set[k-1]+1;
			}
		} while(j>=0);
	}
	free((void*) set);
	free((void*)comp);
	return lattice;
}

/*free lattice*/
void freeLattice(TypeLattice *lattice) {
	int level, levelComp, i, j, k, *set;
	if(lattice != NULL) {
		set = (int*) malloc(lattice->length*sizeof(int));
		for(level=0; level<lattice->length; level++) {
			levelComp = lattice->length-level;
			for(i=0; i<level; i++)
				set[i] = i;
			do {
				TypeState st;
				setToState(set, level, &st);
				for(i=0; i<levelComp; i++) {
					free((void*)lattice->trans[st][i]);
					free((void*)lattice->shift[st][i]);
				}
				free((void*)lattice->trans[st]);
				free((void*)lattice->shift[st]);
				for(j=level-1; j>=0 && set[j]>=(lattice->length-level+j); j--);
				if(j>=0) {
					set[j]++;
					for(k=j+1; k<level; k++)
						set[k] = set[k-1]+1;
				}
			} while(j>=0);
		}
		free((void*)lattice->trans);
		free((void*)lattice->shift);
		free((void*)lattice);
		free((void*)set);
	}
}

/*print the set/state st in text format*/
void fprintState(FILE *f, TypeState st, int length) {
	int i;
	fprintf(f, " ");
	for(i=0; i<length; i++)
		if(isIn(i, st))
			fprintf(f, "%d ", i);
}

/*print set/state st in TeX format*/
void fprintStateTex(FILE *f, TypeState st, int length) {
	int i, l=0;
	for(i=0; i<length; i++)
		if(isIn(i, st))
			l++;
	if(l==0)
		fprintf(f, "$\\emptyset$");
	else {
		fprintf(f, "$\\{");
		for(i=0; !isIn(i, st); i++);
		fprintf(f, "%d", i);
		i++;
		for(; i<length; i++)
			if(isIn(i, st))
				fprintf(f, ", %d", i);
		fprintf(f, "\\}$");
	}
}

/*print lattice in text format*/
void fprintLattice(FILE *f, TypeLattice *lattice, char *alphabet) {
	TypeState size, st;
	int i, *comp, levelComp;
	TypeCardinal a;
	
	size = (1<<lattice->length) - 1;
	comp = (int*) malloc(lattice->length*sizeof(int));
	for(st=0; st<size; st++) {
		fprintState(f, st, lattice->length);
		fprintf(f, "\n");
		stateToSet(~st, lattice->length, comp, &levelComp);
		for(i=0; i<levelComp; i++)
			for(a=0; a<lattice->cardinal; a++) {
				if(alphabet)
					fprintf(f, "\t-%d,%c-> ", comp[i], alphabet[a]);
				else
					fprintf(f, "\t-%d,%d-> ", comp[i], a);
				fprintState(f, lattice->trans[st][i][a], lattice->length);
				fprintf(f, " (shift %d)\n", lattice->shift[st][i][a]);
			}
		fprintf(f, "-------------------------------------------------------------------\n");
	}
	free((void*)comp);
}

/*print lattice in dot format*/
void fprintLatticeDot(FILE *f, TypeLattice *lattice, char *alphabet) {
	TypeState size, st;
	int i, *comp, *set, level, levelComp;
	TypeCardinal a;
	
	size = (1<<lattice->length) - 1;
	set = (int*) malloc(lattice->length*sizeof(int));
	comp = (int*) malloc(lattice->length*sizeof(int));
	fprintf(f, "digraph matching_machine {\nrankdir=TB;");
	for(st=0; st<size; st++) {
		fprintf(f, "node [shape=ellipse, label = \"{");
		fprintState(f, st, lattice->length);
		fprintf(f, "}\"] %lu;\n",  (unsigned long) st);
	}
	for(level=0; level<(lattice->length-1); level++) {
		int j, k;
		fprintf(f, "{ rank = same;");
		for(i=0; i<level; i++)
			set[i] = i;
		do{
			setToState(set, level, &st);
			fprintf(f, " %lu;",  (unsigned long) st);
			for(j=level-1; j>=0 && set[j]>=(lattice->length-level+j); j--);
			if(j>=0) {
				set[j]++;
				for(k=j+1; k<level; k++)
					set[k] = set[k-1]+1;
			}
		} while(j>=0);
		fprintf(f, "}\n");
	}
	for(st=0; st<size; st++) {
		stateToSet(~st, lattice->length, comp, &levelComp);
		for(i=0; i<levelComp; i++)
			for(a=0; a<lattice->cardinal; a++) {
				fprintf(f, "%lu->%lu [ label=\"%d/",  (unsigned long) st,  (unsigned long) lattice->trans[st][i][a], comp[i]);
				if(alphabet)
					fprintf(f, "%c", alphabet[a]);
				else
					fprintf(f, "%d", a);
				fprintf(f, "/%d\"", lattice->shift[st][i][a]);
				fprintf(f, " ]\n");
			}
	}
	free((void*)set);
	free((void*)comp);
	fprintf(f, "}\n");
}

/*return the factorial of n*/
int factorial(int n) {
	int i, res = 1;
	for(i=2; i<=n; i++)
		res *= i;
	return res;
}

#define LATT_SPACE 50

/*print lattice in GasTeX format*/
void fprintLatticeGasTex(FILE *f, TypeLattice *lattice, char *alphabet) {
	TypeState size, st;
	int i, *comp, *set, level, levelComp, maxline, center;
	TypeCardinal a;
	
	maxline = factorial(lattice->length)/(factorial(lattice->length/2)*factorial(lattice->length-lattice->length/2));
	size = (1<<lattice->length) - 1;
	set = (int*) malloc(lattice->length*sizeof(int));
	comp = (int*) malloc(lattice->length*sizeof(int));
	fprintf(f, "\\begin{picture}(%d,%d)(0,0)", (maxline)*LATT_SPACE, (lattice->length)*LATT_SPACE);
	center = ((maxline-1)*LATT_SPACE)/2;
	for(level=0; level<lattice->length; level++) {
		int j, k, size, start;
		size = factorial(lattice->length)/(factorial(lattice->length-level)*factorial(level));
		start = center-((size-1)*LATT_SPACE)/2;
		for(i=0; i<level; i++)
			set[i] = i;
		i=0;
		do{
			setToState(set, level, &st);
			fprintf(f, "\\node[Nframe=n](n%lu)(%d,%d){",  (unsigned long) st, start+i*LATT_SPACE, (lattice->length-level)*LATT_SPACE);
			i++;
			fprintStateTex(f, st, lattice->length);
			fprintf(f, "}\n");
			for(j=level-1; j>=0 && set[j]>=(lattice->length-level+j); j--);
			if(j>=0) {
				set[j]++;
				for(k=j+1; k<level; k++)
					set[k] = set[k-1]+1;
			}
		} while(j>=0);
	}
	for(st=0; st<size; st++) {
		stateToSet(~st, lattice->length, comp, &levelComp);
		for(i=0; i<levelComp; i++)
			for(a=0; a<lattice->cardinal; a++) {
				if(st != lattice->trans[st][i][a])
					fprintf(f, "\\drawedge(n%lu,n%lu)",  (unsigned long)st,  (unsigned long)lattice->trans[st][i][a]);
				else
						fprintf(f, "\\drawloop(n%lu)",  (unsigned long)st);
				fprintf(f, "{$\\scriptstyle %d/", comp[i]);				
				if(alphabet)
					fprintf(f, "%c", alphabet[a]);
				else
					fprintf(f, "%d", a);
				fprintf(f, "/%d$}\n", lattice->shift[st][i][a]);
			}
	}
	free((void*)set);
	free((void*)comp);
	fprintf(f, "\\end{picture}\n");
}

/*print transition in text format*/
void fprintTransition(FILE *f, TypeState st, TypeLattice *lattice) {
	int i, *comp, levelComp;
	TypeCardinal a;
	comp = (int*) malloc(lattice->length*sizeof(int));
	{
		fprintState(f, st, lattice->length);
		fprintf(f, "\n");
		stateToSet(~st, lattice->length, comp, &levelComp);
		for(i=0; i<levelComp; i++)
			for(a=0; a<lattice->cardinal; a++) {
				fprintf(f, "\t-%d,%d-> ", comp[i], a);
				fprintState(f, lattice->trans[st][i][a], lattice->length);
				fprintf(f, " (shift %d)\n", lattice->shift[st][i][a]);
			}
		fprintf(f, "-------------------------------------------------------------------\n");
	}
	free((void*)comp);
}



int *getBorderSpe(TypeSymbol *pattern, int length) {
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

void fprintPartialLattice(FILE *f, TypePartialLattice *lattice, char *alphabet) {
	TypeState st;
	int i;
	TypeCardinal a;
	
	fprintf(f, "\n\n\nPartial lattice\n");
	for(st=0; st<lattice->size; st++) {
		fprintf(f, "%lu\n", (unsigned long)st);
		for(i=0; i<lattice->ntrans[st]; i++) {
			for(a=0; a<lattice->cardinal; a++) {
				if(alphabet)
					fprintf(f, "\t-%d,%c-> ", lattice->next[st][i], alphabet[a]);
				else
					fprintf(f, "\t-%d,%d-> ", lattice->next[st][i], a);
				fprintf(f, "%lu (shift %d)\n", (unsigned long) lattice->trans[st][i][a], lattice->shift[st][i][a]);
			}
		}
		fprintf(f, "-------------------------------------------------------------------\n");
	}
	fprintf(f, "\nPre-match states\n");
	for(st=0; st<lattice->sizeTerm; st++)
		fprintf(f, "%lu (%d)\n", (unsigned long) lattice->term[st], lattice->xterm[st]);
}

/*print partial lattice in dot format*/
void fprintPartialLatticeDot(FILE *f, TypePartialLattice *lattice, char *alphabet) {
	TypeState st;
	int i;
	TypeCardinal a;
	
	fprintf(f, "digraph matching_machine {\nrankdir=TB;");
	for(st=0; st<lattice->size; st++) {
		fprintf(f, "node [shape=ellipse, label = \"{");
		fprintState(f, st, lattice->length);
		fprintf(f, "}\"] %lu;\n",  (unsigned long) st);
	}
	for(st=0; st<lattice->size; st++)
		for(i=0; i<lattice->ntrans[st]; i++)
			for(a=0; a<lattice->cardinal; a++) {
				fprintf(f, "%lu->%lu [ label=\"%d/",  (unsigned long) st,  (unsigned long) lattice->trans[st][i][a], lattice->next[st][i]);
				if(alphabet)
					fprintf(f, "%c", alphabet[a]);
				else
					fprintf(f, "%d", a);
				fprintf(f, "/%d\"", lattice->shift[st][i][a]);
				fprintf(f, " ]\n");
			}
	fprintf(f, "}\n");
}

void freePartialLattice(TypePartialLattice *lattice) {
	if(lattice != NULL) {
		TypeState st;
		for(st=0; st<lattice->size; st++) {
			int i;
			for(i=0; i<lattice->ntrans[st]; i++) {
				free((void*)lattice->trans[st][i]);
				free((void*)lattice->shift[st][i]);
			}
			free((void*)lattice->next[st]);
			free((void*)lattice->trans[st]);
			free((void*)lattice->shift[st]);
		}
		free((void*)lattice->ntrans);
		free((void*)lattice->next);
		free((void*)lattice->trans);
		free((void*)lattice->shift);
		free((void*)lattice->term);
		free((void*)lattice->xterm);
		free((void*)lattice);
	}
}

TypePartialLattice *getPartialLatticePowerSet(TypeCardinal cardinal, TypeSymbol *pattern, int length, int K) {
	TypeCardinal a;
	TypePartialLattice *lattice;
	TypeState stMatch;
	int **occPrec, sing, pos, level, l, *comp, shMatch, *bord;
	TypePowerSet *ps;
	lattice = (TypePartialLattice*) malloc(sizeof(TypePartialLattice));
	lattice->length = length;
	lattice->cardinal = cardinal;
	ps = newKPowerSet(lattice->length, K);
	lattice->size = ps->size;
	lattice->term = (TypeState*) malloc(lattice->length*sizeof(TypeState));
	lattice->xterm = (TypeSymbol*) malloc(lattice->length*sizeof(TypeSymbol));
	lattice->sizeTerm = 0;
	occPrec = (int**) malloc(lattice->length*sizeof(int*));
	for(pos=0; pos<lattice->length; pos++) {
		occPrec[pos] = (int*) malloc(lattice->cardinal*sizeof(int));
	}
	for(a=0; a<lattice->cardinal; a++)
		occPrec[0][a] = -1;
	for(pos=1; pos<lattice->length; pos++) {
		for(a=0; a<pattern[pos-1]; a++)
			occPrec[pos][a] = occPrec[pos-1][a];
		occPrec[pos][pattern[pos-1]] = pos-1;
		for(a=pattern[pos-1]+1; a<lattice->cardinal; a++)
			occPrec[pos][a] = occPrec[pos-1][a];
	}
	bord = getBorderSpe(pattern, lattice->length);
	stMatch = getPrefixIndexKPowerSet(bord[lattice->length], ps);
	shMatch = lattice->length-bord[lattice->length];
	comp = (int*) malloc((lattice->length+1)*sizeof(int));
	lattice->trans = (TypeState***) malloc(lattice->size*sizeof(TypeState**));
	lattice->shift = (int***) malloc(lattice->size*sizeof(int**));
	lattice->ntrans = (int*) malloc(lattice->size*sizeof(int*));
	lattice->next = (int**) malloc(lattice->size*sizeof(int*));
/*empty state*/
	if(lattice->length == 1) {
		lattice->term[lattice->sizeTerm] = 0;
		lattice->xterm[lattice->sizeTerm] = pattern[0];
		lattice->sizeTerm++;
	}
	lattice->ntrans[0] = lattice->length;
	lattice->next[0] = (int*) malloc(lattice->ntrans[0]*sizeof(int));
	lattice->trans[0] = (TypeState**) malloc(lattice->ntrans[0]*sizeof(TypeState*));
	lattice->shift[0] = (int**) malloc(lattice->ntrans[0]*sizeof(int*));
	for(pos=0; pos<lattice->length; pos++) {
		lattice->next[0][pos] = pos;
		lattice->trans[0][pos] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState));
		lattice->shift[0][pos] = (int*) malloc(lattice->cardinal*sizeof(int));
		for(a=0; a<pattern[pos]; a++) {
			if(occPrec[pos][a] >= 0)
				lattice->trans[0][pos][a] = addPosKPowerSet(0L, occPrec[pos][a], ps);
			else
				lattice->trans[0][pos][a] = 0L;
			lattice->shift[0][pos][a] = pos-occPrec[pos][a];
		}
		if(lattice->length>1) {
			lattice->trans[0][pos][pattern[pos]] = addPosKPowerSet(0L, pos, ps);
			lattice->shift[0][pos][pattern[pos]] = 0;
		} else {
			lattice->trans[0][pos][pattern[pos]] = stMatch;
			lattice->shift[0][pos][pattern[pos]] = shMatch;
		}
		for(a=pattern[pos]+1; a<lattice->cardinal; a++) {
			if(occPrec[pos][a] >= 0)
				lattice->trans[0][pos][a] = addPosKPowerSet(0L, occPrec[pos][a], ps);
			else
				lattice->trans[0][pos][a] = 0L;
			lattice->shift[0][pos][a] = pos-occPrec[pos][a];
		}
	}
	for(l=1; l<lattice->length; l++) {
		TypeState st, bd;
		st = getPrefixIndexKPowerSet(l, ps);
		bd = getPrefixIndexKPowerSet(bord[l], ps);
//printf("A st %lu\n", (unsigned long) st);
		if(l == lattice->length-1) {
			lattice->term[lattice->sizeTerm] = st;
			lattice->xterm[lattice->sizeTerm] = pattern[lattice->length-1];
			lattice->sizeTerm++;
		}
		lattice->ntrans[st] = lattice->length-l;
		lattice->next[st] = (int*) malloc(lattice->ntrans[st]*sizeof(int));
		lattice->trans[st] = (TypeState**) malloc(lattice->ntrans[st]*sizeof(TypeState*));
		lattice->shift[st] = (int**) malloc(lattice->ntrans[st]*sizeof(int*));
		for(pos=l; pos<lattice->length; pos++) {
			int index = pos-l;
			lattice->next[st][index] = pos;
			lattice->trans[st][index] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState));
			lattice->shift[st][index] = (int*) malloc(lattice->cardinal*sizeof(int));
			for(a=0; a<pattern[pos]; a++) {
				lattice->trans[st][index][a] = lattice->trans[bd][pos-l][a];
				lattice->shift[st][index][a] = l-bord[l]+lattice->shift[bd][pos-l][a];
			}
			if(l<lattice->length-1) {
				lattice->trans[st][index][pattern[pos]] = addPosKPowerSet(st, pos, ps);
				lattice->shift[st][index][pattern[pos]] = 0;
			} else {
				lattice->trans[st][index][pattern[l]] = stMatch;
				lattice->shift[st][index][pattern[l]] = shMatch;
			}
			for(a=pattern[pos]+1; a<lattice->cardinal; a++) {
				lattice->trans[st][index][a] = lattice->trans[bd][pos-l][a];
				lattice->shift[st][index][a] = l-bord[l]+lattice->shift[bd][pos-l][a];
			}
		}
	}
/*singletons*/
	if(K>1) {
		for(sing=1; sing<lattice->length; sing++) {
			TypeState st = addPosKPowerSet(0L, sing, ps), sty;
			int shift;
//printf("B st %lu\n", (unsigned long) st);
			if(lattice->length == 2) {
				lattice->term[lattice->sizeTerm] = st;
				lattice->xterm[lattice->sizeTerm] = pattern[0];
				lattice->sizeTerm++;
			}
			lattice->ntrans[st] = lattice->length-1;
			lattice->next[st] = (int*) malloc(lattice->ntrans[st]*sizeof(int));
			lattice->trans[st] = (TypeState**) malloc(lattice->ntrans[st]*sizeof(TypeState*));
			lattice->shift[st] = (int**) malloc(lattice->ntrans[st]*sizeof(int*));
			for(pos=0; pos<sing; pos++) {
				lattice->next[st][pos] = pos;
				lattice->trans[st][pos] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState*));
				lattice->shift[st][pos] = (int*) malloc(lattice->cardinal*sizeof(int*));
				for(a=0; a<pattern[pos]; a++) {
					shift = lattice->shift[0][pos][a];
					sty = lattice->trans[0][pos][a];
					if(sty>0) {
						lattice->trans[st][pos][a] = lattice->trans[sty][sing-shift-1][pattern[sing]];
						lattice->shift[st][pos][a] = lattice->shift[sty][sing-shift-1][pattern[sing]]+shift;
					} else {
						lattice->trans[st][pos][a] = lattice->trans[sty][sing-shift][pattern[sing]];
						lattice->shift[st][pos][a] = lattice->shift[sty][sing-shift][pattern[sing]]+shift;
					}
				}
				if(lattice->length>2) {
					lattice->trans[st][pos][pattern[pos]] = addPosKPowerSet(st, pos, ps);
					lattice->shift[st][pos][pattern[pos]] = 0;
				} else {
					lattice->trans[st][pos][pattern[pos]] = stMatch;
					lattice->shift[st][pos][pattern[pos]] = shMatch;
				}
				for(a=pattern[pos]+1; a<lattice->cardinal; a++) {
					shift = lattice->shift[0][pos][a];
					sty = lattice->trans[0][pos][a];
					if(sty>0) {
						lattice->trans[st][pos][a] = lattice->trans[sty][sing-shift-1][pattern[sing]];
						lattice->shift[st][pos][a] = lattice->shift[sty][sing-shift-1][pattern[sing]]+shift;
					} else {
						lattice->trans[st][pos][a] = lattice->trans[sty][sing-shift][pattern[sing]];
						lattice->shift[st][pos][a] = lattice->shift[sty][sing-shift][pattern[sing]]+shift;
					}
				}
			}
			for(pos=sing+1; pos<lattice->length; pos++) {
				lattice->next[st][pos-1] = pos;
				lattice->trans[st][pos-1] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState*));
				lattice->shift[st][pos-1] = (int*) malloc(lattice->cardinal*sizeof(int*));
				for(a=0; a<pattern[pos]; a++) {
					shift = sing-occPrec[sing][pattern[sing]];
					if(shift <= sing)
						sty = addPosKPowerSet(0L, sing-shift, ps);
					else
						sty = 0L;
					if(sty>0) {
						lattice->trans[st][pos-1][a] = lattice->trans[sty][pos-shift-1][a];
						lattice->shift[st][pos-1][a] = lattice->shift[sty][pos-shift-1][a]+shift;
					} else {
						lattice->trans[st][pos-1][a] = lattice->trans[sty][pos-shift][a];
						lattice->shift[st][pos-1][a] = lattice->shift[sty][pos-shift][a]+shift;
					}
				}
				if(lattice->length>2) {
					lattice->trans[st][pos-1][pattern[pos]] = addPosKPowerSet(st, pos, ps);
					lattice->shift[st][pos-1][pattern[pos]] = 0;
				} else {
					lattice->trans[st][pos-1][pattern[pos]] = stMatch;
					lattice->shift[st][pos-1][pattern[pos]] = shMatch;
				}
				for(a=pattern[pos]+1; a<lattice->cardinal; a++) {
					shift = sing-occPrec[sing][pattern[sing]];
					if(shift <= sing)
						sty = addPosKPowerSet(0L, sing-shift, ps);
					else
						sty = 0L;
					if(sty>0) {
						lattice->trans[st][pos-1][a] = lattice->trans[sty][pos-shift-1][a];
						lattice->shift[st][pos-1][a] = lattice->shift[sty][pos-shift-1][a]+shift;
					} else {
						lattice->trans[st][pos-1][a] = lattice->trans[sty][pos-shift][a];
						lattice->shift[st][pos-1][a] = lattice->shift[sty][pos-shift][a]+shift;
					}
				}
			}
		}
		for(l=1; l<lattice->length-1; l++) {
			TypeState st, bd;
			int shift;
			for(sing=l+1; sing<lattice->length; sing++) {
				st = addPosKPowerSet(getPrefixIndexKPowerSet(l, ps), sing, ps);
				bd = lattice->trans[getPrefixIndexKPowerSet(bord[l], ps)][sing-l][pattern[sing]];
//printf("C st %lu\n", (unsigned long) st);
				shift = l-bord[l]+lattice->shift[getPrefixIndexKPowerSet(bord[l], ps)][sing-l][pattern[sing]];
				if(l == lattice->length-2) {
					lattice->term[lattice->sizeTerm] = st;
					lattice->xterm[lattice->sizeTerm] = pattern[lattice->length-2];
					lattice->sizeTerm++;
				}
				lattice->ntrans[st] = lattice->length-l-1;
				lattice->next[st] = (int*) malloc(lattice->ntrans[st]*sizeof(int));
				lattice->trans[st] = (TypeState**) malloc(lattice->ntrans[st]*sizeof(TypeState*));
				lattice->shift[st] = (int**) malloc(lattice->ntrans[st]*sizeof(int*));
				for(pos=l; pos<sing; pos++) {
					int index = pos-l, indexBis = pos-MAX(l, shift);
					lattice->next[st][index] = pos;
					lattice->trans[st][index] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState));
					lattice->shift[st][index] = (int*) malloc(lattice->cardinal*sizeof(int));
					for(a=0; a<pattern[pos]; a++) {
						if(pos>shift) {
							lattice->trans[st][index][a] = lattice->trans[bd][indexBis][a];
							lattice->shift[st][index][a] = lattice->shift[bd][indexBis][a]+shift;
						} else {
							lattice->trans[st][index][a] = bd;
							lattice->shift[st][index][a] = shift;
						}
					}
					if(l<lattice->length-2) {
						lattice->trans[st][index][pattern[pos]] = addPosKPowerSet(st, pos, ps);
						lattice->shift[st][index][pattern[pos]] = 0;
					} else {
						lattice->trans[st][index][pattern[pos]] = stMatch;
						lattice->shift[st][index][pattern[pos]] = shMatch;
					}
					for(a=pattern[pos]+1; a<lattice->cardinal; a++) {
						if(pos>shift) {
							lattice->trans[st][index][a] = lattice->trans[bd][indexBis][a];
							lattice->shift[st][index][a] = lattice->shift[bd][indexBis][a]+shift;
						} else {
							lattice->trans[st][index][a] = bd;
							lattice->shift[st][index][a] = shift;
						}
					}
				}
				for(pos=sing+1; pos<lattice->length; pos++) {
					int index = pos-l-1, indexBis;
					if(shift>sing)
						indexBis = pos-shift;
					else
						indexBis = pos-1-MAX(l, shift);
					lattice->next[st][index] = pos;
					lattice->trans[st][index] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState));
					lattice->shift[st][index] = (int*) malloc(lattice->cardinal*sizeof(int));
					for(a=0; a<pattern[pos]; a++) {
						lattice->trans[st][index][a] = lattice->trans[bd][indexBis][a];
						lattice->shift[st][index][a] = lattice->shift[bd][indexBis][a]+shift;
					}
					if(l<lattice->length-2) {
						lattice->trans[st][index][pattern[pos]] = addPosKPowerSet(st, pos, ps);
						lattice->shift[st][index][pattern[pos]] = 0;
					} else {
						lattice->trans[st][index][pattern[pos]] = stMatch;
						lattice->shift[st][index][pattern[pos]] = shMatch;
					}
					for(a=pattern[pos]+1; a<lattice->cardinal; a++) {
						lattice->trans[st][index][a] = lattice->trans[bd][indexBis][a];
						lattice->shift[st][index][a] = lattice->shift[bd][indexBis][a]+shift;
					}
				}
			}
		}
	}
	for(pos=0; pos<lattice->length; pos++)
		free((void*)occPrec[pos]);
	free((void*)occPrec);
/*subsets with more than 1 element*/
	for(level=2; level<K; level++) {
		TypeState n, st, stx, sty;
		int shift, ind;
		for(l=0; l<lattice->length-level; l++) {
			for(n=ps->begin[l+1][level]; n<END_STATE; n=ps->next[n]) {
				int index;
				st = n+ps->start[l];
				stx = removeLastKPowerSet(st, ps);
				getCompKPowerSet(comp, st, ps);
//printf("D st %lu\n", (unsigned long) st);
				if(l+level == lattice->length-1) {
					lattice->term[lattice->sizeTerm] = st;
					lattice->xterm[lattice->sizeTerm] = pattern[comp[0]];
					lattice->sizeTerm++;
				}
				lattice->ntrans[st] = lattice->length-l-level;
				lattice->next[st] = (int*) malloc(lattice->ntrans[st]*sizeof(int));
				lattice->trans[st] = (TypeState**) malloc(lattice->ntrans[st]*sizeof(TypeState*));
				lattice->shift[st] = (int**) malloc(lattice->ntrans[st]*sizeof(int*));
				for(index=0; index<lattice->ntrans[st]; index++) {
					lattice->next[st][index] = comp[index];
					lattice->trans[st][index] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState));
					lattice->shift[st][index] = (int*) malloc(lattice->cardinal*sizeof(int));
					for(a=0; a<pattern[lattice->next[st][index]]; a++) {
						if(lattice->next[st][index]>ps->dict->node[n].pos) {
							sty = lattice->trans[stx][index+1][a];
							shift = lattice->shift[stx][index+1][a];
						} else {
							sty = lattice->trans[stx][index][a];
							shift = lattice->shift[stx][index][a];
						}
						if(ps->dict->node[n].pos>=shift) {
							ind = getTransIndexPowerSet(ps->dict->node[n].pos-shift, sty, ps);
							lattice->trans[st][index][a] = lattice->trans[sty][ind][pattern[ps->dict->node[n].pos]];
							lattice->shift[st][index][a] = shift+lattice->shift[sty][ind][pattern[ps->dict->node[n].pos]];
						} else {
							lattice->trans[st][index][a] = sty;
							lattice->shift[st][index][a] = shift;
						}
					}
					if(l+level<lattice->length-1) {
						lattice->trans[st][index][pattern[lattice->next[st][index]]] = addPosKPowerSet(st, lattice->next[st][index], ps);
						lattice->shift[st][index][pattern[lattice->next[st][index]]] = 0;
					} else {
						lattice->trans[st][index][pattern[lattice->next[st][index]]] = stMatch;
						lattice->shift[st][index][pattern[lattice->next[st][index]]] = shMatch;
					}
					for(a=pattern[lattice->next[st][index]]+1; a<lattice->cardinal; a++) {
						if(lattice->next[st][index]>ps->dict->node[n].pos) {
							sty = lattice->trans[stx][index+1][a];
							shift = lattice->shift[stx][index+1][a];
						} else {
							sty = lattice->trans[stx][index][a];
							shift = lattice->shift[stx][index][a];
						}
						if(ps->dict->node[n].pos>=shift) {
							ind = getTransIndexPowerSet(ps->dict->node[n].pos-shift, sty, ps);
							lattice->trans[st][index][a] = lattice->trans[sty][ind][pattern[ps->dict->node[n].pos]];
							lattice->shift[st][index][a] = shift+lattice->shift[sty][ind][pattern[ps->dict->node[n].pos]];
						} else {
							lattice->trans[st][index][a] = sty;
							lattice->shift[st][index][a] = shift;
						}
					}
				}
			}
		}
	}
	if(K<lattice->length) {
		TypeState n, st, stx, sty;
		int shift, ind;
		for(l=0; l<lattice->length-K; l++) {
			for(n=ps->begin[l+1][K]; n<END_STATE; n=ps->next[n]) {
				st = n+ps->start[l];
				stx = removeLastKPowerSet(st, ps);
//printf("E st %lu\n", (unsigned long) st);
				if(l+K == lattice->length-1) {
					lattice->term[lattice->sizeTerm] = st;
					lattice->xterm[lattice->sizeTerm] = pattern[l];
					lattice->sizeTerm++;
				}
				lattice->ntrans[st] = 1;
				lattice->next[st] =  (int*) malloc(sizeof(int));;
				lattice->trans[st] = (TypeState**) malloc(sizeof(TypeState*));
				lattice->shift[st] = (int**) malloc(sizeof(int*));
				lattice->next[st][0] = l;
				lattice->trans[st][0] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState));
				lattice->shift[st][0] = (int*) malloc(lattice->cardinal*sizeof(int));
				for(a=0; a<pattern[l]; a++) {
					shift = lattice->shift[stx][0][a];
					sty = lattice->trans[stx][0][a];
					if(ps->dict->node[n].pos >= shift) {
						ind = getTransIndexPowerSet(ps->dict->node[n].pos-shift, sty, ps);
						lattice->trans[st][0][a] = lattice->trans[sty][ind][pattern[ps->dict->node[n].pos]];
						lattice->shift[st][0][a] = lattice->shift[sty][ind][pattern[ps->dict->node[n].pos]]+shift;
					} else {
						lattice->trans[st][0][a] = sty;
						lattice->shift[st][0][a] = shift;
					}
				}
				if(l+K<lattice->length-1) {
					lattice->trans[st][0][pattern[l]] = addPosKPowerSet(st, l, ps);
					lattice->shift[st][0][pattern[l]] = 0;
				} else {
					lattice->trans[st][0][pattern[l]] = stMatch;
					lattice->shift[st][0][pattern[l]] = shMatch;
				}
				for(a=pattern[l]+1; a<lattice->cardinal; a++) {
					shift = lattice->shift[stx][0][a];
					sty = lattice->trans[stx][0][a];
					if(ps->dict->node[n].pos >= shift) {
						ind = getTransIndexPowerSet(ps->dict->node[n].pos-shift, sty, ps);
						lattice->trans[st][0][a] = lattice->trans[sty][ind][pattern[ps->dict->node[n].pos]];
						lattice->shift[st][0][a] = lattice->shift[sty][ind][pattern[ps->dict->node[n].pos]]+shift;
					} else {
						lattice->trans[st][0][a] = sty;
						lattice->shift[st][0][a] = shift;
					}
				}
			}
		}
	}
	lattice->term = (TypeState*) realloc(lattice->term, lattice->sizeTerm*sizeof(TypeState*));
	lattice->xterm = (TypeSymbol*) realloc(lattice->xterm, lattice->sizeTerm*sizeof(TypeSymbol));
	free((void*)bord);
	free((void*)comp);
	freeKPowerSet(ps);
	return lattice;
}

