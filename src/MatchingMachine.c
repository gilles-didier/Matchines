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

#include "MatchingMachine.h"
#include "Utils.h"
#include "PowerSet.h"

#define ST_FREE -1
#define ST_STACKED -2
#define INC_STRATEGY_TMP_SIZE_DICT 1000



typedef struct STRATEGY_TMP_STATE {
	TypeState state;
	int *mem;
} TypeStateTmp;

/*desallocate matchine*/
void freeMatchingMachine(TypeMatchingMachine *matchine) {
	if(matchine != NULL) {
		int i;
		for(i=0; i<matchine->size; i++) {
			free((void*)matchine->trans[i]);
			free((void*)matchine->shift[i]);
		}
		free((void*)matchine->trans);
		free((void*)matchine->shift);
		free((void*)matchine->term);
		free((void*)matchine->xterm);
		free((void*)matchine->next);
		free((void*)matchine);
	}
}

/*print matchine as text*/
void fprintMatchingMachine(FILE *f, TypeMatchingMachine *matchine, char *alphabet) {
	TypeState st;
	TypeCardinal a;
	
	for(st=0; st<matchine->size; st++) {
		fprintf(f, "%lu\n",  (unsigned long) st);
		for(a=0; a<matchine->cardinal; a++) {
			if(alphabet)
				fprintf(f, "\t-%d,%c-> ", matchine->next[st], alphabet[a]);
			else
				fprintf(f, "\t-%d,%d-> ", matchine->next[st], a);
			fprintf(f, "%lu (shift %d)\n",  (unsigned long) matchine->trans[st][a], matchine->shift[st][a]);
		}
		fprintf(f, "-------------------------------------------------------------------\n");
	}
}
 
#define LABEL "label_"

/*print C source implementing matchine*/
void fprintMatchingMachineCode(FILE *f, TypeMatchingMachine *matchine) {
	TypeState *isTerm, st;
	
	isTerm = (TypeState*) malloc(matchine->size*sizeof(TypeState));
	for(st=0; st<matchine->size; st++)
		isTerm[st] = SINK;
	for(st=0; st<matchine->sizeTerm; st++)
		isTerm[matchine->term[st]] = st;
	fprintf(f, "\n\nint foo(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n) {\n\tint p=0, count=0;\n\tTypeSymbol c;\n");
	for(st=0; st<matchine->size; st++) {
		TypeSymbol d;
		fprintf(f, "\n\t%s%lu:\n", LABEL, (unsigned long) st);
		fprintf(f,"\tswitch(y[p+%d]){\n", matchine->next[st]);
		if(isTerm[st] != SINK) {
			for(d=0; d<matchine->xterm[isTerm[st]]; d++)
				fprintf(f,"\t\tcase %u:\n\t\t\tp += %d;\n\t\t\tgoto %s%lu;\n", d, matchine->shift[st][d], LABEL, (unsigned long) matchine->trans[st][d]);
			fprintf(f,"\t\tcase %u:\n\t\t\tp += %d;\n\t\t\tif(p<n) {\n\t\t\t\tcount++;\n\t\t\t\tgoto %s%lu;\n\t\t\t} else\n\t\t\t\treturn count;\n", d, matchine->shift[st][d], LABEL, (unsigned long) matchine->trans[st][d]);
			for(d=matchine->xterm[isTerm[st]]+1; d<matchine->cardinal; d++)
				fprintf(f,"\t\tcase %u:\n\t\t\tp += %d;\n\t\t\tgoto %s%lu;\n", d, matchine->shift[st][d], LABEL, (unsigned long) matchine->trans[st][d]);
		} else {
			for(d=0; d<matchine->cardinal; d++)
				fprintf(f,"\t\tcase %u:\n\t\t\tp += %d;\n\t\t\tgoto %s%lu;\n", d, matchine->shift[st][d], LABEL, (unsigned long) matchine->trans[st][d]);
		}
		fprintf(f,"\t}\n");
	}
	fprintf(f,"\treturn count;\n}\n");
	free((void*)isTerm);
}

/*print matchine in .gv format*/
void fprintMatchingMachineDot(FILE *f, TypeMatchingMachine *matchine, char *alphabet) {
	TypeState st;
	TypeCardinal a;
	TypeState *isTerm;
	
	isTerm = (TypeState*) malloc(matchine->size*sizeof(TypeState));
	for(st=0; st<matchine->size; st++)
		isTerm[st] = SINK;
	for(st=0; st<matchine->sizeTerm; st++)
		isTerm[matchine->term[st]] = st;
	fprintf(f, "digraph matching_machine {\n");
	if(isTerm[0] != SINK)
		fprintf(f, "node [shape=diamond, style=filled, color=lightgrey, label = \"state %lu (%d)\"] %lu;\n",  0L, matchine->next[0],  0L);
	else
		fprintf(f, "node [shape=diamond, style=solid, color=black, label = \"state %lu (%d)\"] %lu;\n",  0L, matchine->next[0],  0L);
	for(st=1; st<matchine->size; st++)
		if(isTerm[st] != SINK)
			fprintf(f, "node [shape=ellipse, style=filled, color=lightgrey, label = \"state %lu (%d)\"] %lu;\n",  (unsigned long) st, matchine->next[st],  (unsigned long) st);
		else
			fprintf(f, "node [shape=ellipse, style=solid, color=black, label = \"state %lu (%d)\"] %lu;\n",  (unsigned long) st, matchine->next[st],  (unsigned long) st);
	for(st=0; st<matchine->size; st++) {
		if(isTerm[st] != SINK) {
			for(a=0; a<matchine->cardinal; a++) {
				if(matchine->trans[st][a] != SINK) {
					fprintf(f, "%lu->%lu [ label=\"",  (unsigned long) st,  (unsigned long) matchine->trans[st][a]);
					if(alphabet)
						fprintf(f, "%c", alphabet[a]);
					else
						fprintf(f, "%d", a);
					fprintf(f, "/%d\"", matchine->shift[st][a]);
					if(a == matchine->xterm[isTerm[st]])
						fprintf(f, ", style=bold, color=blue");
//						fprintf(f, ",penwidth = 3");
					fprintf(f, " ]\n");
				}
			}
		} else {
			for(a=0; a<matchine->cardinal; a++) {
				if(matchine->trans[st][a] != SINK) {
					fprintf(f, "%lu->%lu [ label=\"",  (unsigned long) st,  (unsigned long) matchine->trans[st][a]);
					if(alphabet)
						fprintf(f, "%c", alphabet[a]);
					else
						fprintf(f, "%d", a);
					fprintf(f, "/%d\" ]\n", matchine->shift[st][a]);
				}
			}
		}
	}
	fprintf(f, "}\n");
	free((void*) isTerm);
}

/*get bounds of the matchine, in particular its order*/
void getMatchingMachineBounds(TypeMatchingMachine *matchine, int *minNext, int *maxNext, int *minShift, int *maxShift) {
	TypeState st;
	TypeCardinal a;
	*minNext = matchine->next[0];
	*maxNext = matchine->next[0];
	*minShift = matchine->shift[0][0];
	*maxShift = matchine->shift[0][0];
	for(a=1; a<matchine->cardinal; a++) {
		if(matchine->shift[0][a]<*minShift)
			*minShift = matchine->shift[0][a];
		if(matchine->shift[0][a]>*maxShift)
			*maxShift = matchine->shift[0][a];
	}
	for(st=1; st<matchine->size; st++) {
		if(matchine->next[st]<*minNext)
			*minNext = matchine->next[st];
		if(matchine->next[st]>*maxNext)
			*maxNext = matchine->next[st];
		for(a=0; a<matchine->cardinal; a++) {
			if(matchine->shift[st][a]<*minShift)
				*minShift = matchine->shift[st][a];
			if(matchine->shift[st][a]>*maxShift)
				*maxShift = matchine->shift[st][a];
		}
	}
}
		
/*return the full memory extension of matchine*/
TypeMatchingMachine *expandMatchingMachine(TypeMatchingMachine *matchine) {
	TypeMatchingMachine *res;
	TypeState st, sizeStack, sizeBufferState, sizeBufferStack, incBuffer;
	TypeLexiTreeSet **lexi;
	TypeStateTmp *stack;
	int i, minNext, maxNext, minShift, maxShift, length;

	getMatchingMachineBounds(matchine, &minNext, &maxNext, &minShift, &maxShift);
	length = maxNext+1;
	if(minShift<0) {
		printf("Warning negative shift (%d)\n", minShift);
	}
	if(minNext<0) {
		printf("Warning negative next (%d)\n", minNext);
	}
	res = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	res->cardinal = matchine->cardinal;
	res->size = 0;
	incBuffer = matchine->size;
	sizeBufferState = incBuffer;
	res->next = (int*) malloc(sizeBufferState*sizeof(int));
	res->shift = (int**) malloc(sizeBufferState*sizeof(int*));
	res->trans = (TypeState**) malloc(sizeBufferState*sizeof(TypeState*));
	res->sizeTerm = matchine->sizeTerm;
	lexi = (TypeLexiTreeSet**) malloc(matchine->size*sizeof(TypeLexiTreeSet*));
	for(st=0; st<matchine->size; st++)
		lexi[st] = newLexiTreeSet(length);
	sizeBufferStack = incBuffer;
	stack = (TypeStateTmp*) malloc(sizeBufferStack*sizeof(TypeStateTmp));
	stack[0].state = 0;
	stack[0].mem = (int*) malloc((length)*sizeof(int));
	for(i=0; i<length; i++)
		stack[0].mem[i] = matchine->cardinal;
	addSetLexiTreeSet(stack[0].mem, res->size, lexi[stack[0].state]);
	res->size++;
	sizeStack = 1;
	while(sizeStack>0) {
		TypeStateTmp cur;
		TypeState st;
		TypeCardinal a;
		cur = stack[--sizeStack];
		st = findSetLexiTreeSet(cur.mem, lexi[cur.state]);
		res->next[st] = matchine->next[cur.state];
		res->trans[st] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		res->shift[st] = (int*) malloc(matchine->cardinal*sizeof(int));
		if(cur.mem[matchine->next[cur.state]] == matchine->cardinal) {
			for(a=0; a<matchine->cardinal; a++) {
				TypeState fol;
				int *mem;
				int i, j;
				res->shift[st][a] = matchine->shift[cur.state][a];
				cur.mem[matchine->next[cur.state]] = a;
				mem = (int*) malloc((length)*sizeof(int));
				for(i=res->shift[st][a],j=0; i<length; i++,j++)
					mem[j] = cur.mem[i];
				for(;j<length;j++)
					mem[j] = matchine->cardinal;
				fol = findSetLexiTreeSet(mem, lexi[matchine->trans[cur.state][a]]);
				if(fol == SINK) {
					if(sizeStack>=sizeBufferStack) {
						sizeBufferStack += incBuffer;
						stack = (TypeStateTmp*) realloc((void*) stack, sizeBufferStack*sizeof(TypeStateTmp));	
					}
					if(res->size>=sizeBufferState) {
						sizeBufferState += incBuffer;
						res->next = (int*) realloc((void*) res->next, sizeBufferState*sizeof(int));
						res->shift = (int**) realloc((void*) res->shift, sizeBufferState*sizeof(int*));
						res->trans = (TypeState**) realloc((void*) res->trans, sizeBufferState*sizeof(TypeState*));
					}
					addSetLexiTreeSet(mem, res->size, lexi[matchine->trans[cur.state][a]]);
					fol = res->size;
					res->size++;
					stack[sizeStack].state = matchine->trans[cur.state][a];
					stack[sizeStack].mem = mem;
					sizeStack++;
				} else
					free((void*)mem);
				res->trans[st][a] = fol;
			}
		} else {
			TypeState fol;
			int *mem;
			int i, j, next;
			for(a=0; a<matchine->cardinal; a++) {
				res->trans[st][a] = SINK;
				res->shift[st][a] = 0;
			}
			next = matchine->next[cur.state];
			res->shift[st][cur.mem[next]] = matchine->shift[cur.state][cur.mem[next]];
			mem = (int*) malloc((length)*sizeof(int));
			for(i=res->shift[st][cur.mem[next]],j=0; i<length && j<length; i++,j++)
				mem[j] = cur.mem[i];
			for(;j<length;j++)
				mem[j] = matchine->cardinal;
			fol = findSetLexiTreeSet(mem, lexi[matchine->trans[cur.state][cur.mem[next]]]);
			if(fol == SINK) {
				if(sizeStack>=sizeBufferStack) {
					sizeBufferStack += incBuffer;
					stack = (TypeStateTmp*) realloc((void*) stack, sizeBufferStack*sizeof(TypeStateTmp));	
				}
				if(res->size>=sizeBufferState) {
					sizeBufferState += incBuffer;
					res->next = (int*) realloc((void*) res->next, sizeBufferState*sizeof(int));
					res->shift = (int**) realloc((void*) res->shift, sizeBufferState*sizeof(int*));
					res->trans = (TypeState**) realloc((void*) res->trans, sizeBufferState*sizeof(TypeState*));
				}
				addSetLexiTreeSet(mem, res->size, lexi[matchine->trans[cur.state][cur.mem[next]]]);
				fol = res->size;
				res->size++;
				stack[sizeStack].state = matchine->trans[cur.state][cur.mem[next]];
				stack[sizeStack].mem = mem;
				sizeStack++;
			} else
				free((void*)mem);
			res->trans[st][cur.mem[next]] = fol;
		}
		free((void*)cur.mem);
	}
	res->next = (int*) realloc((void*) res->next, res->size*sizeof(int));
	res->shift = (int**) realloc((void*) res->shift, res->size*sizeof(int*));
	res->trans = (TypeState**) realloc((void*) res->trans, res->size*sizeof(TypeState*));
	res->sizeTerm = matchine->sizeTerm;
	res->sizeTerm = 0;
	for(st=0; st<matchine->sizeTerm; st++)
		res->sizeTerm += lexi[matchine->term[st]]->size;
	res->term = (TypeState*) malloc(res->sizeTerm*sizeof(TypeState));
	res->xterm = (TypeSymbol*) malloc(res->sizeTerm*sizeof(TypeSymbol));
	res->sizeTerm = 0;
	for(st=0; st<matchine->sizeTerm; st++) {
		TypeState *tmp, size, x;
		getTermIndex(lexi[matchine->term[st]], &tmp, &size);
		for(x=0; x<size; x++) {
			res->term[res->sizeTerm] = tmp[x];
			res->xterm[res->sizeTerm] = matchine->xterm[st];
			res->sizeTerm++;
		}
		free((void*)tmp);
	}
	res->term = (TypeState*) realloc(res->term, res->sizeTerm*sizeof(TypeState));
	res->xterm = (TypeSymbol*) realloc(res->xterm, res->sizeTerm*sizeof(TypeSymbol));
	for(st=0; st<matchine->size; st++)
		freeLexiTreeSet(lexi[st]);
	free((void*)lexi);
	free((void*)stack);
	return res;
}

/*simplify matchine*/
TypeMatchingMachine *simplify(TypeMatchingMachine *matchine) {
	TypeMatchingMachine *res;
	TypeState *stack, *index, st;
	int last=-1, *isTerm;
	stack = (TypeState*) malloc(matchine->size*sizeof(TypeState));
	index = (TypeState*) malloc(matchine->size*sizeof(TypeState));
	for(st=0; st<matchine->size; st++)
		index[st] = NO_STATE;
	res = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	res->size = 0;
	res->next = (int*) malloc(matchine->size*sizeof(int));
	res->trans = (TypeState**) malloc(matchine->size*sizeof(TypeState*));
	res->shift = (int**) malloc(matchine->size*sizeof(int*));
	res->cardinal = matchine->cardinal;
	res->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	res->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	isTerm = (int*) malloc(matchine->size*sizeof(int));
	for(st=0; st<matchine->size; st++)
		isTerm[st] = 0;
	for(st=0; st<matchine->sizeTerm; st++)
		isTerm[matchine->term[st]] = 1+matchine->xterm[st];
	res->sizeTerm = 0;
	stack[++last] = 0;
	index[0] = res->size++;
	while(last>=0) {
		TypeCardinal a;
		TypeState cur = stack[last--];
		res->next[index[cur]] = matchine->next[cur];
		res->trans[index[cur]] = (TypeState*) malloc(matchine->cardinal*sizeof(TypeState));
		res->shift[index[cur]] = (int*) malloc(matchine->cardinal*sizeof(int));
		for(a=0; a<matchine->cardinal; a++) {
			if(index[matchine->trans[cur][a]] == NO_STATE) {
				index[matchine->trans[cur][a]] = res->size++;
				if(isTerm[matchine->trans[cur][a]]) {
					res->term[res->sizeTerm] = index[matchine->trans[cur][a]];	
					res->xterm[res->sizeTerm] = isTerm[matchine->trans[cur][a]]-1;
					res->sizeTerm++;
				}
				stack[++last] = matchine->trans[cur][a];
			}
			res->trans[index[cur]][a] = index[matchine->trans[cur][a]];
			res->shift[index[cur]][a] = matchine->shift[cur][a];
		}
	}
	free((void*)isTerm);
	res->next = (int*) realloc((void*)res->next, res->size*sizeof(int));
	res->trans = (TypeState**) realloc((void*)res->trans, res->size*sizeof(TypeState*));
	res->shift = (int**) realloc((void*)res->shift, res->size*sizeof(int*));
	res->term = (TypeState*) realloc((void*)res->term, res->sizeTerm*sizeof(TypeState));
	res->xterm = (TypeSymbol*) realloc((void*)res->xterm, res->sizeTerm*sizeof(TypeSymbol));
	free((void*)stack);
	free((void*)index);
	return res;
}


TypeMatchingMachine *getNewKMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length, TypeBernoulli *bernoulli, int K) {
	TypeMatchingMachine *matchine;
	int *index, *remain, *next, levelComp, i, k, *comp, is;
	TypeCardinal a;
	TypeLattice *lattice = getLattice(cardinal, pattern, length);
	TypeState *status, st, *stack, size;
	double *expShiftA, *expShiftB, *tmp;
	
	size = (1<<lattice->length) - 1;
	expShiftA = (double*) malloc(size*sizeof(double));
	expShiftB = (double*) malloc(size*sizeof(double));
	index = (int*) malloc(size*sizeof(int));
	remain = (int*) malloc(size*sizeof(int));
	next = (int*) malloc(size*sizeof(int));
	comp = (int*) malloc(lattice->length*sizeof(int));
	for(st=0; st<size; st++) {
		stateToSet(~st, lattice->length, comp, &(remain[st]));
		expShiftA[st] = 0.;
		index[st] = 0;
		next[st] = comp[0];
	}
	for(k=0; k<K; k++) {
		for(st=0; st<size; st++) {
			double tmp;
			stateToSet(~st, lattice->length, comp, &levelComp);
			tmp = 0.;
			for(a=0; a<lattice->cardinal; a++)
				tmp += bernoulli->prob[a]*(lattice->shift[st][0][a]+expShiftA[lattice->trans[st][0][a]]);
			expShiftB[st] = tmp;
			index[st] = 0;
			next[st] = comp[0];
			for(i=1; i<levelComp; i++) {
				tmp = 0.;
				for(a=0; a<lattice->cardinal; a++)
					tmp += bernoulli->prob[a]*(lattice->shift[st][i][a]+expShiftA[lattice->trans[st][i][a]]);
				if(tmp>expShiftB[st]) {
					expShiftB[st] = tmp;
					index[st] = i;
					next[st] = comp[i];
				}
			}
		}
		tmp = expShiftA;
		expShiftA = expShiftB;
		expShiftB = tmp;
	}
	free((void*) comp);
	free((void*) expShiftA);
	free((void*) expShiftB);
	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	matchine->cardinal = cardinal;
	matchine->size = 0;
	matchine->trans = (TypeState**) malloc(size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(size*sizeof(int*));
	matchine->next = (int*) malloc(size*sizeof(int));
	stack = (TypeState*) malloc(size*sizeof(TypeState));
	status = (TypeState*) malloc(size*sizeof(TypeState));
	for(st=0; st<size; st++)
		status[st] = ST_FREE;
	is = 0;
	stack[is] = 0;
	status[0] = matchine->size++;
	is++;
	while(is>0) {
		st = stack[--is];
		matchine->trans[status[st]] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState));
		matchine->shift[status[st]] = (int*) malloc(lattice->cardinal*sizeof(int));
		matchine->next[status[st]] = next[st];
		for(a=0; a<lattice->cardinal; a++) {
			if(status[lattice->trans[st][index[st]][a]] == ST_FREE) {
				status[lattice->trans[st][index[st]][a]] = matchine->size++;
				stack[is++] = lattice->trans[st][index[st]][a];
			}
			matchine->trans[status[st]][a] = status[lattice->trans[st][index[st]][a]];
			matchine->shift[status[st]][a] = lattice->shift[st][index[st]][a];
		}
	}
	matchine->sizeTerm = 0;
	for(st=0; st<size; st++)
		if(status[st] != ST_FREE && remain[st] == 1)
			matchine->sizeTerm++;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	matchine->sizeTerm = 0;
	for(st=0; st<size; st++)
		if(status[st] != ST_FREE && remain[st] == 1) {
			matchine->term[matchine->sizeTerm] = status[st];
			matchine->xterm[matchine->sizeTerm] = pattern[matchine->next[status[st]]];
			matchine->sizeTerm++;
		}
	free((void*) status);
	free((void*) stack);
	free((void*) index);
	free((void*) remain);
	free((void*) next);
	freeLattice(lattice);
	matchine->trans = (TypeState**) realloc(matchine->trans, matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) realloc(matchine->shift, matchine->size*sizeof(int*));
	return matchine;
}

TypeMatchingMachine *getReverseMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length) {
	TypeMatchingMachine *matchine;
	int *status, *index, *remain, *next, level, levelComp, i, j, k, *set, *comp, is;
	TypeCardinal a;
	TypeLattice *lattice = getLattice(cardinal, pattern, length);
	TypeState st, *stack, size;
	
	size = (1<<lattice->length) - 1;
	remain = (int*) malloc(size*sizeof(int));
	index = (int*) malloc(size*sizeof(int));
	next = (int*) malloc(size*sizeof(int));
	set = (int*) malloc(lattice->length*sizeof(int));
	comp = (int*) malloc(lattice->length*sizeof(int));
	for(level=lattice->length-1; level>=0; level--) {
		levelComp = lattice->length-level;
		for(i=0; i<level; i++)
			set[i] = i;
		do {
			TypeState st;
			setToState(set, level, &st);
			stateToSet(~st, lattice->length, comp, &i); //we don't care about i
			index[st] = levelComp-1;
			remain[st] = i;
			next[st] = comp[levelComp-1];
			for(j=level-1; j>=0 && set[j]>=(lattice->length-level+j); j--);
			if(j>=0) {
				set[j]++;
				for(k=j+1; k<level; k++)
					set[k] = set[k-1]+1;
			}
		} while(j>=0);
	}
	free((void*) set);
	free((void*) comp);
	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	matchine->cardinal = cardinal;
	matchine->size = 0;
	matchine->trans = (TypeState**) malloc(size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(size*sizeof(int*));
	matchine->next = (int*) malloc(size*sizeof(int));
	stack = (TypeState*) malloc(size*sizeof(TypeState));
	status = (int*) malloc(size*sizeof(int));
	for(st=0; st<size; st++)
		status[st] = ST_FREE;
	is = 0;
	stack[is] = 0;
	status[0] = matchine->size++;
	is++;
	while(is>0) {
		st = stack[--is];
		matchine->trans[status[st]] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState));
		matchine->shift[status[st]] = (int*) malloc(lattice->cardinal*sizeof(int));
		matchine->next[status[st]] = next[st];
		for(a=0; a<lattice->cardinal; a++) {
			if(status[lattice->trans[st][index[st]][a]] == ST_FREE) {
				status[lattice->trans[st][index[st]][a]] = matchine->size++;
				stack[is++] = lattice->trans[st][index[st]][a];
			}
			matchine->trans[status[st]][a] = status[lattice->trans[st][index[st]][a]];
			matchine->shift[status[st]][a] = lattice->shift[st][index[st]][a];
		}
	}
	matchine->sizeTerm = 0;
	for(st=0; st<size; st++)
		if(status[st] != ST_FREE && remain[st] == 1)
			matchine->sizeTerm++;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	matchine->sizeTerm = 0;
	for(st=0; st<size; st++)
		if(status[st] != ST_FREE && remain[st] == 1) {
			matchine->term[matchine->sizeTerm] = status[st];
			matchine->xterm[matchine->sizeTerm] = pattern[matchine->next[status[st]]];
			matchine->sizeTerm++;
		}
	free((void*) remain);
	free((void*) status);
	free((void*) stack);
	free((void*) index);
	free((void*) next);
	freeLattice(lattice);
	matchine->trans = (TypeState**) realloc(matchine->trans, matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) realloc(matchine->shift, matchine->size*sizeof(int*));
	return matchine;
}



TypeMatchingMachine *getPolyKMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length, TypeBernoulli *bernoulli, int K) {
	TypeMatchingMachine *matchine;
	TypePartialLattice *lattice;
	if(K<1)
		K = 1;
	lattice = getPartialLatticePowerSet(cardinal, pattern, length, K);
	matchine = getBestKFromPartial(lattice, bernoulli, K+10);
	freePartialLattice(lattice);
	return matchine;
}



TypeMatchingMachine *getBestKFromPartial(TypePartialLattice *lattice, TypeBernoulli *bernoulli, int K) {
	TypeMatchingMachine *matchine;
	int *index, *next, i, k, is;
	TypeCardinal a;
	TypeState *status, st, *stack;
	double *expShiftA, *expShiftB, *tmp;
	
	matchine = (TypeMatchingMachine*) malloc(sizeof(TypeMatchingMachine));
	expShiftA = (double*) malloc(lattice->size*sizeof(double));
	expShiftB = (double*) malloc(lattice->size*sizeof(double));
	index = (int*) malloc(lattice->size*sizeof(int));
	next = (int*) malloc(lattice->size*sizeof(int));	
	for(st=0; st<lattice->size; st++) {
		expShiftA[st] = 0.;
		index[st] = 0;
		next[st] = lattice->next[st][0];
	}
	for(k=0; k<K; k++) {
		for(st=0; st<lattice->size; st++) {
			double tmp;
			tmp = 0.;
			for(a=0; a<lattice->cardinal; a++) {
				tmp += bernoulli->prob[a]*(lattice->shift[st][0][a]+expShiftA[lattice->trans[st][0][a]]);
			}
			expShiftB[st] = tmp;
			index[st] = 0;
			next[st] = lattice->next[st][0];
			for(i=1; i<lattice->ntrans[st]; i++) {
				tmp = 0.;
				for(a=0; a<lattice->cardinal; a++)
					tmp += bernoulli->prob[a]*(lattice->shift[st][i][a]+expShiftA[lattice->trans[st][i][a]]);
				if(tmp>=expShiftB[st]) {
					expShiftB[st] = tmp;
					index[st] = i;
					next[st] = lattice->next[st][i];
				}
			}
		}
		tmp = expShiftA;
		expShiftA = expShiftB;
		expShiftB = tmp;
	}
	free((void*) expShiftA);
	free((void*) expShiftB);
	matchine->cardinal = lattice->cardinal;
	matchine->size = 0;
	matchine->trans = (TypeState**) malloc(lattice->size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(lattice->size*sizeof(int*));
	matchine->next = (int*) malloc(lattice->size*sizeof(int));
	stack = (TypeState*) malloc(lattice->size*sizeof(TypeState));
	status = (TypeState*) malloc(lattice->size*sizeof(TypeState));
	for(st=0; st<lattice->size; st++)
		status[st] = ST_FREE;
	is = 0;
	stack[is] = 0;
	status[0] = matchine->size++;
	is++;
	while(is>0) {
		st = stack[--is];
		matchine->trans[status[st]] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState));
		matchine->shift[status[st]] = (int*) malloc(lattice->cardinal*sizeof(int));
		matchine->next[status[st]] = next[st];
		for(a=0; a<lattice->cardinal; a++) {
			if(status[lattice->trans[st][index[st]][a]] == ST_FREE) {
				status[lattice->trans[st][index[st]][a]] = matchine->size++;
				stack[is++] = lattice->trans[st][index[st]][a];
			}
			matchine->trans[status[st]][a] = status[lattice->trans[st][index[st]][a]];
			matchine->shift[status[st]][a] = lattice->shift[st][index[st]][a];
		}
	}
	matchine->sizeTerm = 0;
	matchine->term = (TypeState*) malloc(lattice->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(lattice->sizeTerm*sizeof(TypeSymbol));
	matchine->sizeTerm = 0;
	for(st=0; st<lattice->sizeTerm; st++)
		if(status[lattice->term[st]] != ST_FREE) {
			matchine->term[matchine->sizeTerm] = status[lattice->term[st]];
			matchine->xterm[matchine->sizeTerm] = lattice->xterm[st];
			matchine->sizeTerm++;
		}
	matchine->term = (TypeState*) realloc(matchine->term, matchine->sizeTerm*sizeof(TypeState*));
	matchine->xterm = (TypeSymbol*) realloc(matchine->xterm, matchine->sizeTerm*sizeof(TypeSymbol));
	free((void*) status);
	free((void*) stack);
	free((void*) index);
	free((void*) next);
	matchine->trans = (TypeState**) realloc(matchine->trans, matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) realloc(matchine->shift, matchine->size*sizeof(int*));
	return matchine;
}
