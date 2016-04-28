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

#include "Explore.h"
#include "AsymptoticSpeed.h"

/*return the next matching machine from current*/
TypeMatchingMachine *getMatchingMachineFromLatticeNext(TypeLattice *lattice, TypeSymbol *pattern, int *current) {
	TypeMatchingMachine *matchine, *simplerMatchingMachine;
	TypeState st, size;
	int *comp;
	
	comp = (int*) malloc(lattice->length*sizeof(int));
	size = (1<<lattice->length)-1;
	matchine = (TypeMatchingMachine*)  malloc(sizeof(TypeMatchingMachine));
	matchine->cardinal = lattice->cardinal;
	matchine->size = size;
	matchine->trans = (TypeState**) malloc(size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(size*sizeof(int*));
	matchine->next = (int*) malloc(size*sizeof(int));
	matchine->term = (TypeState*) malloc(size*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(size*sizeof(TypeSymbol));
	matchine->sizeTerm = 0;
	for(st=0; st<size; st++) {
		int i, a;
		stateToSet(~st, lattice->length, comp, &i);
		matchine->next[st] = comp[current[st]];
		matchine->trans[st] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState));
		matchine->shift[st] = (int*) malloc(lattice->cardinal*sizeof(int));
		if(i==1) {
			matchine->term[matchine->sizeTerm] = st;
			matchine->xterm[matchine->sizeTerm] = pattern[matchine->next[matchine->term[matchine->sizeTerm]]];
			matchine->sizeTerm++;
		}
		for(a=0; a<lattice->cardinal; a++) {
			matchine->trans[st][a] = lattice->trans[st][current[st]][a];
			matchine->shift[st][a] = lattice->shift[st][current[st]][a];
		}
	}
	free((void*)comp);
	simplerMatchingMachine = simplify(matchine);
	freeMatchingMachine(matchine);
	return simplerMatchingMachine;
}

/*return the fastest matching machine wrt bernoulli model*/
TypeMatchingMachine *explore(TypeCardinal cardinal, TypeSymbol *pattern, int length, TypeBernoulli *bernoulli) {
	TypeMatchingMachine *matchine;
	int *current, *max;
	TypeLattice *lattice = getLattice(cardinal, pattern, length);
	TypeState st, size;
	double maxSpeed = -1.;
	
	size = (1<<lattice->length) - 1;
	current = (int*) malloc(size*sizeof(int));
	max = (int*) malloc(size*sizeof(int));
	for(st=0; st<size; st++) {
		max[st] = lattice->length-stateCardinal(st, lattice->length);
		current[st] = 0;
	}
	matchine = getMatchingMachineFromLatticeNext(lattice, pattern, current);
	maxSpeed = getAsymptoticSpeed(matchine, bernoulli);
	current[0]++;
	do {
		int cont;
		TypeMatchingMachine *matchineTmp = getMatchingMachineFromLatticeNext(lattice, pattern, current);
		double tmpSpeed;
/*compute Speed*/
		tmpSpeed = getAsymptoticSpeed(matchineTmp, bernoulli);
		if(tmpSpeed > maxSpeed) {
			freeMatchingMachine(matchine);
			matchine = matchineTmp;
			maxSpeed = tmpSpeed;
		} else
			freeMatchingMachine(matchineTmp);
/*build the next table for the next matchine*/
		st = 0;
		do {
			cont = 0;
			current[st]++;
			if(current[st] >= max[st]) {
				current[st] = 0;
				cont = 1;
				st++;
			}
		} while(cont && st<size);
	} while(st<size);
	freeLattice(lattice);
	free((void*)current);
	free((void*)max);
	return matchine;
}


/*return the next matching machine from current*/
TypeMatchingMachine *getMatchingMachineFromPartialLatticeNext(TypePartialLattice *lattice, int *current) {
	TypeMatchingMachine *matchine, *simplerMatchingMachine;
	TypeState st;
	
	matchine = (TypeMatchingMachine*)  malloc(sizeof(TypeMatchingMachine));
	matchine->cardinal = lattice->cardinal;
	matchine->size = lattice->size;
	matchine->trans = (TypeState**) malloc(matchine->size*sizeof(TypeState*));
	matchine->shift = (int**) malloc(matchine->size*sizeof(int*));
	matchine->next = (int*) malloc(matchine->size*sizeof(int));
	matchine->sizeTerm = lattice->sizeTerm;
	matchine->term = (TypeState*) malloc(matchine->sizeTerm*sizeof(TypeState));
	matchine->xterm = (TypeSymbol*) malloc(matchine->sizeTerm*sizeof(TypeSymbol));
	for(st=0; st<lattice->size; st++) {
		TypeCardinal a;
		matchine->next[st] = lattice->next[st][current[st]];
		matchine->trans[st] = (TypeState*) malloc(lattice->cardinal*sizeof(TypeState));
		matchine->shift[st] = (int*) malloc(lattice->cardinal*sizeof(int));
		for(a=0; a<lattice->cardinal; a++) {
			matchine->trans[st][a] = lattice->trans[st][current[st]][a];
			matchine->shift[st][a] = lattice->shift[st][current[st]][a];
		}
	}
	for(st=0; st<matchine->sizeTerm; st++) {
		matchine->term[st] = lattice->term[st];
		matchine->xterm[st] = lattice->xterm[st];
	}
	simplerMatchingMachine = simplify(matchine);
	freeMatchingMachine(matchine);
	return simplerMatchingMachine;
}

/*return the fastest matching machine from the partial lattice wrt bernoulli model*/
TypeMatchingMachine *explorePartialLattice(TypePartialLattice *lattice, TypeBernoulli *bernoulli) {
	TypeMatchingMachine *matchine;
	int *current, *max;
	TypeState st;
	double maxSpeed = -1.;
	
	current = (int*) malloc(lattice->size*sizeof(int));
	max = (int*) malloc(lattice->size*sizeof(int));
	for(st=0; st<lattice->size; st++) {
		max[st] = lattice->ntrans[st];
		current[st] = 0;
	}
	matchine = getMatchingMachineFromPartialLatticeNext(lattice, current);
	maxSpeed = getAsymptoticSpeed(matchine, bernoulli);
	current[0]++;
	do {
		int cont;
		TypeMatchingMachine *matchineTmp = getMatchingMachineFromPartialLatticeNext(lattice, current);
		double tmpSpeed;
/*compute Speed*/
		tmpSpeed = getAsymptoticSpeed(matchineTmp, bernoulli);
		if(tmpSpeed > maxSpeed) {
			freeMatchingMachine(matchine);
			matchine = matchineTmp;
			maxSpeed = tmpSpeed;
		} else
			freeMatchingMachine(matchineTmp);
/*build the next table for the next matchine*/
		st = 0;
		do {
			cont = 0;
			current[st]++;
			if(current[st] >= max[st]) {
				current[st] = 0;
				cont = 1;
				st++;
			}
		} while(cont && st<lattice->size);
	} while(st<lattice->size);
	free((void*)current);
	free((void*)max);
	return matchine;
}
