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




#ifndef LatticeF
#define LatticeF

#include <limits.h>
#include <stdio.h>
#include "Text.h"
#include "PowerSet.h"

#define SINK UINT_MAX




typedef struct PARTIAL_LATTICE {
	TypeCardinal cardinal;
	TypeSymbol *xterm;
	int length, ***shift, *ntrans, **next;
	TypeState ***trans, size, *term, sizeTerm;
} TypePartialLattice;

typedef struct LATTICE {
	TypeCardinal cardinal;
	int length, ***shift;
	TypeState ***trans, size;
} TypeLattice;

TypeLattice *getLattice(TypeCardinal cardinal, TypeSymbol *pattern, int length);
void freeLattice(TypeLattice *lattice);

void fprintLattice(FILE *f, TypeLattice *lattice, char *alphabet);
void fprintLatticeDot(FILE *f, TypeLattice *lattice, char *alphabet);
void fprintLatticeGasTex(FILE *f, TypeLattice *lattice, char *alphabet);
void fprintTransition(FILE *f, TypeState st, TypeLattice *lattice);
void fprintState(FILE *f, TypeState st, int length);
void fprintPartialLattice(FILE *f, TypePartialLattice *lattice, char *alphabet);
/*print partial lattice in dot format*/
void fprintPartialLatticeDot(FILE *f, TypePartialLattice *lattice, char *alphabet);
void setToState(int *set, int size, TypeState *st);
void stateToSet(TypeState st, int length, int *set, int *size);
int stateCardinal(TypeState st, int length);

TypePartialLattice *getPartialLatticePowerSet(TypeCardinal cardinal, TypeSymbol *pattern, int length, int K);
void freePartialLattice(TypePartialLattice *lattice);

#endif
