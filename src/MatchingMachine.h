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




#ifndef MatchingMachineF
#define MatchingMachineF

#include <limits.h>
#include <stdio.h>
#include "Bernoulli.h"
#include "Lattice.h"
#include "PowerSet.h"

#define SINK UINT_MAX

#define FASTEST "Fastest"
#define UNIFORM "Unif."
#define HEURISTIC "-Heuristic"

typedef struct MATCHING_MACHINE {
	TypeCardinal cardinal;
	TypeSymbol *xterm;
	int *next, **shift;
	TypeState size, **trans, *term, sizeTerm;
} TypeMatchingMachine;

TypeMatchingMachine *getNewKMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length, TypeBernoulli *bernoulli, int K);
TypeMatchingMachine *getReverseMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length);
void freeMatchingMachine(TypeMatchingMachine *matchine);

void fprintMatchingMachine(FILE *f, TypeMatchingMachine *matchine, char *alphabet);
void fprintMatchingMachineDot(FILE *f, TypeMatchingMachine *matchine, char *alphabet);

TypeMatchingMachine *simplify(TypeMatchingMachine *matchine);

TypeMatchingMachine *expandMatchingMachine(TypeMatchingMachine *matchine);
void getMatchingMachineBounds(TypeMatchingMachine *matchine, int *minNext, int *maxNext, int *minShift, int *maxShift);
void fprintMatchingMachineCode(FILE *f, TypeMatchingMachine *matchine);
TypeMatchingMachine *getBestKFromPartial(TypePartialLattice *lattice, TypeBernoulli *bernoulli, int K);
TypeMatchingMachine *getPolyKMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length, TypeBernoulli *bernoulli, int K);

#endif
