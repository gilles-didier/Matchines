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




#ifndef ExploreF
#define ExploreF

#include "MatchingMachine.h"
#include "Bernoulli.h"


/*return the next matching machine from current*/
TypeMatchingMachine *getMatchingMachineFromLatticeNext(TypeLattice *lattice, TypeSymbol *pattern, int *current);
/*return the fastest matching machine wrt bernoulli model*/
TypeMatchingMachine *explore(TypeCardinal cardinal, TypeSymbol *pattern, int length, TypeBernoulli *bernoulli);
/*return the next matching machine from current*/
TypeMatchingMachine *getMatchingMachineFromPartialLatticeNext(TypePartialLattice *lattice, int *current);
/*return the fastest matching machine from the partial lattice wrt bernoulli model*/
TypeMatchingMachine *explorePartialLattice(TypePartialLattice *lattice, TypeBernoulli *bernoulli);


#endif
