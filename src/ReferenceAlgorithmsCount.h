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




#ifndef ReferenceAlgorithmsCountF
#define ReferenceAlgorithmsCountF
#include "MatchingMachine.h"

int genericAlgorithmCount(TypeMatchingMachine *strat, int m, TypeSymbol *t, int n);
int searchNaiveCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n);
int searchHASHqCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n);
int searchKMPCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n);
int searchMPCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n);
int searchSACount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n);
int searchFJSCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n);
int searchTVSBSCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n);
int searchEBOMCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n);
int searchHorspoolCount(TypeCardinal cardinal, TypeSymbol *P, int m, TypeSymbol *T, int n);
int searchQuicksearchCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n);

#endif
