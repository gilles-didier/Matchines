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




#ifndef PatternMatchingAlgorithmsF
#define PatternMatchingAlgorithmsF

#include "MatchingMachine.h"

TypeMatchingMachine *getNaiveMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length);

TypeMatchingMachine *getMorrisPrattMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length);
TypeMatchingMachine *getKnuthMorrisPrattMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length);

TypeMatchingMachine *getMatchAutomataMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length);

TypeMatchingMachine *getQuickSearchMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length);
TypeMatchingMachine *getHorspoolMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length);

TypeMatchingMachine *getSAMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length);
TypeMatchingMachine *getFJSMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length);
TypeMatchingMachine *getTVSBMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length);
TypeMatchingMachine *getEBOMMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length);
TypeMatchingMachine *getBOMMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length);
TypeMatchingMachine *getHASHqMatchingMachine(TypeCardinal cardinal, TypeSymbol *pattern, int length);
#endif
