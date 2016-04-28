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




#ifndef PowerSetF
#define PowerSetF

#include <stdio.h>
#include <limits.h>

#define END_SET UINT_MAX
#define END_STATE ULONG_MAX
#define NO_STATE ULONG_MAX

typedef unsigned long TypeState;

typedef struct LEXI_TREE_SET_NODE {
	int pos;
	TypeState child, sibling;
} TypeLexiTreeSetNode;

typedef struct LEXI_TREE_SET {
	int length;
	TypeLexiTreeSetNode *node;
	TypeState root, size, sizeBuf;
} TypeLexiTreeSet;


typedef struct LEXI_TREE_SET_SPECIAL_NODE {
	int cardinal, pos, min; /*cardinal is the size of the set, pos its max and min its min*/
	TypeState index, ancestor, *link, suff; /*link[j] is the index of the union with {j}, suff is the set remaining by adding a prefix of size min*/
} TypeLexiTreeSetSpecialNode;

typedef struct LEXI_TREE_SET_SPECIAL {
	TypeLexiTreeSetSpecialNode *node;
	TypeState size;
	int length;
} TypeLexiTreeSetSpecial;

typedef struct POWER_SET {
	TypeLexiTreeSetSpecial *dict;
	int *prefix, length, K;
	TypeState *start, *next, **begin, size;
} TypePowerSet;









int compareState(const void* a, const void* b);
int compareStateDec(const void* a, const void* b);


void getTermIndexRec(TypeState n, int depth, TypeLexiTreeSet *dict, TypeState *list, TypeState *size);
void getTermIndex(TypeLexiTreeSet *dict, TypeState **list, TypeState *size);
TypeState findSetLexiTreeSet(int *w, TypeLexiTreeSet *dict);
/*add word w in dict. If w is already in, then it returns its index (stored in the child field of a leaf - labelled by '\0'), -1 otherwise*/
int addSetLexiTreeSet(int *w, TypeState index, TypeLexiTreeSet *dict);
void initLexiTreeSetNode(int pos, TypeLexiTreeSetNode *n);
TypeLexiTreeSet *newLexiTreeSet(int length);
void freeLexiTreeSet(TypeLexiTreeSet *dict);



void fprintLexiTreeSetSpecialRec(FILE *f, int *tmp, int i, TypeState n, TypeLexiTreeSetSpecial *dict);
void fprintLexiTreeSetSpecial(FILE *f, TypeLexiTreeSetSpecial *dict);
void completeLexiPos(TypeState n, int K, TypeLexiTreeSetSpecial *dict);
TypeLexiTreeSetSpecial *newLexiTreeSetSpecial(int length, TypeState size);
TypeState findWordLexiPos(int *w, TypeLexiTreeSetSpecial *dict);
void addLexiTreeSetSpecial(int *w, TypeState index, TypeLexiTreeSetSpecial *dict);
int indexWordLexiPos(int *w, TypeLexiTreeSetSpecial *dict);
void initLexiNodePos(int pos, TypeLexiTreeSetSpecialNode *n);
void fillLexiTreeSetSpecial(int *set, TypeState n, TypeLexiTreeSetSpecial *dict);
TypeState shortenPrefixKPowerSet(TypeState st, TypePowerSet *ps);
TypeState increasePrefixKPowerSet(TypeState st, TypePowerSet *ps);
TypeState getStateKPowerSet(int prefix, int *set, TypePowerSet *ps);
TypeState addPosKPowerSet(TypeState st, int pos, TypePowerSet *ps);
TypeState removeLastKPowerSet(TypeState st, TypePowerSet *ps);
TypeState getPrefixIndexKPowerSet(int pos, TypePowerSet *ps);
void getCompKPowerSet(int *comp, TypeState st, TypePowerSet *ps);
int getTransIndexPowerSet(int pos, TypeState st, TypePowerSet *ps);
void fprintSet(FILE *f, int *set);
void fprintKPowerSet(FILE *f, TypePowerSet *ps);
void fprintStateKPowerSet(FILE *f, TypeState st, TypePowerSet *ps);
void freeKPowerSet(TypePowerSet *ps);
TypePowerSet *newKPowerSet(int length, int K);
void fprintSetLexiTreeSetSpecial(FILE *f, TypeState n, TypeLexiTreeSetSpecial *dict);

#endif
