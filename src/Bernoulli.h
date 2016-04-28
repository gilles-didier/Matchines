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




#ifndef BernoulliF
#define BernoulliF

#include <stdio.h>
#include "Text.h"

typedef struct BERNOUILLI {
	TypeCardinal cardinal;
	double *prob;
} TypeBernoulli;

/*estimate Bernoulli model from text*/
TypeBernoulli *getBernoulli(TypeCardinal cardinal, TypeSymbol *text, int length);

/*return the uniform Bernoulli model on an alphabet of size cardinal*/
TypeBernoulli *getUniformBernoulli(TypeCardinal cardinal);

/*read a Bernoulli text file*/
TypeBernoulli *freadBernoulli(FILE *f, char *alphabet);

/*write a Bernoulli model as text*/
void fprintBernoulli(FILE *f, TypeBernoulli *bernoulli, char *alphabet);

/*return the probability of pattern m under model*/
double getBernoulliProb(TypeSymbol *m, int length, TypeBernoulli *model);

/*write a Bernoulli model as text*/
TypeSymbol getRandomSymbol(TypeBernoulli *bernoulli);

/*return a random text of length size under model*/
TypeSymbol *randomText(int size, TypeBernoulli *model);

/*fill text with a random text of length size under model*/
void fillRandomText(TypeSymbol *text, int size, TypeBernoulli *model);

/*desallocate a Bernoulli model*/
void freeBernoulli(TypeBernoulli *bernoulli);

#endif
