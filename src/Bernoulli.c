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




#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "Bernoulli.h"
#include "Utils.h"

#define MAX_SIZE 200

/*return a random text of length size under model*/
TypeSymbol *randomText(int size, TypeBernoulli *model) {
	int i;
	TypeSymbol *text = (TypeSymbol*) malloc(size*sizeof(TypeSymbol));
	for(i=0; i<size; i++)
		text[i] = getRandomSymbol(model);
	return text;
}

/*fill text with a random text of length size under model*/
void fillRandomText(TypeSymbol *text, int size, TypeBernoulli *model) {
	int i;
	for(i=0; i<size; i++)
		text[i] = getRandomSymbol(model);
}

/*estimate Bernoulli model from text*/
TypeBernoulli *getBernoulli(TypeCardinal cardinal, TypeSymbol *text, int length) {
	int i;
	TypeBernoulli *bernoulli;
	bernoulli = (TypeBernoulli*) malloc(sizeof(TypeBernoulli));
	bernoulli->cardinal = cardinal;
	bernoulli->prob = (double*) malloc(bernoulli->cardinal*sizeof(double));
	for(i=0; i<bernoulli->cardinal; i++)
		bernoulli->prob[i] = 0.;
	for(i=0; i<length; i++)
		if(text[i] < 0 || text[i] >= bernoulli->cardinal) {
			printf("Bad character %d/%d!\n", text[i], bernoulli->cardinal);
			exit(1);
		} else
			bernoulli->prob[text[i]]++;
	for(i=0; i<bernoulli->cardinal; i++)
		bernoulli->prob[i] /= (double) length;
	return bernoulli;
}

/*return the uniform Bernoulli model on an alphabet of size cardinal*/
TypeBernoulli *getUniformBernoulli(TypeCardinal cardinal) {
	int i;
	TypeBernoulli *bernoulli;
	if(cardinal == 0) {
		fprintf(stderr, "Uniform Bernoulli model with empty alphabet\n");
		exit(1);
	}
	bernoulli = (TypeBernoulli*) malloc(sizeof(TypeBernoulli));
	bernoulli->cardinal = cardinal;
	bernoulli->prob = (double*) malloc(bernoulli->cardinal*sizeof(double));
	bernoulli->prob[0] = 1./((double)bernoulli->cardinal);
	for(i=1; i<bernoulli->cardinal; i++)
		bernoulli->prob[i] = bernoulli->prob[0];
	return bernoulli;
}

/*read a Bernoulli text file*/
TypeBernoulli *freadBernoulli(FILE *f, char *alphabet) {
	int code[256], i, n=0;
	int c;
	char tmp[MAX_SIZE];
	TypeBernoulli *bernoulli;
	bernoulli = (TypeBernoulli*) malloc(sizeof(TypeBernoulli));
	bernoulli->cardinal = strlen(alphabet);
	bernoulli->prob = (double*) malloc(bernoulli->cardinal*sizeof(double));
	for(i=0; i<256; i++)
		code[i] = -1;
	for(i=0; alphabet[i] != '\0'; i++)
		code[(unsigned char) alphabet[i]] = i;
	do {
		int index;
		for(c=fgetc(f); c != EOF && isspace(c); c=fgetc(f));
		if(c != EOF) {
			index = code[c];
			code[c] = -1;
			if(index<0) {
				printf("bad or duplicated character\n");
				exit(1);
			}
			c = fgetc(f);
			if(!issep(c)) {
				printf("lines have to start with a single character\n");
				exit(1);
			}
			for(; c != EOF && issep(c); c = fgetc(f));
			for(i=0; c != EOF && !isspace(c) && i<MAX_SIZE-1; i++) {
				tmp[i] = c;
				c = fgetc(f);
			}
			if(i<MAX_SIZE-1)
				tmp[i++] = '\0';
			else {
				printf("number too long\n");
				exit(1);
			}
			bernoulli->prob[index] = atof(tmp);
			n++;
		}
	} while(c != EOF);
	if(n<bernoulli->cardinal) {
		printf("missing character(s)\n");
		exit(1);
	}
	return bernoulli;
}

/*write a Bernoulli model as text*/
void fprintBernoulli(FILE *f, TypeBernoulli *bernoulli, char *alphabet) {
	int i;
	for(i=0; i<bernoulli->cardinal; i++)
		if(alphabet != NULL)
			fprintf(f, "%c\t%lf\n", alphabet[i], bernoulli->prob[i]);
		else
			fprintf(f, "%d\t%lf\n", i, bernoulli->prob[i]);
}

/*desallocate a Bernoulli model*/
void freeBernoulli(TypeBernoulli *bernoulli) {
	if(bernoulli != NULL) {
		free((void*) bernoulli->prob);
		free((void*) bernoulli);
	}
}
		
/*return a symbol drawn according to model*/	
TypeSymbol getRandomSymbol(TypeBernoulli *model) {
	TypeCardinal c;
	double x = ((double) rand())/((double)RAND_MAX), y = 0;
	for(c=0; c<model->cardinal; c++) {
		y+=model->prob[c];
		if(x<y)
			return c;
	}
	return model->cardinal-1;
}

/*return the probability of pattern m under model*/
double getBernoulliProb(TypeSymbol *m, int length, TypeBernoulli *model) {
	double lprob = 0.;
	int i;
	for(i=0; i<length; i++)
		lprob += log(model->prob[m[i]]);
	return exp(lprob);
}
