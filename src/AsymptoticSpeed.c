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




#include "AsymptoticSpeed.h"
# include <superlu/slu_ddefs.h>
#include "PowerSet.h"
#include "Utils.h"


typedef struct ASYMP_SPEED_TRANS {
	double trans;
	TypeState state, next;
} TypeAsymptoticSpeedTrans;

/* acc[i] ith value of the matrix from left to right and top to bottom
 * icc[i] row of the ith value of the matrix
 * ccc[c] index of the first value of the cth column
 * 
 * A[j,i] = p(i->j)
 */
static void parseTarjan(TypeState st, TypeState *n, TypeState *class, TypeState c, TypeState *stack, TypeState *is, int *inStack, TypeState *num, TypeState *numAcc, TypeMatchingMachine *matchine, TypeBernoulli *bernoulli);

/*depth first parsing of the matching machine to get the communicating classes*/
void parseTarjan(TypeState st, TypeState *n, TypeState *class, TypeState c, TypeState *stack, TypeState *is, int *inStack, TypeState *num, TypeState *numAcc, TypeMatchingMachine *matchine, TypeBernoulli *bernoulli) {
	TypeCardinal a;
	num[st] = *n;
	numAcc[st] = *n;
	(*n)++;
	stack[++(*is)] = st;
	inStack[st] = 1;
	for(a=0; a<matchine->cardinal; a++) {
		if(matchine->trans[st][a] != SINK && bernoulli->prob[a]>0.) {
			if(num[matchine->trans[st][a]] == SINK) {
				parseTarjan(matchine->trans[st][a], n, class, c, stack, is, inStack,  num,  numAcc, matchine, bernoulli);
				numAcc[st] = MIN(numAcc[st], numAcc[matchine->trans[st][a]]);
			} else {
				if(inStack[matchine->trans[st][a]])
					numAcc[st] = MIN(numAcc[st], num[matchine->trans[st][a]]);
			}
		}
	}
	if(num[st] == numAcc[st]) {
		TypeState sti;
		c++;
		do {
			sti = stack[(*is)--];
			inStack[sti] = 0;
			class[sti] = c;
		} while(sti != st);
	}
}

/*return a table C where C[st] is the communicating class of st if C[st]>0 and st is transient if C[st]=0*/
TypeState *getCommunicatingClasses(TypeMatchingMachine *matchine, TypeBernoulli *bernoulli) {
	TypeState *stack, is = -1, *num, *numAcc, *class, st, n = 0;
	int *inStack;
	
	class = (TypeState*) malloc(matchine->size*sizeof(TypeState));
	stack = (TypeState*) malloc(matchine->size*sizeof(TypeState));
	num = (TypeState*) malloc(matchine->size*sizeof(TypeState));
	numAcc = (TypeState*) malloc(matchine->size*sizeof(TypeState));
	inStack = (int*) malloc(matchine->size*sizeof(int));
	for(st=0; st<matchine->size; st++) {
		num[st] = SINK;
		numAcc[st] = SINK;
		inStack[st] = 0;
		class[st] = 0;
	}
	parseTarjan(0, &n, class, 0, stack, &is, inStack, num, numAcc, matchine, bernoulli);
	free((void*)stack);
	free((void*)num);
	free((void*)numAcc);
	free((void*)inStack);
	return class;
}

/*return the asymptotic speed of matchine wrt model bernoulli*/
double getAsymptoticSpeed(TypeMatchingMachine *matchine, TypeBernoulli *bernoulli) {
	TypeState *class, st, max, min;
	double expect;
	class = getCommunicatingClasses(matchine, bernoulli);
	max = class[0];
	min = class[0];
	for(st=1; st<matchine->size; st++) {
		if(class[st]>max)
			max = class[st];
		if(class[st]<min)
			min = class[st];
	}
	if(max == 1) {
		TypeState *communicating, *index, effC = 0;
		communicating = (TypeState*) malloc(matchine->size*sizeof(TypeState));
		index = (TypeState*) malloc(matchine->size*sizeof(TypeState));
		for(st=0; st<matchine->size; st++)
			if(class[st] == 1) {
				index[st] = effC;
				communicating[effC] = st;
				effC++;
			}
		expect = getAsymptoticSpeedCommunicatingClass(communicating, effC, index, matchine, bernoulli);
		free((void*)communicating);
		free((void*)index);
	} else {
		TypeState *effC, effT, c, **communicating, *transient, *index;
		SuperMatrix A, B, L, U;
		superlu_options_t options;
		SuperLUStat_t stat;
		int *perm_c, *perm_r, info, *icc, *ccc;
		double *f, *transTmp, *acc;
		TypeState *target, nnz;

		if(class[0] != 0) {
			fprintf(stderr, "Execution error: not irreductible with 0 not transient\n");
			exit(1);
		}
		effC = (TypeState*) malloc(max*sizeof(TypeState));
		index = (TypeState*) malloc(matchine->size*sizeof(TypeState));
		for(c=0; c<max; c++)
			effC[c] = 0;
		effT = 0;
		for(st=0; st<matchine->size; st++)
			if(class[st]>0)
				index[st] = effC[class[st]-1]++;
			else
				index[st] = effT++;
		communicating = (TypeState**) malloc(max*sizeof(TypeState*));
		for(c=0; c<max; c++)
			communicating[c] = (TypeState*) malloc(effC[c]*sizeof(TypeState));
		transient = (TypeState*) malloc(effT*sizeof(TypeState));
		for(c=0; c<max; c++)
			effC[c] = 0;
		effT = 0;
		for(st=0; st<matchine->size; st++)
			if(class[st]>0)
				communicating[class[st]-1][effC[class[st]-1]++] = st;
			else
				transient[effT++] = st;
		f = (double*) malloc(max*c*sizeof(double*));
		for(st=0; st<effT*max; st++)
			f[st] = 0.;
		transTmp = (double*) malloc(effT*sizeof(double));
		for(st=0; st<effT; st++)
			transTmp[st] = 0.;
		target = (TypeState*) malloc((matchine->cardinal+1)*sizeof(TypeState));
		acc = (double*) malloc(effT*(matchine->cardinal+2)*sizeof(double));
		icc = (int*) malloc(effT*(matchine->cardinal+2)*sizeof(int));
		ccc = (int*) malloc((effT+1)*sizeof(int));
		nnz = 0L;
		ccc[st] = 0;
		for(st=0; st<effT; st++) {
			TypeCardinal a, size = 0;
			double sum = 0.;
			for(a=0; a<matchine->cardinal; a++)
				if(matchine->trans[transient[st]][a] != SINK && bernoulli->prob[a]>0.) {
					sum += bernoulli->prob[a];
					if(class[matchine->trans[transient[st]][a]]>0)
						f[(class[matchine->trans[transient[st]][a]]-1)*effT+st] += bernoulli->prob[a];
					else {
						target[size++] = index[matchine->trans[transient[st]][a]];
						transTmp[index[matchine->trans[transient[st]][a]]] -= bernoulli->prob[a];
					}
				}
			qsort(target, size, sizeof(TypeState), compareState);
			if(size>0) {
				TypeCardinal tmp = 1;
				for(a=1; a<size; a++)
					if(target[a-1] != target[a])
						target[tmp++] = target[a];
				size = tmp;
				for(a=0; a<size; a++)
					transTmp[target[a]] /= sum;
			}
			for(a=0; a<size && target[a]<st; a++)
				;
			if(a==size)
				target[size++] = st;
			else {
				if(target[a]>st) {
					TypeCardinal b;
					for(b=size-1; b<size; b--) {
						target[b+1] = target[b];
					}
					target[a] = st;
					size++;
				}
			}		
			transTmp[st] += 1.;	
			for(a=0; a<size; a++) {
				acc[nnz] = transTmp[target[a]];
				icc[nnz] = target[a];
				nnz++;
			}
			ccc[st+1] = nnz;
			for(a=0; a<size; a++)
				transTmp[target[a]] = 0.;
		}
		free((void*)target);
		free((void*)transTmp);
		free((void*)index);
		free((void*)transient);
		dCreate_CompCol_Matrix (&A, effT, effT, nnz, acc, icc, ccc, SLU_NC, SLU_D, SLU_GE);
		dCreate_Dense_Matrix (&B, effT, max, f, effT, SLU_DN, SLU_D, SLU_GE);
		perm_r = (int *) malloc (effT*sizeof(int));
		perm_c = (int *) malloc (effT*sizeof(int));
		set_default_options ( &options );
		options.ColPerm = NATURAL;
		StatInit(&stat);
		dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
		expect = 0.;
		for(c=0; c<max; c++)
			expect += f[c*effT]*getAsymptoticSpeedCommunicatingClass(communicating[c], effC[c], index, matchine, bernoulli);
		free((void*)effC);
		for(c=0; c<max; c++)
			free((void*)communicating[c]);
		free((void*)communicating);
		free((void*)f);
		free((void*)acc);
		free((void*)ccc);
		free((void*)icc);
		free((void*)perm_c);
		free((void*)perm_r);
		Destroy_SuperMatrix_Store(&A);
		Destroy_SuperMatrix_Store(&B);
		Destroy_SuperNode_Matrix(&L);
		Destroy_CompCol_Matrix(&U);
		StatFree(&stat);
	}
	free((void*)class);
	return expect;
}
		
/*return the asymptotic speed of a given communicating class of matchine wrt model bernoulli*/
double getAsymptoticSpeedCommunicatingClass(TypeState *comm, TypeState sizeComm, TypeState *index, TypeMatchingMachine *matchine, TypeBernoulli *bernoulli) {
	double *b, expect;
	SuperMatrix A, B, L, U;
	superlu_options_t options;
	int *perm_c, *perm_r, info;
	SuperLUStat_t stat;
	TypeCardinal a;
	double *transTmp, *acc;
	TypeState st, *target, nnz;
	int *icc, *ccc;


	transTmp = (double*) malloc(sizeComm*sizeof(double));
	target = (TypeState*) malloc((matchine->cardinal+1)*sizeof(TypeState));
	for(st=0; st<sizeComm; st++)
		transTmp[st] = 0.;
	acc = (double*) malloc(sizeComm*(matchine->cardinal+2)*sizeof(double));
	icc = (int*) malloc(sizeComm*(matchine->cardinal+2)*sizeof(int));
	ccc = (int*) malloc((sizeComm+1)*sizeof(int));
	nnz = 0L;
	ccc[0] = 0;
	for(st=0; st<sizeComm; st++) {
		TypeCardinal a, size = 0;
		double sum = 0.;
		for(a=0; a<matchine->cardinal; a++) {
			if(matchine->trans[comm[st]][a] != SINK && bernoulli->prob[a] > 0.) {
				transTmp[index[matchine->trans[comm[st]][a]]] += bernoulli->prob[a];
				target[size++] = index[matchine->trans[comm[st]][a]];
				sum += bernoulli->prob[a];
			}
		}
		if(size>0) {
			TypeCardinal tmp = 1;
			qsort(target, size, sizeof(TypeState), compareState);
			for(a=1; a<size; a++)
				if(target[a-1] != target[a])
					target[tmp++] = target[a];
			size = tmp;
			for(a=0; a<size; a++)
				transTmp[target[a]] /= sum;
		}
		for(a=0; a<size && target[a]<st; a++)
			;
		if(a==size)
			target[size++] = st;
		else {
			if(target[a]>st) {
				TypeCardinal b;
				for(b=size-1; b>=a && b<size; b--) {
					target[b+1] = target[b];
				}
				target[a] = st;
				size++;
			}
		}
		transTmp[st] -= 1.;	
		for(a=0; a<size; a++) {
			if(target[a] != sizeComm-1) {
				acc[nnz] = transTmp[target[a]];
				icc[nnz] = target[a];
				nnz++;
			}
		}
		acc[nnz] = 1.;
		icc[nnz] = sizeComm-1;
		nnz++;
		ccc[st+1] = nnz;
		for(a=0; a<size; a++)
			transTmp[target[a]] = 0.;
	}
	free((void*)target);
	free((void*)transTmp);
	dCreate_CompCol_Matrix (&A, sizeComm, sizeComm, nnz, acc, icc, ccc, SLU_NC, SLU_D, SLU_GE);
	b = (double*) malloc(sizeComm*sizeof(double));
	for(st=0; st<sizeComm-1; st++)
		b[st] = 0.;
	b[sizeComm-1] = 1;
	dCreate_Dense_Matrix (&B, sizeComm, 1, b, sizeComm, SLU_DN, SLU_D, SLU_GE);
	perm_r = (int *) malloc (sizeComm*sizeof(int));
	perm_c = (int *) malloc (sizeComm*sizeof(int));
	set_default_options ( &options );
	options.ColPerm = NATURAL;
	StatInit(&stat);
	dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
	expect = 0.;
/*for(st=0; st<sizeComm; st++)
	expect += b[st];
printf("exp %lf\n", expect);
*/	for(st=0; st<sizeComm; st++) {
		double sum = 0.;
		for(a=0; a<matchine->cardinal; a++)
			if(matchine->trans[comm[st]][a] != SINK)
				sum += bernoulli->prob[a];
		for(a=0; a<matchine->cardinal; a++) {
			if(matchine->trans[comm[st]][a] != SINK) {
				expect += ((double)matchine->shift[comm[st]][a])*b[st]*bernoulli->prob[a]/sum;
			}
		}
	}
	free((void*)acc);
	free((void*)ccc);
	free((void*)icc);
	free((void*)b);
	free((void*)perm_c);
	free((void*)perm_r);
	Destroy_SuperMatrix_Store(&A);
	Destroy_SuperMatrix_Store(&B);
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);
	StatFree(&stat);
	return expect;
}
