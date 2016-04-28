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
#include <string.h>

#include "ReferenceAlgorithmsCount.h"
#include "Utils.h"


static int fillHASHqBadChar(TypeCardinal cardinal, TypeSymbol *pattern, int m, int *shift);
static void preKmp(TypeSymbol *x, int m, int kmpNext[]);
static void preMP(TypeSymbol *motif, int length, int border[]);
static void preSA(TypeCardinal cardinal, TypeSymbol *x, int m, unsigned int S[]);
static void preQsBc(TypeCardinal cardinal, TypeSymbol *x, int m, int *qbc);
static void TVSBSpreBrBc(TypeCardinal cardinal, TypeSymbol *x, int m, int **brBc);
static void Pre_Horspool(TypeCardinal cardinal, TypeSymbol *P, int m, int hbc[]);

int genericAlgorithmCount(TypeMatchingMachine *strat, int m, TypeSymbol *t, int n) {
	int count=0, q=0, p=0;
	while(p <= n-m && p+strat->next[q] < n) {
		int qnext;
		qnext = strat->trans[q][t[p+strat->next[q]]];
		p += strat->shift[q][t[p+strat->next[q]]];
		q = qnext;
		count++;
	}
	return count;
}


/* HASHq Lecroq Wu Manber */

#define RANK 3
#define WSIZE 255

int fillHASHqBadChar(TypeCardinal cardinal, TypeSymbol *pattern, int m, int *shift) {
	int i, r, sh1;
	TypeSymbol h;
	
	for(i=0; i<WSIZE; ++i)
		shift[i] = m-RANK+1;
	h = pattern[0];
	for(r=1; r<RANK; r++)
		h = ((h<<1) + pattern[r]);
	shift[h] = m-RANK;
	for(i=RANK; i<m-1; ++i) {
		h = pattern[i-RANK+1];
		for(r=1; r<RANK; r++)
			h = ((h<<1) + pattern[r+i-RANK+1]);
		shift[h%WSIZE] = m-1-i;
	}
	h = pattern[m-RANK];
	for(r=1; r<RANK; r++)
		h = ((h<<1) + pattern[r+m-RANK]);
	sh1 = shift[h%WSIZE];
	shift[h%WSIZE] = 0;
	if(sh1==0)
		sh1 = 1;
	return sh1;
}

int searchHASHqCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n) {
	int count, j, i, sh1, *shift;
	TypeSymbol h;
	if (m<3) return -1;
	count = 0;

	/* Preprocessing */
	shift = (int *) malloc (WSIZE * sizeof(int)); 
	sh1 = fillHASHqBadChar(cardinal, x, m, shift);

	/* Searching */
	i = m-1;

	while (i<n) {
		h = y[i-2];
		h = ((h<<1) + y[i-1]);
		h = ((h<<1) + y[i]);
		count += 3;
		if (shift[h] == 0) {
			for(j=0; j<m && x[j]==y[i-m+1+j]; j++);
			count += j;
			if (j <m) {
				count++;
			}
			i += sh1;
		} else {
			i += shift[h];
		}
	}
	return count;
} 


/*Naive*/

int searchNaiveCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n) {
	int i, j, count = 0;

	/* Searching */
	for(i=0; i<=n-m; i++) {
		for(j=0; j<m && x[j] == y[i+j]; j++)
			count++;
		if(j<m)
			count++;
	}
	return count;
}

/*KMP*/

void preKmp(TypeSymbol *x, int m, int kmpNext[]) {
   int i, j;

   i = 0;
   j = kmpNext[0] = -1;
   while(i<m) {
      while(j>-1 && x[i] != x[j])
         j = kmpNext[j];
      i++;
      j++;
      if(i<m && x[i] == x[j])
         kmpNext[i] = kmpNext[j];
      else
         kmpNext[i] = j;
   }
}


int searchKMPCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n) {
	int i, j, *kmpNext, count = 0;

	/* Preprocessing */
	kmpNext = (int*) malloc((m+1)*sizeof(int));
	preKmp(x, m, kmpNext);

	i = j = 0;
	while ((j-i) <= n-m) {
		while (i > -1 && (j-i) <= (n-m) && x[i] != y[j]) {
			i = kmpNext[i];
			count++;
		}
		if(i>-1 && (j-i) <= (n-m)) {
			count++;
		}
		i++;
		j++;
		if (i >= m) {
			i = kmpNext[i];
		}
	}

	free((void*)kmpNext);
	return count;
}


void preMP(TypeSymbol *motif, int length, int border[]) {
	int i, j;
	i = 0; 
	j = -1;
	border[0] = -1;
	while(i<length) {
		while(j>=0 && motif[i] != motif[j])
			j = border[j];
		border[++i] = ++j;
	}
}

int searchMPCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n) {
	int i, j, *kmpNext, count = 0;

	/* Preprocessing */
	kmpNext = (int*) malloc((m+1)*sizeof(int));
	preMP(x, m, kmpNext);
	/* Searching */
	i = j = 0;
	while ((j-i) <= n-m) {
		while (i > -1 && (j-i) <= (n-m) && x[i] != y[j]) {
			i = kmpNext[i];
			count++;
		}
		if(i>-1 && (j-i) <= (n-m)) {
			count++;
		}
		i++;
		j++;
		if (i >= m) {
			i = kmpNext[i];
		}
	}
	free((void*)kmpNext);
	return count;
}

/*SA*/

void preSA(TypeCardinal cardinal, TypeSymbol *x, int m, unsigned int S[]) { 
   unsigned int j; 
   int i; 
   for(i=0; i<cardinal; ++i) S[i] = 0; 
   for(i=0, j=1; i<m; ++i, j <<= 1) { 
      S[x[i]] |= j; 
   } 
} 

int searchSACount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n) { 
   unsigned int D; 
   unsigned int *S; 
   int j; 

	S = (unsigned int*) malloc(cardinal*sizeof(unsigned int));
   /* Preprocessing */ 
   preSA(cardinal, x, m, S); 

   /* Searching */ 
   for(D=0, j=0; j<n; ++j)
      D = ((D<<1) | 1) & S[y[j]];
   return n;
} 


/*FJS*/

void preQsBc(TypeCardinal cardinal, TypeSymbol *x, int m, int *qbc) {
	int i;
	for (i = 0; i < cardinal; i++)
		qbc[i] = m+1;
	for (i = 0; i < m; i++)
		qbc[x[i]] = m-i;
}

int searchFJSCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n) {
	int i, s, count, *qsbc, *kmp;
	
	qsbc = (int*) malloc(cardinal*sizeof(int));
	kmp = (int*) malloc((m+1)*sizeof(int));
	
   /* Preprocessing */
	preQsBc(cardinal, x,m,qsbc);
	preKmp(x,m,kmp);

	/* Searching */
	s = 0;
	count = 0;
	while(s<=n-m) {
		while(s<=n-m && x[m-1]!=y[s+m-1]) {
			if(s<n-m) {
				s+=qsbc[y[s+m]] ;
				count+=2;
			} else {
				s+=1;
				count++;
			}
		}
		if(s<=n-m) {
			count++;
			i=0; 
			while(i<m && x[i]==y[s+i]){
				i++;
			}
			count += i;
			if(i<m)
				count++;
			s+=(i-kmp[i]);
		}
	}
	free((void*)kmp);
	free((void*)qsbc);
	return count;
}


void TVSBSpreBrBc(TypeCardinal cardinal, TypeSymbol *x, int m, int **brBc) {
   int i;
   TypeCardinal a, b;
   for(a=0; a<cardinal; ++a)
      for (b=0; b<cardinal; ++b)
         brBc[a][b] = m+2;
   for(a=0; a<cardinal; ++a)
      brBc[a][x[0]] = m+1;
   for(i=0; i<m-1; ++i)
      brBc[x[i]][x[i+1]] = m-i;
   for(a=0; a<cardinal; ++a)
      brBc[x[m-1]][a] = 1;
}


int searchTVSBSCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *z, int n){
	int count,i,j =0, **BrBc;
	TypeSymbol firstCh, lastCh, *y;
	count = 0;
	BrBc = (int**) malloc(cardinal*sizeof(int*));
	for(i=0; i<cardinal; i++)
		BrBc[i] = (int*) malloc(cardinal*sizeof(int));
	TVSBSpreBrBc(cardinal, x, m, BrBc);
	firstCh = x[0];
	lastCh = x[m-1];
	y = (TypeSymbol *) malloc((n+2*m)*sizeof(TypeSymbol)); 
	memcpy((void*)y, (void*)z, n*sizeof(TypeSymbol));
	memcpy((void*)(y+n), (void*)x, m*sizeof(TypeSymbol));
	memcpy((void*)(y+n+m), (void*)x, m*sizeof(TypeSymbol));
//	for(i=0; i<m; i++)
//		y[n+i] = y[n+m+i] = x[i];
	while(j<=n-m){
		if(lastCh == y[j+m-1]) {
			count++;
			if (firstCh==y[j]) {
				count++;
				for (i=m-2; i>0 && x[i]==y[j+i]; i--)
					count++;
				if(i>0)
					count++;
			} else {
				count++;
			}
		} else {
			count++;
		}
		j += BrBc[y[j+m]][y[j+m+1]];
		count += 2;
	}
	for(i=0; i<cardinal; i++)
		free((void*)BrBc[i]);
	free((void*)BrBc);
	free((void*)y);
	return count;
 }

/*EBOM*/

int searchEBOMCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *z, int n) {
	int *S, **FT;
	int **trans;
	int i, j, p, q=1;
	int iMinus1, mMinus1, count;
	TypeCardinal a, b, c;
	TypeSymbol *y;
	count = 0;

	y = (TypeSymbol *) malloc((n+m+1)*sizeof(TypeSymbol)); 
	memcpy((void*)y, z, n*sizeof(TypeSymbol));
	memcpy((void*)(y+n), x, m*sizeof(TypeSymbol));
	S = (int*) malloc((m+1)*sizeof(int));
	FT = (int**) malloc(cardinal*sizeof(int*));
	for(a=0; a<cardinal; a++)
		FT[a] = (int*) malloc(cardinal*sizeof(int));
   // Preprocessing
	trans = (int**) malloc((m+2)*sizeof(int*));
	for(i=0; i<=m+1; i++) {
		trans[i] = (int *)malloc (sizeof(int)*(cardinal));
	}
	for(i=0; i<=m+1; i++) {
		for (a=0; a<cardinal; a++) {
			trans[i][a]=SPECIAL;
		}
	}
	S[m] = m+1;
	for(i=m; i>0; --i) {
		iMinus1 = i - 1;
		c = x[iMinus1];
		trans[i][c] = iMinus1;
 		p = S[i];
		while(p <= m && (q = trans[p][c]) ==  SPECIAL) {
			trans[p][c] = iMinus1;
			p = S[p];
		}
		S[iMinus1] = (p == m + 1 ? m : q);
	}

   /* Construct the FirstTransition table */
	for(a=0; a<cardinal; a++) {
		q = trans[m][a];
		for(b=0; b<cardinal; b++)
         if (q>=0) {
			FT[a][b] = trans[q][b];
         } else {
			FT[a][b] = SPECIAL;
		}
	}
   /* Searching */
	for(i=0; i<m; i++) y[n+i]=x[i]; 
	if(!memcmp(x,y,m)) count++;
	j=m;
	mMinus1 = m-1;
	while(j<n) {
		while(FT[y[j]][y[j-1]] == SPECIAL) {
			count += 2;
			j+=mMinus1;
		}
		count += 2;
		i = j-2;
		p = FT[y[j]][y[j-1]];
		while((p = trans[p][y[i]]) != SPECIAL) {
			count++;
			i--;
		}
		count++;
		if(i<j-mMinus1 && j<n) {
			i++;
		}
		j = i + m;
	}
    free((void*)S);
	for(i=0; i<cardinal; i++)
		free((void*)FT[i]);
	free((void*)FT);
	for(i=0; i<=m+1; i++)
		free((void*)trans[i]);
	free((void*)trans);
	free((void *)y);
	return count;
}


void Pre_Horspool(TypeCardinal cardinal, TypeSymbol *P, int m, int hbc[]) {
   int i;
   TypeCardinal a;
   for(a=0; a<cardinal; a++) hbc[a]=m;
   for(i=0; i<m-1; i++) hbc[P[i]]=m-i-1;
}


int searchHorspoolCount(TypeCardinal cardinal, TypeSymbol *P, int m, TypeSymbol *T, int n) {
	int i, s, count, *hbc;
	hbc = (int*) malloc(cardinal*sizeof(int));
	Pre_Horspool(cardinal, P, m, hbc);
	/* Searching */
	s = 0;
	count = 0;
	while(s<=n-m) {
		i=m-1;
		while(i>=0 && ++count && P[i]==T[s+i]) {
			i--;
		}
		s+=hbc[T[s+m-1]];
	}
	free((void*)hbc);
	return count;
}




int searchQuicksearchCount(TypeCardinal cardinal, TypeSymbol *x, int m, TypeSymbol *y, int n) {
	int i, j, count = 0, *qsBc;
	qsBc = (int*) malloc(cardinal*sizeof(int));

	/* Preprocessing */
	preQsBc(cardinal, x, m, qsBc);
	/* Searching */
	j = 0;
	while(j<=n-m) {
		i=0;
		while(i<m && ++count && x[i]==y[j+i]) {
			i++;
		}
		if (j < n-m) {
			j+=qsBc[y[j+m]];
			count++;
		} else {
			j+=1;
		}		
	}
	free((void*)qsBc);
	return count;
}



