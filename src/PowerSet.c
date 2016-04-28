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

#include "PowerSet.h"
#include "Utils.h"




#define INC_LEXI_TREE_SET 200;

/*return the sum from i=0 to k of binomial coeff (n i)*/
static TypeState sumbincoef(int n, int k);
/*return the binomial coeff (n k)*/
static TypeState bincoef(int n,int k);


/*compare function of TypeState for qsort*/
int compareState(const void* a, const void* b) {
	if(*((TypeState*)a)>*((TypeState*)b))
		return 1;
	if(*((TypeState*)a)<*((TypeState*)b))
		return -1;
	return 0;
}

/*compare function of TypeState for qsort in decreasing order*/
int compareStateDec(const void* a, const void* b) {
	if(*((TypeState*)a)<*((TypeState*)b))
		return 1;
	if(*((TypeState*)a)>*((TypeState*)b))
		return -1;
	return 0;
}

/*recursively fill the terminal states list*/
void getTermIndexRec(TypeState n, int depth, TypeLexiTreeSet *dict, TypeState *list, TypeState *size) {
	if(depth==dict->length)
		list[(*size)++] = dict->node[n].child;
	else {
		TypeState c;
		for(c=dict->node[n].child; c!=END_SET; c=dict->node[c].sibling)
			getTermIndexRec(c, depth+1, dict, list, size);
	}
}

/*return all the terminal states in table list*/
void getTermIndex(TypeLexiTreeSet *dict, TypeState **list, TypeState *size) {
	*size = 0;
	*list = (TypeState*) malloc(dict->size*sizeof(TypeState));
	getTermIndexRec(dict->root, 0, dict, *list, size);
	*list = (TypeState*) realloc(*list, *size*sizeof(TypeState));
}

/*return the index of the set w in dict*/
TypeState findSetLexiTreeSet(int *w, TypeLexiTreeSet *dict) {
	int i;
	TypeState n;
	n = dict->root;
	for(i=0; i<dict->length;) {
		TypeState c;
		c = dict->node[n].child;
		while(c!=END_SET && dict->node[c].pos<w[i])
			c = dict->node[c].sibling;
		if(c!=END_SET && dict->node[c].pos == w[i]) {
			n = c;
			i++;
		} else
			return END_SET;
	}
	return dict->node[n].child;
}

/*add word w in dict. If w is already in, then it returns its index (stored in the child field of a leaf - labelled by '\0'), -1 otherwise*/
int addSetLexiTreeSet(int *w, TypeState index, TypeLexiTreeSet *dict) {
	int i;
	TypeState n;
	n = dict->root;
	for(i=0; i<dict->length;) {
		TypeState *prec, c;
		prec = &(dict->node[n].child);
		c = dict->node[n].child;
		while(c!=END_SET && dict->node[c].pos<w[i]) {
			prec = &(dict->node[c].sibling);
			c = dict->node[c].sibling;
		}
		if(c!=END_SET && dict->node[c].pos == w[i]) {
			n = c;
			i++;
		} else {
			*prec = dict->size;
			if(dict->size >= dict->sizeBuf) {
				dict->sizeBuf += INC_LEXI_TREE_SET;
				dict->node = (TypeLexiTreeSetNode*) realloc((void*)dict->node, dict->sizeBuf*sizeof(TypeLexiTreeSetNode));
			}
			initLexiTreeSetNode(w[i], &(dict->node[dict->size]));
			dict->node[dict->size].sibling = c;
			n = dict->size;
			dict->size++;
			i++;
			while(i<dict->length) {
				if(dict->size >= dict->sizeBuf) {
					dict->sizeBuf += INC_LEXI_TREE_SET;
					dict->node = (TypeLexiTreeSetNode*) realloc((void*)dict->node, dict->sizeBuf*sizeof(TypeLexiTreeSetNode));
				}
				initLexiTreeSetNode(w[i], &(dict->node[dict->size]));
				dict->node[n].child = dict->size;
				n = dict->node[n].child;
				dict->size++;
				i++;
			}
		}
	}
	if(dict->node[n].child == END_SET)
		dict->node[n].child = index;
	else {
		if(dict->node[n].child != index) {
			printf("Duplicate index in LexiTreeSet\n");
			exit(1);
		}
	}
	return 0;
}

/*basic init of node n of dict*/
void initLexiTreeSetNode(int pos, TypeLexiTreeSetNode *n) {
	n->pos = pos;
	n->child = END_SET;
	n->sibling = END_SET;
}

/*allocate a "lexi tree set"*/
TypeLexiTreeSet *newLexiTreeSet(int length) {
	TypeLexiTreeSet* dict;
	dict = (TypeLexiTreeSet*) malloc(sizeof(TypeLexiTreeSet));
	dict->sizeBuf = INC_LEXI_TREE_SET;
	dict->node = (TypeLexiTreeSetNode*) malloc(dict->sizeBuf*sizeof(TypeLexiTreeSetNode));
	dict->length = length;
	dict->root = 0;
	initLexiTreeSetNode(0,&(dict->node[dict->root]));
	dict->size = 1;
	return dict;
}

/*desallocate a "lexi tree set"*/
void freeLexiTreeSet(TypeLexiTreeSet *dict) {
	if(dict == NULL)
		return;
	if(dict->node != NULL)
		free((void*)dict->node);
	free((void*) dict);
}

/*return the sum from i=0 to k of binomial coeff (n i)*/
TypeState sumbincoef(int n, int k) {
	TypeState i, bincoef = 1, sum = 1;
	for(i=1; i<=k ; ++i) {
		bincoef = bincoef*(n-i+1)/i;
		sum += bincoef;
	}
	return sum;
}

/*return the binomial coeff (n k)*/
TypeState bincoef(int n,int k) {
	TypeState ans=1;
    k=k>n-k?n-k:k;
    int j;
    for(j=1; j<=k; j++, n--) {
		if(n%j==0) {
            ans*=n/j;
		} else {
			if(ans%j==0) {
				ans=ans/j*n;
			} else {
				ans=(ans*n)/j;
			}
		}
	}
    return ans;
}

/*recursive function of print dict*/
void fprintLexiTreeSetSpecialRec(FILE *f, int *tmp, int i, TypeState n, TypeLexiTreeSetSpecial *dict) {
	int j;
	for(j=0; j<i; j++)
		fprintf(f, "%d, ", tmp[j]);
	fprintf(f, " |%ld| card %d, min %d, suff %lu\n", n, dict->node[n].cardinal, dict->node[n].min, dict->node[n].suff);
	if(dict->node[n].link != NULL)
		for(j=dict->node[n].pos+1; j<dict->length; j++) {
			tmp[i] = dict->node[dict->node[n].link[j]].pos;
			fprintLexiTreeSetSpecialRec(f, tmp, i+1, dict->node[n].link[j], dict);
		}
}

/*print dict*/
void fprintLexiTreeSetSpecial(FILE *f, TypeLexiTreeSetSpecial *dict) {
	int tmp[MAX_DICT_LENGTH];
	fprintLexiTreeSetSpecialRec(f, tmp, 0, 0L, dict);
	TypeState n;
	for(n=0; n<dict->size; n++) {
		int i;
		fprintf(f, "node %lu\t", n);
		if(dict->node[n].link != NULL) {
			for(i=1; i<dict->length; i++)
				fprintf(f, " %lu", dict->node[n].link[i]);
		} else
			fprintf(f, "empty");
		fprintf(f, "\n");
	}	
}

/*print the set n of dict*/
void fprintSetLexiTreeSetSpecial(FILE *f, TypeState n, TypeLexiTreeSetSpecial *dict) {
	int tmp[MAX_DICT_LENGTH], i;
	tmp[dict->node[n].cardinal] = END_SET;
	for(i=dict->node[n].cardinal-1; i>=0; i--) {
		tmp[i] = dict->node[n].pos;
		n = dict->node[n].ancestor;
	}
	for(i=0; tmp[i] != END_SET; i++)
		fprintf(f, "%d ", tmp[i]);
}

/*return the index of the set w*/
TypeState findWordLexiPos(int *w, TypeLexiTreeSetSpecial *dict) {
	int i;
	TypeState n;
	n = 0L;
	for(i=0; w[i] != END_SET; i++)
		n = dict->node[n].link[w[i]];
	return n;
}

/*add word w in dict. If w is already in, then it returns its index (stored in the child field of a leaf - labelled by '\0'), -1 otherwise*/
void addSetLexiPos(int *w, TypeState index, TypeLexiTreeSetSpecial *dict) {
	int i;
	TypeState n = 0L;
	if(w[0] == END_SET)
		return;
	for(i=0; w[i+1] != END_SET; i++) /* be careful we don't check if the lexicographic is OK*/
		n = dict->node[n].link[w[i]];
	if(dict->node[n].link == NULL) {
		dict->node[n].link = (TypeState*) malloc(dict->length*sizeof(TypeState));
		if(dict->node[n].pos>=0)
			dict->node[n].link[dict->node[n].pos] = n;
	}
	dict->node[n].link[w[i]] = index;
	dict->node[index].ancestor = n;
	dict->node[index].pos = w[i];
	dict->node[index].cardinal = i+1;
	if(n == 0L)
		dict->node[index].min = w[i];
	else
		dict->node[index].min = dict->node[n].min;
	if(n == 0L || (dict->node[n].suff == 0L && dict->node[index].pos == dict->node[n].pos+1))
		dict->node[index].suff = 0L;
	else
		dict->node[index].suff = dict->node[dict->node[n].suff].link[dict->node[index].pos];
}

/*fill the set 'set' with the positions of index n*/
void fillSetLexiPos(int *set, TypeState n, TypeLexiTreeSetSpecial *dict) {
	int i;
	set[dict->node[n].cardinal] = END_SET;
	for(i=dict->node[n].cardinal-1; i>=0; i--, n=dict->node[n].ancestor)
		set[i] = dict->node[n].pos;
}

/*complete all the links/transitions from al sets n*/
void completeLexiPos(TypeState n, int K, TypeLexiTreeSetSpecial *dict) {
	if(dict->node[n].cardinal >= K) {
		return;
	} else {
		int j;
		if(dict->node[n].link == NULL) {
			dict->node[n].link = (TypeState*) malloc(dict->length*sizeof(TypeState));
			if(dict->node[n].pos >= 0)
				dict->node[n].link[dict->node[n].pos] = n;
		}
		for(j=0; j<dict->node[n].pos; j++)
			dict->node[n].link[j] = dict->node[dict->node[dict->node[n].ancestor].link[j]].link[dict->node[n].pos];
		for(j=dict->node[n].pos+1; j<dict->length; j++)
			completeLexiPos(dict->node[n].link[j], K, dict);
	}
}

/*allocate a lexi tree set of subset of {0,1,..., length-1}*/
TypeLexiTreeSetSpecial *newLexiTreeSetSpecial(int length, TypeState size) {
	TypeLexiTreeSetSpecial* dict;
	TypeState index;
	dict = (TypeLexiTreeSetSpecial*) malloc(sizeof(TypeLexiTreeSetSpecial));
	dict->size = size;
	dict->node = (TypeLexiTreeSetSpecialNode*) malloc(dict->size*sizeof(TypeLexiTreeSetSpecialNode));
	dict->length = length;
	dict->node[0].cardinal = 0;
	dict->node[0].pos = 0;
	dict->node[0].ancestor = 0L;
	dict->node[0].index = 0L;
	dict->node[0].min = 0;
	dict->node[0].suff = 0L;
	for(index=0L; index<dict->size; index++)
		dict->node[index].link = NULL;
	return dict;
}

/*free dict*/
void freeLexiTreeSetSpecial(TypeLexiTreeSetSpecial *dict) {
	if(dict != NULL) {
		if(dict->node != NULL) {
			TypeState n;
			for(n=0L; n<dict->size; n++)
				if(dict->node[n].link != NULL)
					free((void*)dict->node[n].link);			
			free((void*)dict->node);
		}
		free((void*) dict);
	}
}

/*print the set/table set*/
void fprintSet(FILE *f, int *set) {
	int i;
	for(i=0; set[i]!=END_SET; i++)
		fprintf(f, "%d ", set[i]);
}

/*print the K-power set ps*/
void fprintKPowerSet(FILE *f, TypePowerSet *ps) {
	TypeState st;
	for(st=0L; st<ps->size; st++) {
		fprintf(f, "%ld\tp %d\t|\t{", st, ps->prefix[st]);
		fprintSetLexiTreeSetSpecial(f, st-ps->start[ps->prefix[st]], ps->dict);
		fprintf(f, "}\n");
	}
}

/*print the state st of ps*/
void fprintStateKPowerSet(FILE *f, TypeState st, TypePowerSet *ps) {
	fprintf(f, "%lu\tp %d {", st, ps->prefix[st]);
	fprintSetLexiTreeSetSpecial(f, st-ps->start[ps->prefix[st]], ps->dict);
	fprintf(f, "}\n");
}

/*return the index of the set/table set*/
TypeState getStateKPowerSet(int prefix, int *set, TypePowerSet *ps) {
	return ps->start[prefix]+findWordLexiPos(set, ps->dict);
}

/*return a set where the last position of prefix is removed (always works)*/
TypeState shortenPrefixKPowerSet(TypeState st, TypePowerSet *ps) {
	if(ps->prefix[st]>0)
		return st-ps->start[ps->prefix[st]]+ps->start[ps->prefix[st]-1];
	return 0L;
}

/*return a set by adding a position to the prefix - does not work in all case be careful by using it*/
TypeState increasePrefixKPowerSet(TypeState st, TypePowerSet *ps) {
	if(ps->prefix[st]<ps->length-1)
		return st-ps->start[ps->prefix[st]]+ps->start[ps->prefix[st]+1];
	return 0L;
}

/*return the set st union{pos}*/
TypeState addPosKPowerSet(TypeState st, int pos, TypePowerSet *ps) {
	TypeState n = st-ps->start[ps->prefix[st]];
	if(pos > ps->prefix[st]) {
		return ps->dict->node[n].link[pos]+ps->start[ps->prefix[st]];
	} else {
		if(pos == ps->prefix[st]) {
			if(ps->dict->node[n].min == pos+1)
				return ps->dict->node[n].suff+ps->start[ps->prefix[st]+ps->dict->node[n].cardinal-ps->dict->node[ps->dict->node[n].suff].cardinal+1];
			else
				return n+ps->start[ps->prefix[st]+1];
		} else
			return st;
	}
	return st;
}

/*return the (index of) the set st minus its last position, if st not empty*/
TypeState removeLastKPowerSet(TypeState st, TypePowerSet *ps) {
	TypeState n = st-ps->start[ps->prefix[st]];
	if(n > 0L)
		return ps->dict->node[n].ancestor+ps->start[ps->prefix[st]];
	else {
		if(ps->prefix[st]>0)
			return ps->start[ps->prefix[st]-1];
		else
			return st;
	}
	return st;
}

/*return the index of the prefix set of length pos*/
TypeState getPrefixIndexKPowerSet(int pos, TypePowerSet *ps) {
	return ps->start[pos];
}

/*return the complementary set of st*/
void getCompKPowerSet(int *comp, TypeState st, TypePowerSet *ps) {
	int i, j, end;
	TypeState n = st-ps->start[ps->prefix[st]];
	end = ps->length-1;
	j = ps->length-ps->prefix[st]-ps->dict->node[n].cardinal;
	comp[j--]= END_SET;
	for(; n!=0L; n=ps->dict->node[n].ancestor) {
		for(i=end; i>ps->dict->node[n].pos; i--)
			comp[j--] = i;
		end = ps->dict->node[n].pos-1;
	}
	for(i=end; i>=ps->prefix[st]; i--)
		comp[j--] = i;
}

/*return the relative index of pos wrt, i.e its order in the complementary of st*/
int getTransIndexPowerSet(int pos, TypeState st, TypePowerSet *ps) {
	TypeState n;
	int minus;
	if(pos<ps->prefix[st])
		return pos;
	n = st-ps->start[ps->prefix[st]];
	minus = ps->dict->node[n].cardinal;
	for(; n!=0L && ps->dict->node[n].pos>pos; n=ps->dict->node[n].ancestor)
		minus--;
	return pos-minus-ps->prefix[st];
}

/*free ps*/
void freeKPowerSet(TypePowerSet *ps) {
	int i;
	free((void*)ps->prefix);
	free((void*)ps->start);
	free((void*)ps->next);
	for(i=1; i<ps->length; i++)
		free((void*)ps->begin[i]);
	free((void*)ps->begin);
	freeLexiTreeSetSpecial(ps->dict);
	free((void*)ps);
}

/*return the K-power set of a pattern og length 'length'*/
TypePowerSet *newKPowerSet(int length, int K) {
	TypePowerSet *ps;
	int *set, i, k, end, cont = 1;
	TypeState ind, tot, *tmp;
	ps = (TypePowerSet *) malloc(sizeof(TypePowerSet));
	ps->K = K;
	ps->length = length;
	ps->start = (TypeState*) malloc((length+1)*sizeof(TypeState*));
	ps->size = 0L;
	ps->start[0] = 0L;
	for(i=1; i<length; i++)
		ps->start[i] = ps->start[i-1]+sumbincoef(length-i, MIN(length-i, K));
	ps->start[length] = ps->start[length-1]+1;
	ps->size = ps->start[length];
	ps->prefix = (int*) malloc(ps->size*sizeof(int));
	for(i=0; i<length; i++)
		for(ind=ps->start[i]; ind<ps->start[i+1]; ind++)
			ps->prefix[ind] = i;
	tot = 0L;
	for(k=1; k<=K; k++) {
		tot += (k+1)*bincoef(length-1, k);
	}
	tot++;
	ps->dict = newLexiTreeSetSpecial(length, ps->start[1]);
	ps->next = (TypeState*) malloc(ps->start[1]*sizeof(TypeState));
	ps->begin = (TypeState**) malloc(ps->length*sizeof(TypeState*));
	for(i=1; i<ps->length; i++) {
		ps->begin[i] = (TypeState*) malloc((MIN(K,ps->length-i)+1)*sizeof(TypeState));
		for(k=0; k<=MIN(K,ps->length-i); k++)
			ps->begin[i][k] = END_STATE;
	}
	tmp = (TypeState*) malloc((K+1)*sizeof(TypeState));
	for(k=0; k<=K; k++)
		tmp[k] = END_STATE;
	set = (int*) malloc(length*sizeof(int));
	end = 0;
	ind = 1;
	set[0] = length-1;
	set[1] = END_SET;
	if(K>0) {
		do {
			int j;
			addSetLexiPos(set, ind, ps->dict);
			ps->next[ind] = tmp[end+1];
			tmp[end+1] = ind;
			if(set[end]<length-1 && end < K-1) {
				set[++end] = length-1;
			} else {
				for(j=end; j>0 && set[j]==(set[j-1]+1); j--);
				if(j>0 || set[0]>1) {
					set[j]--;
					end=j;
				} else
					cont = 0;
			}
			set[end+1] = END_SET;
			ind++;
		} while(cont);
	}
	free((void*)set);
	completeLexiPos(0L, K, ps->dict);
	for(k=1; k<=K && tmp[k]!=END_STATE; k++) {
		TypeState n;
		int min = ps->dict->node[tmp[k]].min;
		ps->begin[min][k] = tmp[k];
		for(n=tmp[k]; n<END_STATE; n=ps->next[n]) {
			if(ps->dict->node[n].min>min) {
				min = ps->dict->node[n].min;
				ps->begin[min][k] = n;
			}
		}
	}
	free((void*)tmp);
	return ps;
}
