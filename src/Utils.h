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




#ifndef UtilsF
#define UtilsF


#include <stdio.h>
#include <sys/time.h>

typedef /*unsigned*/ long TypePosition;

typedef unsigned int TypeNumber;

typedef enum TS {
	ErrorArgument=0,
	ErrorInit,
	ErrorReading,
	ErrorWriting,
	ErrorMemory,
	ErrorExec,
	ErrorArgs,
	ExitOk
} TypeExit;

#define INC_SIZE_DICT 100
#define MAX_DICT_LENGTH 1000

typedef struct LIST_DOUBLE {
	int size;
	double *val;
} TypeListDouble;

/*if symbol ==-1 then child contains an index*/
typedef struct LEXI_NODE {
	char symbol;
	int child, sibling;
} TypeLexiNode;

typedef struct LEXI_TREE {
	TypeLexiNode *node;
	int root, size, sizeBuf;
} TypeLexiTree;

typedef struct DICT_NODE {
	char symbol;
	struct DICT_NODE *child;
	struct DICT_NODE *sibling;
	long index;
} TypeDictNode;

typedef struct DICT_TREE {
	TypeDictNode *root;
	long index;
} TypeDictTree;

typedef struct INDEX {
	int size, buffer;
	char **name;
	TypeDictNode *dict;
} TypeIndex;


double difftime_ms(struct timeval t1, struct timeval t0);

char *strdpl(char *src);

void replaceChar(char *s, char a, char n);
void removeSpaces(char *s);
void freeListDouble(TypeListDouble *list);

int findWordLexi(char *w, TypeLexiTree *dict);
int addWordLexi(char *w, int index, TypeLexiTree *dict);
void initLexiNode(char symbol, TypeLexiNode *n);
TypeLexiTree *newLexiTree();
void freeLexiTree(TypeLexiTree *dict);
void fprintLexiTree(FILE *f, TypeLexiTree *dict);
/*return a lexitree from tab name*/
TypeLexiTree *getDictNameTab(char **name, int size);

void initIndex(TypeIndex *species);
void printIndex(FILE *f, TypeIndex *index);
int addIndex(char *name, TypeIndex *species);

void freeDictNode(TypeDictNode *n);
TypeDictNode *newDictNode(char c);
int getIndexString(char *s, TypeDictNode *cur, int *size);
int getIndex(char *s, TypeDictNode *cur);
void exitProg(TypeExit code, char *message);
void *monmalloc(long size);
void *monrealloc(void *in, long size);
int IsSeparator(char c);
int IsItemSeparator(char c);
int IsLineSeparator(char c);
int readLine(FILE *f, char *buffer);
int readItem(FILE *f, char *buffer);
int tokenize(char *src, char **dest);
char nextStartLine(FILE *f, char c);
char nextStartItem(FILE *f, char c);
void skipSeparator(FILE *f);
char passLines(FILE *f, char c);
char passSpaces(FILE *f, char c);
int find(char *src, char **dest, int size);
void fixSpace(char *src);
char *truncFileName(char* name);
char *getExtension(char* name);

char skipLineSpaceComment(FILE *f, char c);
char skipSep(FILE *f, char c);
int isline(char c);
int issep(char c);
int issepline(char c);

int compareInt(const void* a, const void* b);
int compareDouble(const void* a, const void* b);

size_t *qsortindex(void *base, size_t nitems, size_t size, int (*compar)(const void *, const void*));
/*return the ordered table of index of base*/ 
size_t *qsortTable(void *base, size_t nitems, size_t size, int (*compar)(const void *, const void*));

#define SPECIAL UINT_MAX
#define SEP '|'
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define ABS(x) ((x)<(0)?(-x):(x))
#define POS(x) ((x)<(0)?(0):(x))
#define NEG_INFTY -9.999e99
#define POS_INFTY 9.999e99
//#define DEBUGX
#endif
