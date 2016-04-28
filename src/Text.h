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




#ifndef TextF
#define TextF
#include <stdio.h>

typedef unsigned char TypeSymbol;
typedef unsigned int TypeCardinal;

/*return the text read from f*/
char *readText(FILE *f);
/*return the alphabet of chains text and motif (no problem if one of them is NULL)*/
char *getAlphabet(char *text, char *motif);
TypeSymbol *toSymbolSequence(char *text, int length, char *alphabet);
char *toCharSequence(TypeSymbol *text, int length, char *alphabet);
/*print the symbol sequence text*/
void fprintSymbolSequence(FILE *f, TypeSymbol *text, int length, char *alphabet);
/*print hexadecimal*/
void fprintSymbolSequenceHexa(FILE *f, TypeSymbol *text, int length, int size);
/*sprint hexadecimal dest have to be allocated*/
char *sprintSymbolSequenceHexa(char *dest, TypeSymbol *text, int length, int size);
/*fill text with a random text of length size under model*/
void fillRandomPos(TypeSymbol *text, int size, TypeSymbol *textL, int sizeL);
/*return the symbols read from text file f*/
TypeSymbol *readFromTextFile(FILE *f, int *length, char **alphabet);
/*return the symbols read from binary file f*/
TypeSymbol *readFromBinaryFile(FILE *f, int *length);
/* get File size : could be used for text files as well */
int getBinaryFileSize(FILE *f);
/* read the binary content from file f: could be used for text files as well */
unsigned char *readBinary(FILE *f, int size);

#endif
