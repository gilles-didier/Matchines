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
#include <sys/stat.h>

#include "Text.h"
#include "Utils.h"

#define INC_BUFFER_CHAR 500


/* get File size : could be used for text files as well */
int getBinaryFileSize(FILE *f) {
	struct stat st;
	fstat(fileno(f), &st);
	return st.st_size;
}

/* read the binary content from file f: could be used for text files as well */
unsigned char *readBinary(FILE *f, int size) {
	unsigned char *buffer;
	buffer = (unsigned char*) malloc(size*sizeof(unsigned char));
	size = fread(buffer, sizeof(char), size, f);
	return buffer;
}

/*return the text read from f*/
char *readText(FILE *f) {
	int sizeBuffer, size;
	char c, *text;
	sizeBuffer = INC_BUFFER_CHAR;
	text = (char*) malloc(sizeBuffer*sizeof(char));
	size = 0;
	while((c=getc(f)) != EOF) {
		if(size>=sizeBuffer) {
			sizeBuffer += INC_BUFFER_CHAR;
			text = (char*) realloc(text, sizeBuffer*sizeof(char));
		}
		text[size++] = c;
	}
	if(size>=sizeBuffer) {
		sizeBuffer += INC_BUFFER_CHAR;
		text = (char*) realloc(text, sizeBuffer*sizeof(char));
	}
	text[size++] = '\0';
	text = (char*) realloc(text, size*sizeof(char));
	return text;
}

/*return the symbols read from text file f*/
TypeSymbol *readFromTextFile(FILE *f, int *length, char **alphabet) {
	int sizeBuffer, size;
	char c, *text;
	TypeSymbol *symb;
	sizeBuffer = INC_BUFFER_CHAR;
	text = (char*) malloc(sizeBuffer*sizeof(char));
	size = 0;
	while((c=getc(f)) != EOF && (*length<0 || size<*length)) {
		if(size>=sizeBuffer) {
			sizeBuffer += INC_BUFFER_CHAR;
			text = (char*) realloc(text, sizeBuffer*sizeof(char));
		}
		text[size++] = c;
	}
	*length = size;
	if(size>=sizeBuffer) {
		sizeBuffer ++;
		text = (char*) realloc(text, sizeBuffer*sizeof(char));
	}
	text[size++] = '\0';
	text = (char*) realloc(text, size*sizeof(char));
	if(*alphabet == NULL) {
		*alphabet = getAlphabet(text, NULL);
		fprintf(stderr, "No alphabet provided. Get it from the text.\n");
	}
	symb = toSymbolSequence(text, *length,*alphabet);
	free((void*)text);
	return symb;
}

/*return the symbols read from binary file f*/
TypeSymbol *readFromBinaryFile(FILE *f, int *length) {
	TypeSymbol *symb;
	unsigned char *binary_data;
	int i;
	*length = MIN(*length, getBinaryFileSize(f));
	binary_data = readBinary(f, *length);
	symb = (TypeSymbol*) malloc(*length*sizeof(TypeSymbol));
	for(i=0; i<*length; i++)
		symb[i] = (TypeSymbol) binary_data[i];
	free((void*)binary_data);
	return symb;
}

/*return the alphabet of chains text and motif (no problem if one of them is NULL)*/
char *getAlphabet(char *text, char *motif) {
	int i, size = 0;
	char present[256], *alphabet;
	for(i=0; i<256; i++)
		present[i] = 0;
	if(text != NULL)
		for(i=0; text[i] != '\0'; i++)
			present[text[i]] = 1;
	if(motif != NULL)
		for(i=0; motif[i] != '\0'; i++)
			present[motif[i]] = 1;
	for(i=0; i<256; i++)
		if(present[i])
			size++;
	alphabet = (char*) malloc((size+1)*sizeof(char));
	size = 0;
	for(i=0; i<256; i++)
		if(present[i])
			alphabet[size++] = i;
	alphabet[size++] = '\0';
	return alphabet;
}

/*return the symbol translation of text*/
TypeSymbol *toSymbolSequence(char *text, int length, char *alphabet) {
	int i, cardinal;
	TypeSymbol code[256], *seq;
	for(i=0; i<256; i++)
		code[i] = 255;
	for(cardinal=0; alphabet[cardinal]!='\0'; cardinal++)
		code[alphabet[cardinal]] = cardinal;
	seq = (TypeSymbol*) malloc(length*sizeof(TypeSymbol));
	for(i=0; i<length; i++)
		if(code[text[i]]<cardinal)
			seq[i] = code[text[i]];
		else {
			fprintf(stderr, "Bad character in text:\nCharacter '%c' at position %d replaces by '%c'.\n", text[i], i, alphabet[0]);
			seq[i] = 0;
		}
	return seq;
}

/*print the symbol sequence text*/
void fprintSymbolSequence(FILE *f, TypeSymbol *text, int length, char *alphabet) {
	int i;
	if(alphabet != NULL)
		for(i=0; i<length; i++)
			fprintf(f, "%c", alphabet[text[i]]);
	else
		for(i=0; i<length; i++)
			fprintf(f, "%lu ", (unsigned long) text[i]);
}

/*return the char translation of text*/
char *toCharSequence(TypeSymbol *text, int length, char *alphabet) {
	int i;
	char *seq;
	seq = (char*) malloc((length+1)*sizeof(char));
	for(i=0; i<length; i++)
		if(seq[i]>=0 && seq[i]<strlen(alphabet))
			seq[i] = alphabet[text[i]];
		else
			seq[i] = '?';
	seq[i] = '\0';
	return seq;
}

/*print hexadecimal*/
void fprintSymbolSequenceHexa(FILE *f, TypeSymbol *text, int length, int size) {
	int i;
	char format[10];
	if(length == 0)
		return;
	sprintf(format, "%%0%dx", size);
	fprintf(f, format, (unsigned long) text[0]);
	sprintf(format, " %%0%dx", size);
	for(i=1; i<length; i++)
		fprintf(f, format, (unsigned long) text[i]);
}

/*sprint hexadecimal dest have to be allocated*/
char *sprintSymbolSequenceHexa(char *dest, TypeSymbol *text, int length, int size) {
	int i;
	char format[10];
	if(length == 0) {
		dest[0]= '\0';
		return dest;
	}
	sprintf(format, "%%0%dx", size);
	dest += sprintf(dest, format, (unsigned long) text[0]);
	sprintf(format, " %%0%dx", size);
	for(i=1; i<length; i++)
		dest += sprintf(dest, format, (unsigned long) text[i]);
	return dest;
}

/*fill text with a random text of length size under model*/
void fillRandomPos(TypeSymbol *text, int size, TypeSymbol *textL, int sizeL) {
	int i, pos = (((long)rand())*((long)(sizeL-size)))/((long)RAND_MAX);
	for(i=0; i<size; i++)
		text[i] = textL[pos+i];
}

