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
#include <math.h>

#include "TableUtils.h"


/*replace new lines and tabs of text with spaces*/
void fixString(char *text) {
	int i;
	for(i=0; text[i]!='\0'; i++)
		if(text[i] == '\n' || text[i] == '\r' || text[i] == '\t')
			text[i] = ' ';
}

/*print table*/
void fprintTable(FILE *fo, char **column, char **row, double **val, int nrow, int ncol, int rgrp, int doMax, char type, int sgn) {
	switch(type) {
		case 'l':
			fprintTableLatex(fo, column, row, val, nrow, ncol, rgrp, doMax, sgn);
			break;
		case 'c':
		default:
			fprintTableCSV(fo, column, row, val, nrow, ncol, rgrp, doMax, sgn);
			break;
	}
}

/*print table in Latex format*/
void fprintTableLatex(FILE *fo, char **column, char **row, double **val, int nrow, int ncol, int rgrp, int doMax, int sgn) {
	int i, j, k, n;
	n = nrow/rgrp;
	double prec = pow(10., (double) -sgn);
	char form[50];
	fprintf(fo, "%%\\usepackage{adjustbox}\n%%\\usepackage{array}\n\n\\newcolumntype{R}[2]{>{\\adjustbox{angle=#1,lap=\\width-(#2)}\\bgroup}l<{\\egroup}}\n\\newcommand*\\rot{\\multicolumn{1}{R{45}{1em}}}\n");
	fprintf(fo, "\\begin{tabular}{l");
	if(column != NULL) {
		for(i=0; i<ncol; i++)
			fprintf(fo, "r");
		fprintf(fo, "}\n~");
		for(i=0; i<ncol; i++)
			if(column[i] != NULL) {
				fixString(column[i]);
				fprintf(fo, "& \\rot{%s} ", column[i]);
			}else
				fprintf(fo, "& ");
		fprintf(fo, "\\\\\n\\hline\n");
	}
	for(i=0; i<n; i++) {
		double max;
		if(row != NULL && row[i] != NULL) {
			fixString(row[i]);
			fprintf(fo, "\\texttt{%s} ", row[i]);
		}
		if(doMax) {
			max = 0.;
			for(j=0; j<rgrp; j++)
				for(k=0; k<ncol; k++)
					if(val[i*rgrp+j][k] != TABLE_NO_VALUE && val[i*rgrp+j][k]>=max)
						max = val[i*rgrp+j][k];
			max = prec*(round(max/prec));
			for(j=0; j<rgrp; j++) {
				for(k=0; k<ncol; k++)
					if(val[i*rgrp+j][k] != TABLE_NO_VALUE && prec*(round(val[i*rgrp+j][k]/prec)) == max) {
						sprintf(form, "& \\textbf{%%.%df} ", sgn);
						fprintf(fo, form, val[i*rgrp+j][k]);
					} else {
						if(val[i*rgrp+j][k] != TABLE_NO_VALUE) {
							sprintf(form, "& %%.%df ", sgn);
							fprintf(fo, form, val[i*rgrp+j][k]);
						} else
							fprintf(fo, "& -- ");
					}
				fprintf(fo, "\\\\\n");
			}
		} else {
			for(j=0; j<rgrp; j++) {
				for(k=0; k<ncol; k++)
					if(val[i*rgrp+j][k] != TABLE_NO_VALUE) {
						sprintf(form, "& %%.%df ", sgn);
						fprintf(fo, form, val[i*rgrp+j][k]);
					} else
						fprintf(fo, "& -- ");
				fprintf(fo, "\\\\\n");
			}
		}
		if(rgrp > 1)
			fprintf(fo, "\\hline\n");
	}
	if(rgrp == 1)
		fprintf(fo, "\\hline\n");
	fprintf(fo, "\\end{tabular}\n");
}

/*print table in CSV format*/
void fprintTableCSV(FILE *fo, char **column, char **row, double **val, int nrow, int ncol, int rgrp, int doMax, int sgn) {
	int i, j, k, n;
	char form[50];
	n = nrow/rgrp;
	if(column != NULL) {
		for(j=0; j<ncol; j++)
			if(column[j] != NULL) {
				fixString(column[j]);
				fprintf(fo, "\t%s", column[j]);
			} else
				fprintf(fo, "\t ");
		fprintf(fo, "\n");
	}
	for(i=0; i<n; i++) {
		if(row != NULL && row[i] != NULL) {
			fixString(row[i]);
			fprintf(fo, "%s", row[i]);
		}
		for(j=0; j<rgrp; j++) {
			for(k=0; k<ncol; k++) {
				if(val[i*rgrp+j][k] != TABLE_NO_VALUE) {
					sprintf(form, "\\t%%.%df", sgn);
					fprintf(fo, form, val[i*rgrp+j][k]);
				} else
					fprintf(fo, "\t ");
			}
			fprintf(fo, "\n");
		}
	}
}
