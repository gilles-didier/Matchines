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




#ifndef TableUtilsF
#define TableUtilsF

#include <stdio.h>
#include <limits.h>
#include <float.h>

#define TABLE_NO_VALUE DBL_MAX

/*replace new lines and tabs of text with spaces*/
void fixString(char *text);
/*print table*/
void fprintTable(FILE *fo, char **column, char **row, double **val, int nrow, int ncol, int rgrp, int doMax, char type, int sgn);
void fprintTableLatex(FILE *fo, char **column, char **row, double **val, int nrow, int ncol, int rgrp, int doMax, int sgn);
void fprintTableCSV(FILE *fo, char **column, char **row, double **val, int nrow, int ncol, int rgrp, int doMax, int sgn);

#endif
