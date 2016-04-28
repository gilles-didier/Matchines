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




#ifndef AsymptoticSpeedF
#define AsymptoticSpeedF

#include "MatchingMachine.h"
#include "Bernoulli.h"

/*return a table C where C[st] is the communicating class of st if C[st]>0 and C[st]=0 if st is transient*/
TypeState *getCommunicatingClasses(TypeMatchingMachine *matchine, TypeBernoulli *bernoulli);
/*return the asymptotic speed of matchine wrt model bernoulli*/
double getAsymptoticSpeed(TypeMatchingMachine *matchine, TypeBernoulli *bernoulli);
/*return the asymptotic speed of a given communicating class of matchine wrt model bernoulli*/
double getAsymptoticSpeedCommunicatingClass(TypeState *comm, TypeState sizeComm, TypeState *index, TypeMatchingMachine *matchine, TypeBernoulli *bernoulli);

#endif
