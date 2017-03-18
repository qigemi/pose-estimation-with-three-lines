/* libbairstow.h
 * rebX
 *
*/

#ifndef _BAIRSTOWLIB
#define _BAIRSTOWLIB 1

#include <math.h>

#define DEGREE 15 /* Highest degree of polynomial */
#define EPSILON 0.00001

void getsignal(char *);
char getsign(float);
void getroots(float,float,float);
float getpoly(void);
void getcoeff(int,float []);
void putparam(int,float [],float,int);
int bairstow(int,float [],float,int);
#endif
