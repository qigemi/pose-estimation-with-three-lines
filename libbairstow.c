
/* libbairstow.c

*/

#include <stdio.h>
#include <math.h>

#ifndef _BAIRSTOWLIB
#include "libbairstow.h"
#endif

void getsignal(char *msg) {
    printf("%s",msg);
}

char getsign(float x) {
    if (x>0)
        return '+';
    return  ' ';
}

void getroots(float u,float v,float epsilon) {
    float r1,r2,discr;

    discr=(u*u)-(4*v);
    u=-u;
    /* Make adjustments with respect to epsilon */
    if ((discr-epsilon)<0)
        discr=0;
    if (discr==0) {
        r1=u/2;
        printf("Found identical real roots.\n%f\n%f\n",r1,r1);
    } else if (discr>0) {
        r1=(u+sqrt(discr))/2;
        r2=(u-sqrt(discr))/2;
        printf("Found distinct real roots.\n%f\n%f\n",r1,r2);
    } else {
        u/=2;
        discr=sqrt(-discr)/2;
        printf("Found imaginary roots.\n%f+%fi\n%f-%fi\n"
               ,u,discr,u,discr);
    }
}

float getpoly(void) {
    int n;
    do {
        printf("Enter the degree of the polynomial, n(2<n<%d): ",DEGREE);
        scanf("%d",&n);
    } while (n>DEGREE || n<2);
    return n;
}

void getcoeff(int n,float x[]) {
    int i,decr = n;

    printf("For the coefficients of x, \n");
    printf("x[3] is the coefficient of the x raised to the third power\n");
    printf("x[0] is for the constant\n");

    for (i=0; i<=n; i++) {
        printf("Enter values for x[%d]: ",decr);
        scanf("%f",&x[i]);
        decr--;
    }
}

void putparam(int n,float x[],float epsilon,int maxiter) {
    int i;

    printf("Finding the roots of the polynomial\n");
    for (i=0; i<=(n-2); i++)
        printf("%c%.2fx^%d",getsign(x[i]),x[i],n-i);
    printf("%c%.2fx",getsign(x[n-1]),x[n-1]);
    printf("%c%.2f",getsign(x[n]),x[n]);
    printf(" with epsilon=%f up to %d iterations.\n",epsilon,maxiter);
}

int bairstow(int n,float x[],float epsilon,int maxiter) {
    int i,iter;
    int isclosetozero; /* Boolean */
    float b[DEGREE+1],c[DEGREE+1];
    float u,v,deltau,deltav,denom;

    while (n>=3) {
        /* Step 1 */
        iter=0;
        u=v=0;
        isclosetozero=0;
        do {
            /* Step 2 */
            b[0]=1;
            b[1]=x[1]-u;
            for (i=2; i<=n; i++)
                b[i]=x[i]-(b[i-1]*u)-(b[i-2]*v);
            /* Step 3 */
            c[0]=1;
            c[1]=b[1]-u;
            for (i=2; i<n; i++)
                c[i]=b[i]-(c[i-1]*u)-(c[i-2]*v);
            /* Step 4 */
            denom=(c[n-1]*c[n-3])-(c[n-2]*c[n-2]);
            if (denom==0) {
                getsignal("Error: denom is zero!");
                return 1;
            }
            deltau=((b[n]*c[n-3])-(b[n-1]*c[n-2]))/denom;
            deltav=((c[n-1]*b[n-1])-(c[n-2]*b[n]))/denom;
            /* Step 5 */
            u+=deltau;
            v+=deltav;
            if (iter<maxiter) {
                printf("Iteration number %d (u=%f,v=%f):\n",iter+1,u,v);
                for (i=0; i<=n; i++)
                    printf("x[%d] ",i);
                printf("\n");
                for (i=0; i<=n; i++)
                    printf("%.2f ",x[i]);
                printf("\n");
                for (i=0; i<=n; i++)
                    printf("b[%d] ",i);
                printf("\n");
                for (i=0; i<=n; i++)
                    printf("%.2f ",b[i]);
                printf("\n");
                for (i=0; i<n; i++)
                    printf("c[%d] ",i);
                printf("\n");
                for (i=0; i<n; i++)
                    printf("%.2f ",c[i]);
                iter++;
                getsignal("\nContinuing to next iteration...\n");
            } else {
                iter=0;
                u=x[n-1]/x[n-2];
                v=x[n]/x[n-2];
                getsignal("Limit for number of iterations is reached.\n"
                          "Resetting iteration using different values of u and v.\n");
                continue;
            }
            isclosetozero=(fabs(deltau) + fabs(deltav))<epsilon;
        } while (!isclosetozero);
        getroots(u,v,epsilon);
        getsignal("Preparing to extract roots from the next factor.\n");
        n-=2;
        for (i=0; i<=n; i++)
            x[i]=b[i];
    }
    if (n%2==1)
        printf("Found a real root.\n%f\n",-b[1]);
    else
        getroots(x[1],x[2],epsilon);
    return 0;
}
