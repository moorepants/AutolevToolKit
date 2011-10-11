/* The name of this program is PendulumDynamics.c */
/* Created by Autolev 4.1 on Tue Oct 11 10:20:08 2011 */

#include <ctype.h> 
#include <math.h>  
#include <stdarg.h>
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#define  _NAN      9.99999999999999E+305

void     eqns1     (double T, double VAR[], double VARp[], char boundary);
void     output    (FILE *Fptr[], double T);
void     readf     (FILE *Fp, double *next, ...);
void     pgets     (FILE *Fp, double *x);
void     writef    (FILE *Fp, char format[], ...);
int      kutta     (void(*eqns)(double, double*, double*, char),
                    int numy, double y[], double *t, double integstp,
                    double abserr, double relerr, char com);

double   g,i,l,m;
double   force;
double   omega,theta;
double   omegap,thetap,torque,k,longoutput,p;
double   Pi,DEGtoRAD,RADtoDEG,z[18],A[2][2],B[2][2],C[5][2],D[5][2],Encode[4];

/* ................................ MAIN ............................. */
int      main         (void)
{
FILE     *Fptr[2];
int      iprint, printint, iloop;
double   T, TINITIAL,TFINAL,INTEGSTP,abserr,relerr,_printint;
double   VAR[2];

/* Open input and output files */
for(iloop=0;  iloop<=1;  iloop++)
  {
  char fileName[256];
  if( !iloop ) strcpy(fileName, "PendulumDynamics.in");
  else sprintf(fileName, "PendulumDynamics.%d", iloop);
  if( (Fptr[iloop] = fopen(fileName, iloop ? "w" : "r")) == NULL)
    {printf("Error: unable to open file %s\n", fileName);  exit(0);}
  }
 
/* Read top of input file */
for(iloop=0;  iloop<6;  iloop++) pgets(Fptr[0],NULL);

/* Read values of constants from input file */
readf(Fptr[0],&g,&i,&l,&m,NULL);

/* Read the initial value of each variable from input file */
readf(Fptr[0],&omega,&theta,NULL);

/* Read integration parameters from input file */
readf(Fptr[0],&TINITIAL,&TFINAL,&INTEGSTP,&_printint,&abserr,&relerr,NULL);
printint = (int)_printint;

/* Write heading(s) to output file(s) */
fprintf(stdout,  "%% FILE: PendulumDynamics.1\n%%\n");
fprintf(stdout,  "%%     omega          theta            k              p         longoutput\n"
                 "%% (radian/sec)     (radian)       (joules)       (joules)        (UNITS)\n\n" );
fprintf(Fptr[1], "%% FILE: PendulumDynamics.1\n%%\n");
fprintf(Fptr[1], "%%     omega          theta            k              p         longoutput\n"
                 "%% (radian/sec)     (radian)       (joules)       (joules)        (UNITS)\n\n" );

/* Unit conversions */
  Pi       = 3.141592653589793;
  DEGtoRAD = Pi/180.0;
  RADtoDEG = 180.0/Pi;

/* Evaluate constants */
  z[8] = g*m;
  z[15] = g*l*m;

/* Initialize time, print counter, variables array for integrator */
  T      = TINITIAL;
  iprint = 0;
  VAR[0] = omega;
  VAR[1] = theta;

/* Initalize numerical integrator, with call to eqns1 at T=TINITIAL */
kutta(eqns1, 2, VAR, &T, INTEGSTP, abserr, relerr, 0);

/* Numerically integrate; print results */
while(1)
  {
  if( TFINAL>=TINITIAL && T+.01*INTEGSTP>=TFINAL) iprint=-7;
  if( TFINAL<=TINITIAL && T+.01*INTEGSTP<=TFINAL) iprint=-7;
  if( iprint <= 0 )
    {
    output(Fptr, T);
    if( iprint == -7 ) break;
    iprint = printint;
    }
  if( !kutta(eqns1, 2, VAR, &T, INTEGSTP, abserr, relerr, 2) )
    {
    output(Fptr, T);
    for(iloop=0;  iloop<=1;  iloop++)
      fputs( "\n\tError: Numerical integration failed to converge\n",
              iloop ? Fptr[iloop]:stdout);
    break;
    }
  iprint--;
  }

/* Inform user of input and output filename(s) */
puts( "\n Input is in the file PendulumDynamics.in" );
puts( "\n Output is in the file PendulumDynamics.1" );
puts( "\n The output quantities and associated files are listed in file "
      "PendulumDynamics.dir\n" );
return 0;
}


/* ................................ EQNS1 ............................. */
void     eqns1        (double T, double VAR[], double VARp[], char boundary)
{

/* Update variables after integration step */
  omega = VAR[0];
  theta = VAR[1];

  thetap = omega;
  z[1] = cos(theta);
  z[2] = sin(theta);
  z[3] = pow(z[1],2) + pow(z[2],2);
  z[4] = l*z[3];
  z[10] = i*z[3];
  z[11] = z[3]*z[10] + 0.25*m*pow(z[4],2);

/* Quantities which were specified */
  torque = 10*sin(0.5235987755982988+6.283185307179586*T);

/* Quantities to be specified */
  force = 0;

  z[9] = torque*z[3] - 0.5*z[4]*(force+z[8]*z[2]);
  z[12] = z[9]/z[11];
  omegap = z[12];

/* Update derivative array prior to integration step */
  VARp[0] = omegap;
  VARp[1] = thetap;

}


/* ................................ OUTPUT ............................. */
void     output (FILE *Fptr[], double T)
{
int      i1;

/* Evaluate output quantities */
  k = 0.125*(m*pow(z[4],2)+4*i*pow(z[3],2))*pow(omega,2);
  p = -0.5*g*l*m*z[1];
  longoutput = 2 + 2*theta + p + z[1] + omega + 0.125*(m*pow(z[4],2)+4*i*pow(
  z[3],2))*pow(omega,2);
  z[13] = z[8]*z[1]*z[4]/z[11];
  z[14] = (m*pow(z[4],2)+4*i*pow(z[3],2))*omega;
  z[16] = z[15]*z[2];
  z[17] = 1 + 0.25*(m*pow(z[4],2)+4*i*pow(z[3],2))*omega;

/* Write output to screen and to output file(s) */
  writef(stdout, " %- 14.6E", omega,theta,k,p,longoutput,_NAN);
  writef(Fptr[1]," %- 14.6E", omega,theta,k,p,longoutput,_NAN);
  A[0][0] = 0;
  A[0][1] = -0.5*z[13];
  A[1][0] = 1;
  A[1][1] = 0;
  B[0][0] = z[3]/z[11];
  B[0][1] = -0.5*z[4]/z[11];
  B[1][0] = 0;
  B[1][1] = 0;
  C[0][0] = 1;
  C[0][1] = 0;
  C[1][0] = 0;
  C[1][1] = 1;
  C[2][0] = 0.25*z[14];
  C[2][1] = 0;
  C[3][0] = 0;
  C[3][1] = 0.5*z[16];
  C[4][0] = z[17];
  C[4][1] = 2 + 0.5*z[16] - sin(theta);
  D[0][0] = 0;
  D[0][1] = 0;
  D[1][0] = 0;
  D[1][1] = 0;
  D[2][0] = 0;
  D[2][1] = 0;
  D[3][0] = 0;
  D[3][1] = 0;
  D[4][0] = 0;
  D[4][1] = 0;

  Encode[0] = 0.0;
  Encode[1] = 0.0;
  Encode[2] = 0.0;
  Encode[3] = 0.0;

}


/*................................... READF ................................*/
void  readf( FILE *Fp, double *next, ... )
{
va_list args;                                    /* Variable argument list  */
for( va_start(args,next);  next;  next=va_arg(args,double *) )
  pgets(Fp,next);
va_end(args);                          /* Help function make normal return  */
pgets(Fp,NULL);                        /* Always get a newline at the end   */
}


/*................................... PGETS ................................*/
void  pgets( FILE *Fp, double *x )
{
static long  lineNumber = 0;
char         line[256];

lineNumber++;
if( !fgets(line,256,Fp) )
{
   printf("\n Unable to read line %ld of input file."
          "\n Premature end of file found while reading input file\n", lineNumber);
   exit(0);
}
if( !x ) return;
if( strlen(line) >= 60 )
{
   char *endOfNumber;
   *x = strtod(line+59,&endOfNumber);
   while( isspace(*endOfNumber) )  endOfNumber++;
   if( !*endOfNumber ) return;
}
printf("\n An error occured while reading line %ld of the input file."
       "\n The program was expecting to find one double precision number"
       "\n beginning with the 60th character in the following line:\n\n%s\n",
           lineNumber, line );
exit(0);
}


/* .................................. WRITEF .............................. */
void  writef( FILE *Fp, char format[], ... )
{
va_list  args;                                   /* Variable argument list  */
double   next;                                   /* Current place in list   */

va_start(args,format);                           /* args start after format */
while( (next=va_arg(args,double)) != _NAN )  fprintf(Fp, format, next);
va_end(args);                                    /* End of variable list    */
fprintf(Fp, "\n");                               /* End with newline        */
}


/*****************************************************************************
C**                                                                         **
C** PURPOSE  Solves a set of first order ordinary differential equations    **
C**          of the form dy(i)/dt = F(t,y(1), ..., y(numeqns) (i = 1,       **
C**          ..., numeqns)                                                  **
C**                                                                         **
C** INPUT                                                                   **
C**    eqns: Subroutine that evaluates dy(i)/dt (i = 1, ..., numeqns), the  **
C**          first derivatives of y(1), ..., y(numeqns) with respect to t   **
C**                                                                         **
C** numeqns: The number of differential equations to be solved              **
C**                                                                         **
C**       y: One-dimensional array whose elements are y(1), ..., y(numeqns) **
C**                                                                         **
C**       t: Independent variable                                           **
C**                                                                         **
C** integstp: Maximum integration stepsize                                  **
C**                                                                         **
C**  abserr: Allowable absolute error in y(i)  (i=1, ..., numeqns)          **
C**                                                                         **
C**  relerr: Allowable relative error in y(i)  (i=1, ..., numeqns)          **
C**                                                                         **
C**     com: When com = 2, the Kutta-Merson algorithm (L. Fox, Numerical    **
C**          Solutions of Ordinary and Partial Differential Equations,      **
C**          Palo Alto: Addison-Wesley, 1962, pp. 24-25) is employed to     **
C**          perform the numerical solution of the differential equations.  **
C**          Accordingly, dy(i)/dt (i = 1, ..., numeqns) are evaluated at   **
C**          every integration boundary, including those at Tinitial,       **
C**          Tfinal, and ones created when integstp is halved to satisfy    **
C**          the requirements imposed by abserr and relerr.  Integration    **
C**          is self-starting at each boundary, and the occurrence, at      **
C**          boundaries, of discontinuities in derivatives does not lead    **
C**          to failure of the integration process.                         **
C**                                                                         **
C**          When com = 1, a modified form of the Kutta-Merson algorithm    **
C**          is employed.  It is nearly 20% faster than the one used when   **
C**          com = 2 because no recalculation of derivatives at inte-       **
C**          gration boundaries between Tinitial and Tfinal takes place.    **
C**          Integration is self-starting at Tinitial and Tfinal only.      **
C**          Integration may fail if any of dy(i)/dt (i = 1, ..., numeqns)  **
C**          is discontinuous between Tinitial and Tfinal.                  **
C**                                                                         **
C**          When com = 0, the function eqns is called and dy(i)/dt         **
C**          (i = 1, ..., numeqns) are evaluated, but no integration        **
C**          is performed.                                                  **
C**                                                                         **
C** OUTPUT                                                                  **
C**          The value of t+integstp is returned in t, and the values of    **
C**          y(i) at t+integstp are returned in y.                          **
C**                                                                         **
C** SOURCE                                                                  **
C**          Copyright 1995 by Paul C. Mitiguy, Thomas R. Kane, David A.    **
C**          Levinson, and David B. Schaechter.  Permission is granted      **
C**          to copy, modify, and distribute this subroutine, provided      **
C**          that this copyright notice appear.                             **
C**                                                                         **
C****************************************************************************/
int      kutta        ( void(*eqns)(double, double*, double*, char),
                        int numy, double y[], double *t,
                        double integstp, double abserr, double relerr, char com )
{
static double  f0[100], f1[100], f2[100], y1[100], y2[100];
static int     numcuts = 20;                     /* Max # cuts of integstp  */
static double  hc = 0;                           /* Last value of stepsize  */
char           entry = 1;                        /* Just entered routine    */
char           stepdouble;                       /* Double the stepsize     */
double         tfinal = *t + integstp;           /* Time at end of full step*/
double         tt, h;
int            i;

if(numy >= 100) {printf("\nERROR: INCREASE THE STATIC DOUBLE ARRAY SIZE IN kutta.\n" ); return 0;}
if( !com ) { (*eqns)(*t,y,f0,1);  return 1;}     /* Fill array f0 and return*/
if( numy == 0)  { hc = integstp;  return 1;}     /* Check for initial entry */
if( integstp == 0)  return 0;                    /* Cannot integrate forward*/
if( hc*integstp < 0 ) hc = -hc;                  /* Integrate backward      */
else if( hc == 0 )    hc = integstp;             /* Maybe initial entry     */
h  = hc;                                         /* Current stepsize        */
tt = *t + h;                                     /* Terminal time this step */
*t = tfinal;                                     /* Return updated t value  */

beginning:
while( tt+h != tt )                              /* Check round-off problems*/
  {
  double h2 = h * 0.5;                                     /* Half    of h  */
  double h3 = h / 3.0;                                     /* Third   of h  */
  double h6 = h / 6.0;                                     /* Sixth   of h  */
  double h8 = h * 0.125;                                   /* Eighth  of h  */
  if( com==2 || entry)
 {(*eqns)( tt-h,     y, f0, 1 );   entry=0; }              /* Boundary here */
  for( i=0;  i<numy;  i++)  y1[i] = y[i] + h3*f0[i];
  (*eqns)( tt-2*h3, y1, f1, 0 );
  for( i=0;  i<numy;  i++)  y1[i] = y[i] + h6*(f0[i] + f1[i]);
  (*eqns)( tt-2*h3, y1, f1, 0 );
  for( i=0;  i<numy;  i++)  y1[i] = y[i] + h8*(f0[i] + 3*f1[i]);
  (*eqns)( tt-h2,   y1, f2, 0 );
  for( i=0;  i<numy;  i++)  y1[i] = y[i] + h2*(f0[i] - 3*f1[i] + 4*f2[i]);
  (*eqns)( tt,      y1, f1, 0 );
  for( i=0;  i<numy;  i++)  y2[i] = y[i] + h6*(f0[i] + 4*f2[i] + f1[i]);
  stepdouble = 1;                                /* Assume need to double   */
  for(i=0;  i<numy;  i++)                        /* Check all equations     */
    {
    double error = fabs( y1[i] - y2[i] ) * 0.2;  /* Error in estimate       */
    double test  = fabs( y1[i] ) * relerr;       /* For relative error      */
    if( error >= test  &&  error >= abserr )     /* Error criterion not met?*/
      {                                          /* Restart w/ half stepsize*/
      hc = h = h2;                               /* Halve the stepsize      */
      tt -= h2;                                  /* Change terminal time    */
      if( numcuts-- > 0 )  goto beginning;       /* Back to beginning       */
      printf("\n THE STEPSIZE HAS BEEN HALVED TOO MANY TIMES; T = %15.8E"
             "\n NUMERICAL INTEGRATION FAILED TO CONVERGE.\n", *t=tt-h );
      (*eqns)(*t,y,f0,0);                        /* Fill for error display  */
      return 0;
      }
    if( stepdouble && 64*error>test && 64*error>abserr ) stepdouble = 0;
    }
  for(i=0;  i<numy;  i++)  y[i] = y2[i];         /* Update y for next step  */
  if( stepdouble && fabs(h+h)<=fabs(integstp) && fabs(tt+h+h)<=fabs(tfinal) )
    {hc=(h+=h);  numcuts++;}                     /* Double the stepsize     */
  if( tt == tfinal )                             /* End of integration      */
    { (*eqns)(tfinal,y,f0,2);  return 1;}        /* Derivatives at tfinal   */
  tt += h;                                       /* Update terminal time    */
  /*** If next jump puts tt past or very close to tfinal, adjust h and tt ***/
  if( (h>0 && tt>tfinal-0.1*h) || (h<0 && tt<tfinal-0.1*h) )
    { h = tfinal-(tt-h);  tt = tfinal; }         /* Adjust for last jump    */
  if( com == 1 )                                 /* Approx this derivative  */
    for(i=0;  i<numy;  i++)  f0[i] = f1[i];
  if( com == 3 )                                 /* Calc derivative once    */
    (*eqns)( tt-h, y, f0, 0 );
 }
printf("\nTHE STEPSIZE OF %15.7E IS TOO SMALL RELATIVE TO THE TERMINAL TIME"
       "\nOF %15.7E.  INTEGRATION HALTED BECAUSE OF NUMERICAL ROUND-OFF.\n"
       "\nTHE STEPSIZE MAY HAVE BEEN CUT TOO MANY TIMES.\n\n", h, *t=tt);
(*eqns)(*t,y,f0,0);                              /* Fill for error display  */
return 0;
}



