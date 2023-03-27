
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

void R_init_markovchain(DllInfo* info) {
    R_registerRoutines(info, NULL, NULL, NULL, NULL);
    R_useDynamicSymbols(info, TRUE);
}

/*#define PI 3.141592654*/
#define TOL 1e-20


/* function for binary outcomes*/
void cindexLog(double *prob, int *fstatus, int *n, int *npair,double *cindex){
    int i,j;
    int cc=0,nn=0;
    for(i=0;i< *n;i++)
    for(j=i+1;j< *n;j++){
        if((fstatus[i] == 1) && (fstatus[j] == 0))
        { nn += 1;
        if(prob[i] > prob[j]) cc += 1;
        else {
            if(prob[i] < prob[j]) cc += 0;
            else  cc += 0.5;
        }
        }
        else if((fstatus[i]==0) && (fstatus[j] == 1)){
        nn += 1;
        if(prob[i] < prob[j]) cc  += 1;
        else {
            if( prob[i] > prob[j]) cc  += 0;
            else cc += 0.5;
        }
        }
        else {
        cc += 0;
        nn += 0;
        }
    }
    npair[0] = nn;
    npair[1] = cc;
    *cindex = (double) cc/(double) nn;

}

/* function for ordinary survival outcomes*/

void cindexSurv(double *prob, int *fstatus, double *ftime,int *n, int *npair,double *cindex){
    int i,j;
    int cc=0,nn=0;
    for(i=0;i< *n;i++)
    for(j=i+1;j< *n;j++){
        if(ftime[i] < ftime[j]){
        if((fstatus[i]==1) && (fstatus[j] == 1))
        {   nn +=1;
            if(prob[i] > prob[j]) cc  += 1;
            else {
            if(prob[i] < prob[j]) cc  += 0;
            else cc += 0.5;
            }

        }
        else if((fstatus[i]==1) && (fstatus[j] == 0))
        { nn += 1;
            if(prob[i] > prob[j]) cc += 1;
            else {
            if(prob[i] < prob[j]) cc += 0;
            else cc += 0.5;
            }
        }
        else {
            cc +=0;
            nn += 0;
        }
        }
/* else here means ftime[i] == ftime[j] */
        else {
        if((fstatus[i]==0) && (fstatus[j] == 1))
        { nn += 1;
            if(prob[i] < prob[j]) cc += 1;
            else {
            if(prob[i] > prob[j]) cc += 0;
            else cc += 0.5;
            }
        }
        else if((fstatus[i]==1) && (fstatus[j] == 0))
        { nn += 1;
            if(prob[i] > prob[j]) cc += 1;
            else {
            if(prob[i] < prob[j]) cc += 0;
            else cc += 0.5;
            }
        }
        else {
            cc +=0;
            nn += 0;
        }

        }
    }
    npair[0] = nn;
    npair[1] = cc;
    *cindex = (double) cc/(double) nn;

}

/* function for survival outcomes  with competing risks*/

void cindexCrr(double *prob, int *fstatus, double *ftime,int *n, int *npair,double *cindex){
    int i,j;
    int cc=0,nn=0;
    for(i=0;i< *n;i++)
    for(j=i+1;j< *n;j++){
        if(ftime[i] < ftime[j]){
        if((fstatus[i]==1) && (fstatus[j] == 1))
        {   nn +=1;
            if(prob[i] > prob[j]) cc  += 1;
            else {
            if(prob[i] < prob[j]) cc  += 0;
            else cc += 0.5;
            }
        }
        else if((fstatus[i]==1) && (fstatus[j] == 2))
        { nn += 1;
            if(prob[i] > prob[j]) cc += 1;
            else {
            if(prob[i] < prob[j]) cc  += 0;
            else cc += 0.5;
            }
        }
        else if((fstatus[i]==1) && (fstatus[j] == 0))
        { nn += 1;
            if(prob[i] > prob[j]) cc += 1;
            else {
            if(prob[i] < prob[j]) cc  += 0;
            else cc += 0.5;
            }
        }
        else if((fstatus[i]==2) && (fstatus[j] == 1))
        { nn += 1;
            if(prob[i] < prob[j]) cc += 1;
            else {
            if(prob[i] > prob[j]) cc  += 0;
            else cc += 0.5;
            }
        }
        else {
            cc +=0;
            nn += 0;
        }
        }
/* else here means ftime[i] == ftime[j] since ftime was sorted in ecending order */
        else {
        if((fstatus[i]==1) && (fstatus[j] == 2))
        { nn += 1;
            if(prob[i] > prob[j]) cc += 1;
                else {
            if(prob[i] < prob[j]) cc  += 0;
            else cc += 0.5;
            }
                }
        else if((fstatus[i]==1) && (fstatus[j] == 0))
        { nn += 1;
            if(prob[i] > prob[j]) cc += 1;
            else {
            if(prob[i] < prob[j]) cc  += 0;
            else cc += 0.5;
            }
        }
        else if((fstatus[i]==2) && (fstatus[j] == 1))
        { nn += 1;
            if(prob[i] < prob[j]) cc += 1;
            else {
            if(prob[i] > prob[j]) cc  += 0;
            else cc += 0.5;
            }
        }
        else if((fstatus[i]==0) && (fstatus[j] == 1))
        { nn += 1;
            if(prob[i] < prob[j]) cc += 1;
            else {
            if(prob[i] > prob[j]) cc  += 0;
            else cc += 0.5;
            }
        }
        else {
            cc +=0;
            nn += 0;
        }

        }
    }
    npair[0] = nn;
    npair[1] = cc;
    *cindex = (double) cc/(double) nn;

}

/*
int main(){
    double x;
     printf("Please enter a numeric\n");
    scanf("%lf",&x);
    printf("x = %lf\n",x);
}
*/
