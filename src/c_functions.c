#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>



// Takes the covariates and a beta as inputs and returns a matrix where each row
// corresponds to a vertex and contains the covariates multiplied with beta.
SEXP multiply_covariates(SEXP covariates,SEXP beta)
{
  // Variables
  int p,q,T,i,t,r;
  double* matrix;
  SEXP out;
  
  // Find dimensions
  q=LENGTH(beta);
  p=LENGTH(VECTOR_ELT(VECTOR_ELT(covariates,1),0))/q;
  T=LENGTH(VECTOR_ELT(covariates,0));

  // Prepare output
  matrix=malloc(sizeof(double)*T*p);
  
  // Compute output
  for(t=0;t<T;t++)
  {
    for(r=0;r<p;r++)
    {
      matrix[r+t*p]=0;
      for(i=0;i<q;i++)
        matrix[r+t*p]=matrix[r+t*p]+REAL(beta)[i]*REAL(VECTOR_ELT(VECTOR_ELT(covariates,1),t))[r+i*p];
    }
  }
  
  // Transfer Output to R
  out=PROTECT(allocVector(REALSXP,p*T));
  for(i=0;i<p*T;i++)
    REAL(out)[i]=matrix[i];
  
  // Clean Up
  free(matrix);
  UNPROTECT(1);
  
  return(out);
}

// Locate a given time in an ordered time list
int locate_time(SEXP times,double t0)
{
  int L,lower,upper,range,ind;
  
  L=LENGTH(times);
  lower=0;
  upper=L-1;
  range=upper-lower;
  
  while(range>1)
  {
    if(range%2==0)
    {
      // Range is even
      ind=lower+range/2;
    }
    else
    {
      // Range is odd
      ind=lower+(range-1)/2;
    }
    if(REAL(times)[ind]<=t0)
      lower=ind;
    else
      upper=ind;
    
    range=upper-lower;
  }
  
  return(lower);
}

// Reject Events in the baseliner
SEXP baseline_reject(SEXP proposed_times,SEXP total_event_number,SEXP event_numbers,SEXP randomness,SEXP baselines,SEXP times,SEXP M)
{
  // Variables
  int p,r,i,j,k,list_end,id;
  double test_base;
  double* event_list;
  SEXP out;
  
  // Read Data
  p=LENGTH(event_numbers);
  
  // Allocate Memory
  event_list=malloc(sizeof(double)*4*asInteger(total_event_number));
  list_end=0;
  
  // Perform rejection
  i=0;
  id=1;
  for(r=0;r<p;r++)
  {
    for(k=0;k<INTEGER(event_numbers)[r];k++)
    {
      // Find baseline to test
      test_base=REAL(baselines)[r+locate_time(times,REAL(proposed_times)[i])*p];

      // Perform Test
      if(REAL(randomness)[i]<=test_base/REAL(M)[r])
      {
        event_list[list_end+0*asInteger(total_event_number)]=(double)id;
        event_list[list_end+1*asInteger(total_event_number)]=0;
        event_list[list_end+2*asInteger(total_event_number)]=(double)(r+1);
        event_list[list_end+3*asInteger(total_event_number)]=REAL(proposed_times)[i];
        
        id++;
        list_end++;
      }
      i++;
    }
  }
  
  // Transform to R output
  out=PROTECT(allocVector(REALSXP,list_end*4));
  for(i=0;i<list_end;i++)
  {
    for(j=0;j<4;j++)
      REAL(out)[i+j*list_end]=event_list[i+j*asInteger(total_event_number)];
  }
  
  // Clean Up
  free(event_list);
  UNPROTECT(1);
  
  return(out);
}

// Generate uniform random number
double runif(void)
{
  double U;
  U=((double)rand())/((double)RAND_MAX);
  return(U);
}

// Generate a random exponential number
double rexp(double rate)
{
  double E;
  E=-1/rate*log(1-runif());
  return(E);
}

// Add event to event list and extend list when required
double* add_to_event_list(double* event_list,int* list_length,int* list_end,int current,double id,double parent_id,double process,double time)
{
  double* new_event_list;
  int new_length;
  int i,j;
  
  // Check if new memory needs to be allocated
  if(list_end[0]>=list_length[0])
  {
    // We require new memory
    new_length=(list_end[0]-current)*5+list_length[0];
    new_event_list=malloc(sizeof(double)*new_length*4);
    
    // Copy old list to new memory location
    for(i=0;i<list_length[0];i++)
    {
      for(j=0;j<4;j++)
        new_event_list[i+j*new_length]=event_list[i+j*list_length[0]];
    }
    
    // Free old memory and copy location
    free(event_list);
    event_list=new_event_list;
    list_length[0]=new_length;
  }
  
  // Add new line to event list
  event_list[list_end[0]+0*list_length[0]]=id;
  event_list[list_end[0]+1*list_length[0]]=parent_id;
  event_list[list_end[0]+2*list_length[0]]=process;
  event_list[list_end[0]+3*list_length[0]]=time;
  list_end[0]=list_end[0]+1;
  
  return(event_list);
}

// Generate spin-off events
SEXP generate_spin_off_events(SEXP event_list_R,SEXP C_R,SEXP rand_seed,SEXP T_R,SEXP gamma_R,SEXP first_id)
{
  // Variables
  int r,i,j,p;
  int list_end,list_length;
  int source;
  int id;
  double factor,M,tprop;
  double parent_time;
  double T;
  double gamma;
  double current_intensity;
  double* event_list;
  SEXP out;
  
  // Set seed for random number calculation
  srand(asInteger(rand_seed));
  
  // Read Information
  list_end=LENGTH(event_list_R)/4;
  list_length=list_end;
  T=asReal(T_R);
  gamma=asReal(gamma_R);
  id=asInteger(first_id);
  p=(int)sqrt((double)LENGTH(C_R));
  
  // Allocate Memory
  event_list=malloc(sizeof(double)*4*list_length);
  
  // Copy event list
  for(i=0;i<list_length*4;i++)
    event_list[i]=REAL(event_list_R)[i];
  
  // Generate Spin-Off Events
  for(r=0;r<list_end;r++)
  {
    // Check if Event generates Spin-Off Events
    source=(int)event_list[r+2*list_length]-1;
    for(j=0;j<p;j++)
    {
      if(REAL(C_R)[j+source*p]>0)
      {
        // Spin-Off Events have to be generated
        parent_time=event_list[r+3*list_length];
        factor=REAL(C_R)[j+source*p];
        M=factor;
        tprop=rexp(M);
        
        while(tprop+parent_time<=T)
        {
          // Check if Proposal gets accepted
          current_intensity=factor*exp(-gamma*tprop);
          if(runif()<=current_intensity/M)
          {
            // Accept Event and add to event list
            event_list=add_to_event_list(event_list,&list_length,&list_end,r,(double)id,event_list[r],(double)(j+1),tprop+parent_time);
            M=current_intensity;
            id++;
          }
          
          // Compute new proposition
          tprop=tprop+rexp(M);
        }
      }
    }
  }
  
  // Prepare Output
  out=PROTECT(allocVector(REALSXP,list_end*4));
  for(i=0;i<list_end;i++)
    for(j=0;j<4;j++)
      REAL(out)[i+j*list_end]=event_list[i+j*list_length];
  
  // Clean Up
  free(event_list);
  UNPROTECT(1);
  
  return(out);
}

// Compute the intensity integrals
SEXP compute_intensity_integrals(SEXP event_list,SEXP gamma,SEXP times,SEXP baseline,SEXP C_R)
{
  // Variables
  int L,p;
  int i,k,j;
  int ind;
  double factor;
  double addon;
  double* integrals;
  double* total_intensity;
  SEXP out;
  SEXP out_integrals;
  SEXP out_intensities;
  
  // Read Information
  L=LENGTH(event_list)/4;
  p=(int)sqrt((double)LENGTH(C_R));
  
  // Allocate Memory
  integrals=malloc(sizeof(double)*p*L);
  total_intensity=malloc(sizeof(double)*p*L);
  
  // Compute Intensity Integrals and Total Intensity At time of the first event
  // all integrals are zero and the intensity is just the baseline intensity
  ind=locate_time(times,REAL(event_list)[0+3*L]);
  for(i=0;i<p;i++)
  {
    integrals[i]=0;
    total_intensity[i]=REAL(baseline)[i+ind*p];
  }
  
  // At following time points
  for(k=1;k<L;k++)
  {
    // Locate Time in Baseline
    ind=locate_time(times,REAL(event_list)[k+3*L]);
    
    // Compute attenuating factor
    factor=exp(-asReal(gamma)*(REAL(event_list)[k+3*L]-REAL(event_list)[k-1+3*L]));
    
    for(i=0;i<p;i++)
    {
      if(i+1==(int)REAL(event_list)[k-1+2*L])
        addon=factor;
      else
        addon=0;
      
      integrals[i+k*p]=factor*integrals[i+(k-1)*p]+addon;
    }
    for(i=0;i<p;i++)
    {
      total_intensity[i+k*p]=REAL(baseline)[i+ind*p];
      for(j=0;j<p;j++)
        total_intensity[i+k*p]=total_intensity[i+k*p]+REAL(C_R)[i+j*p]*integrals[j+k*p];
    }
  }
  
  // Prepare Output
  out_integrals=PROTECT(allocVector(REALSXP,p*L));
  out_intensities=PROTECT(allocVector(REALSXP,p*L));
  for(i=0;i<p*L;i++)
  {
    REAL(out_integrals)[i]=integrals[i];
    REAL(out_intensities)[i]=total_intensity[i];
  }
  out=PROTECT(allocVector(VECSXP,2));
  SET_VECTOR_ELT(out,0,out_intensities);
  SET_VECTOR_ELT(out,1,out_integrals);
  
  // Clean Up
  free(total_intensity);
  free(integrals);
  UNPROTECT(3);
  
  return(out);
}

// Compute vector v
SEXP compute_vector_v(SEXP event_list,SEXP nu0,SEXP times)
{
  // Variables
  int p,N,L;
  int i,r;
  int tind;
  double* v;
  SEXP out;
  
  // Get information
  N=LENGTH(times);
  p=LENGTH(nu0)/N;
  L=LENGTH(event_list)/4;
  
  // Allocate Memory
  v=malloc(sizeof(double)*p);
  
  // Compute v
  // Initialize with zeros
  for(i=0;i<p;i++)
    v[i]=0;
  
  // Update
  tind=0;
  for(r=0;r<L;r++)
  {
    while(REAL(times)[tind]<=REAL(event_list)[r+3*L])
      tind++;
    tind--;
    
    i=(int)REAL(event_list)[r+2*L]-1;
    v[i]=v[i]+REAL(nu0)[i+tind*p];
  }
  
  // Prepare Output
  out=PROTECT(allocVector(REALSXP,p));
  for(i=0;i<p;i++)
    REAL(out)[i]=v[i];
  
  // Clean Up
  free(v);
  UNPROTECT(1);
  
  return(out);
}

// Compute derivatives of vector v
SEXP dbeta_vector_v(SEXP event_list,SEXP nu0,SEXP covariates)
{
  // Variables
  int p,q,N,L;
  int i,j,k,r;
  int tind;
  SEXP dbetav;
  SEXP total2deriv;
  SEXP out;
  SEXP help;
  
  // Get information
  N=LENGTH(VECTOR_ELT(covariates,0));
  p=LENGTH(nu0)/N;
  L=LENGTH(event_list)/4;
  q=LENGTH(VECTOR_ELT(VECTOR_ELT(covariates,1),0))/p;
  
  // Allocate Memory
  dbetav =PROTECT(allocVector(REALSXP,p*q));
  total2deriv=PROTECT(allocVector(VECSXP,q));
  for(i=0;i<q;i++)
  {
    help=PROTECT(allocVector(REALSXP,p*q));
    SET_VECTOR_ELT(total2deriv,i,help);
  }

  // Compute v
  // Initialize with zeros
  for(i=0;i<p;i++)
  {
    for(j=0;j<q;j++)
    {
      REAL(dbetav)[i+j*p]=0;
      for(k=0;k<q;k++)
        REAL(VECTOR_ELT(total2deriv,k))[i+j*p]=0;
    }
  }
  
  // Update
  tind=0;
  for(r=0;r<L;r++)
  {
    while(REAL(VECTOR_ELT(covariates,0))[tind]<=REAL(event_list)[r+3*L])
      tind++;
    tind--;
    
    i=(int)REAL(event_list)[r+2*L]-1;
    
    for(k=0;k<q;k++)
    {
      REAL(dbetav)[i+k*p] =REAL(dbetav )[i+k*p]+REAL(VECTOR_ELT(VECTOR_ELT(covariates,1),tind))[i+k*p]                                                       *REAL(nu0)[i+tind*p];
      
      for(j=0;j<q;j++)
        REAL(VECTOR_ELT(total2deriv,j))[i+k*p]=REAL(VECTOR_ELT(total2deriv,j))[i+k*p]+REAL(VECTOR_ELT(VECTOR_ELT(covariates,1),tind))[i+k*p]*REAL(VECTOR_ELT(VECTOR_ELT(covariates,1),tind))[i+j*p]*REAL(nu0)[i+tind*p];
    }
  }
  
  // Prepare Output
  out=PROTECT(allocVector(VECSXP,2));
  SET_VECTOR_ELT(out,0,dbetav );
  SET_VECTOR_ELT(out,1,total2deriv);
  
  // Clean Up
  UNPROTECT(3+q);
  
  return(out);
}


double myabs(double x)
{
  if(x>=0)
    return(x);
  else
    return(-x);
}

double myclip(double x)
{
  if(x>=0)
    return(x);
  else
    return(0);
}

SEXP compute_gamma(SEXP p_R,SEXP event_list,SEXP gammaR,SEXP TR)
{
  // Variables
  int i,j,k,r;
  int p,L;
  double gamma,ti,tj,T;
  SEXP Gamma;
  
  // Read Information
  p=asInteger(p_R);
  L=LENGTH(event_list)/4;
  gamma=asReal(gammaR);
  T=asReal(TR);
  
  // Allocate Memory
  Gamma=PROTECT(allocVector(REALSXP,p*p));
  for(i=0;i<p*p;i++)
    REAL(Gamma)[i]=0;
  
  // Compute Gamma
  for(k=0;k<L;k++)
  {
    for(r=0;r<L;r++)
    {
      i=(int)REAL(event_list)[k+2*L]-1;
      j=(int)REAL(event_list)[r+2*L]-1;
      ti=REAL(event_list)[k+3*L];
      tj=REAL(event_list)[r+3*L];
      
      REAL(Gamma)[i+j*p]=REAL(Gamma)[i+j*p]+1/(2*gamma)*(exp(-gamma*myabs(ti-tj))-exp(-2*gamma*(T-(ti+tj)/2)));
    }
  }
  
  // Clean Up
  UNPROTECT(1);
  
  return(Gamma);
}

// Compute derivatives of Gamma
SEXP compute_gamma_deriv(SEXP p_R,SEXP event_list,SEXP gammaR,SEXP TR)
{
  // Variables
  int i,j,k,r;
  int p,L;
  double gamma,ti,tj,T;
  double abs_titj;
  double avg_titj;
  SEXP dGamma;
  SEXP d2Gamma;
  SEXP out;
  
  // Read Information
  p=asInteger(p_R);
  L=LENGTH(event_list)/4;
  gamma=asReal(gammaR);
  T=asReal(TR);
  
  // Allocate Memory
  dGamma =PROTECT(allocVector(REALSXP,p*p));
  d2Gamma=PROTECT(allocVector(REALSXP,p*p));
  for(i=0;i<p*p;i++)
  {
    REAL(dGamma)[i]=0;
    REAL(d2Gamma)[i]=0;
  }
  
  // Compute derivatives of Gamma
  for(k=0;k<L;k++)
  {
    for(r=0;r<L;r++)
    {
      i=(int)REAL(event_list)[k+2*L]-1;
      j=(int)REAL(event_list)[r+2*L]-1;
      ti=REAL(event_list)[k+3*L];
      tj=REAL(event_list)[r+3*L];
      
      abs_titj=myabs(ti-tj);
      avg_titj=(ti+tj)/2;
      
      REAL(dGamma)[i+j*p] =REAL(dGamma)[i+j*p] +exp(-gamma*abs_titj)*(-1/(2*gamma*gamma)-abs_titj/(2*gamma))+exp(-2*gamma*(T-avg_titj))*(1/(2*gamma*gamma)+(T-avg_titj)/gamma);
      REAL(d2Gamma)[i+j*p]=REAL(d2Gamma)[i+j*p]+exp(-gamma*abs_titj)*(abs_titj/(gamma*gamma)+abs_titj*abs_titj/(2*gamma)+1/(gamma*gamma*gamma))+exp(-2*gamma*(T-avg_titj))*(-2*(T-avg_titj)/(gamma*gamma)-2*(T-avg_titj)*(T-avg_titj)/gamma-1/(gamma*gamma*gamma));
    }
  }
  
  // Prepare Outout
  out=PROTECT(allocVector(VECSXP,2));
  SET_VECTOR_ELT(out,0,dGamma);
  SET_VECTOR_ELT(out,1,d2Gamma);
  
  // Clean Up
  UNPROTECT(3);
  
  return(out);
}

SEXP compute_A(SEXP p_R,SEXP event_list,SEXP gammaR,SEXP TR)
{
  // Variables
  int i,j,k,r;
  int p,L;
  double gamma,ti,tj;
  SEXP A;
  
  // Read Information
  p=asInteger(p_R);
  L=LENGTH(event_list)/4;
  gamma=asReal(gammaR);

  // Allocate Memory
  A=PROTECT(allocVector(REALSXP,p*p));
  for(i=0;i<p*p;i++)
    REAL(A)[i]=0;
  
  // Compute A
  for(k=1;k<L;k++)
  {
    for(r=0;r<k;r++)
    {
      i=(int)REAL(event_list)[k+2*L]-1;
      j=(int)REAL(event_list)[r+2*L]-1;
      ti=REAL(event_list)[k+3*L];
      tj=REAL(event_list)[r+3*L];
      
      REAL(A)[i+j*p]=REAL(A)[i+j*p]+exp(-gamma*(ti-tj));
    }
  }
  
  // Clean Up
  UNPROTECT(1);
  
  return(A);
}

// Compute Derivatives of A
SEXP compute_A_deriv(SEXP p_R,SEXP event_list,SEXP gammaR,SEXP TR)
{
  // Variables
  int i,j,k,r;
  int p,L;
  double gamma,ti,tj;
  SEXP dA;
  SEXP d2A;
  SEXP out;
  
  // Read Information
  p=asInteger(p_R);
  L=LENGTH(event_list)/4;
  gamma=asReal(gammaR);

  // Allocate Memory
  dA=PROTECT(allocVector(REALSXP,p*p));
  d2A=PROTECT(allocVector(REALSXP,p*p));
  for(i=0;i<p*p;i++)
  {
    REAL(dA )[i]=0;
    REAL(d2A)[i]=0;
  }
  
  // Compute A
  for(k=1;k<L;k++)
  {
    for(r=0;r<k;r++)
    {
      i=(int)REAL(event_list)[k+2*L]-1;
      j=(int)REAL(event_list)[r+2*L]-1;
      ti=REAL(event_list)[k+3*L];
      tj=REAL(event_list)[r+3*L];
      
      REAL(dA)[i+j*p] =REAL(dA)[ i+j*p]+        (tj-ti)*exp(-gamma*(ti-tj));
      REAL(d2A)[i+j*p]=REAL(d2A)[i+j*p]+(tj-ti)*(tj-ti)*exp(-gamma*(ti-tj));
    }
  }
  
  // Prepare Output
  out=PROTECT(allocVector(VECSXP,2));
  SET_VECTOR_ELT(out,0,dA);
  SET_VECTOR_ELT(out,1,d2A);
  
  // Clean Up
  UNPROTECT(3);
  
  return(out);
}


SEXP compute_G(SEXP p_R,SEXP event_list,SEXP gammaR,SEXP TR,SEXP baseline,SEXP times)
{
  // Variables
  int i,j,k,r;
  int p,L,N;
  double gamma,tk,tkmin,tj;
  double avg_base;
  SEXP G;
  
  // Read Information
  p=asInteger(p_R);
  L=LENGTH(event_list)/4;
  gamma=asReal(gammaR);
  N=LENGTH(baseline)/p;
  
  // Allocate Memory
  G=PROTECT(allocVector(REALSXP,p*p));
  for(i=0;i<p*p;i++)
    REAL(G)[i]=0;
  
  // Compute G
  for(i=0;i<p;i++)
  {
    for(k=1;k<N;k++)
    {
      avg_base=REAL(baseline)[i+(k-1)*p];
      tk=REAL(times)[k];
      tkmin=REAL(times)[k-1];
      for(r=0;r<L;r++)
      {
        tj=REAL(event_list)[r+3*L];
        if(tj<tk)
        {
          j=(int)REAL(event_list)[r+2*L]-1;
      
          REAL(G)[i+j*p]=REAL(G)[i+j*p]+avg_base/gamma*(exp(-gamma*myclip(tkmin-tj))-exp(-gamma*(tk-tj)));
        }
      }
    }
  }
  
  // Clean Up
  UNPROTECT(1);
  
  return(G);
}

// Compute Derivatives of Matrix G
SEXP compute_deriv_G(SEXP p_R,SEXP event_list,SEXP gammaR,SEXP baseline,SEXP covariates)
{
  // Variables
  int i,j,k,l,r,m;
  int p,q,L,N;
  double gamma,tk,tkmin,tj;
  double avg_base;
  double M,diff;
  SEXP dgamma;
  SEXP d2gamma;
  SEXP dbeta;
  SEXP dgammabeta;
  SEXP d2beta;
  SEXP help1;
  SEXP help2;
  SEXP out;
  
  // Read Information
  p=asInteger(p_R);
  L=LENGTH(event_list)/4;
  gamma=asReal(gammaR);
  N=LENGTH(baseline)/p;
  q=LENGTH(VECTOR_ELT(VECTOR_ELT(covariates,1),0))/p;
  
  // Allocate Memory
  dgamma=PROTECT(allocVector(REALSXP,p*p));
  d2gamma=PROTECT(allocVector(REALSXP,p*p));
  dbeta=PROTECT(allocVector(VECSXP,q));
  dgammabeta=PROTECT(allocVector(VECSXP,q));
  d2beta=PROTECT(allocVector(VECSXP,q));
  
  for(i=0;i<q;i++)
  {
    help1=PROTECT(allocVector(REALSXP,p*p));
    SET_VECTOR_ELT(dbeta,i,help1);
    
    help1=PROTECT(allocVector(REALSXP,p*p));
    SET_VECTOR_ELT(dgammabeta,i,help1);
    
    help1=PROTECT(allocVector(VECSXP,q));
    for(j=0;j<q;j++)
    {
      help2=PROTECT(allocVector(REALSXP,p*p));
      SET_VECTOR_ELT(help1,j,help2);
    }
    SET_VECTOR_ELT(d2beta,i,help1);
  }
  
  // Initialize
  for(i=0;i<p*p;i++)
  {
    REAL( dgamma)[i]=0;
    REAL(d2gamma)[i]=0;
    for(k=0;k<q;k++)
    {
      REAL(VECTOR_ELT(dbeta,k))[i]=0;
      REAL(VECTOR_ELT(dgammabeta,k))[i]=0;
      for(l=0;l<q;l++)
        REAL(VECTOR_ELT(VECTOR_ELT(d2beta,k),l))[i]=0;
    }
  }
  
  // Compute derivatives of G
  for(i=0;i<p;i++)
  {
    for(k=1;k<N;k++)
    {
      avg_base=REAL(baseline)[i+(k-1)*p];
      tk=REAL(VECTOR_ELT(covariates,0))[k];
      tkmin=REAL(VECTOR_ELT(covariates,0))[k-1];
      for(r=0;r<L;r++)
      {
        tj=REAL(event_list)[r+3*L];
        if(tj<tk)
        {
          j=(int)REAL(event_list)[r+2*L]-1;
          
          M=myclip(tkmin-tj);
          diff=tk-tj;
          
          REAL(dgamma )[i+j*p]=REAL(dgamma )[i+j*p]+avg_base*( exp(-gamma*diff)*(                           1/(gamma*gamma)+     diff/gamma)-exp(-gamma*M)*(                        1/(gamma*gamma)+  M/gamma));
          REAL(d2gamma)[i+j*p]=REAL(d2gamma)[i+j*p]+avg_base*(-exp(-gamma*diff)*(2/(gamma*gamma*gamma)+2*diff/(gamma*gamma)+diff*diff/gamma)+exp(-gamma*M)*(2/(gamma*gamma*gamma)+2*M/(gamma*gamma)+M*M/gamma));
          
          for(l=0;l<q;l++)
          {
            REAL(VECTOR_ELT(dbeta     ,l))[i+j*p]=REAL(VECTOR_ELT(dbeta     ,l))[i+j*p]+avg_base*REAL(VECTOR_ELT(VECTOR_ELT(covariates,1),k))[i+l*p]/gamma*(exp(-gamma*M)                          -exp(-gamma*diff));
            REAL(VECTOR_ELT(dgammabeta,l))[i+j*p]=REAL(VECTOR_ELT(dgammabeta,l))[i+j*p]+avg_base*REAL(VECTOR_ELT(VECTOR_ELT(covariates,1),k))[i+l*p]*(     -exp(-gamma*M)*(1/(gamma*gamma)+M/gamma)+exp(-gamma*diff)*(1/(gamma*gamma)+diff/gamma));
            
            for(m=0;m<q;m++)
              REAL(VECTOR_ELT(VECTOR_ELT(d2beta,l),m))[i+j*p]=REAL(VECTOR_ELT(VECTOR_ELT(d2beta,l),m))[i+j*p]+avg_base*REAL(VECTOR_ELT(VECTOR_ELT(covariates,1),k))[i+l*p]*REAL(VECTOR_ELT(VECTOR_ELT(covariates,1),k))[i+m*p]/gamma*(exp(-gamma*M)-exp(-gamma*diff));
          }
        }
      }
    }
  }
  
  // Prepare Output
  out=PROTECT(allocVector(VECSXP,5));
  SET_VECTOR_ELT(out,0,dgamma);
  SET_VECTOR_ELT(out,1,d2gamma);
  SET_VECTOR_ELT(out,2,dbeta);
  SET_VECTOR_ELT(out,3,dgammabeta);
  SET_VECTOR_ELT(out,4,d2beta);
  
  // Clean Up
  UNPROTECT(6+3*q+q*q);
  
  return(out);
}

// Compute derivatives of Matrix V
SEXP dbeta_V(SEXP nu0,SEXP covariates)
{
  // Variables
  int p,q,T;
  int i,j,k,t;
  SEXP dbetaV;
  SEXP total2deriv;
  SEXP help;
  SEXP out;
  
  
  // Read Information
  T=LENGTH(VECTOR_ELT(covariates,0));
  p=LENGTH(nu0)/T;
  q=LENGTH(VECTOR_ELT(VECTOR_ELT(covariates,1),0))/p;
  
  // Allocate Memory
  dbetaV=PROTECT(allocVector(REALSXP,p*q));
  total2deriv=PROTECT(allocVector(VECSXP,q));
  for(i=0;i<q;i++)
  {
    help=PROTECT(allocVector(REALSXP,p*q));
    SET_VECTOR_ELT(total2deriv,i,help);
  }
  
  // Compute Derivatives
  for(k=0;k<q;k++)
  {
    for(i=0;i<p;i++)
    {
      REAL(dbetaV)[i+k*p]=0;
      for(j=0;j<q;j++)
        REAL(VECTOR_ELT(total2deriv,j))[i+k*p]=0;
      
      for(t=0;t<T-1;t++)
      {
        REAL(dbetaV)[i+k*p] =REAL(dbetaV )[i+k*p]+(2*REAL(VECTOR_ELT(VECTOR_ELT(covariates,1),t))[i+k*p]                                                    *REAL(nu0)[i+t*p]*REAL(nu0)[i+t*p])*(REAL(VECTOR_ELT(covariates,0))[t+1]-REAL(VECTOR_ELT(covariates,0))[t]);
        for(j=0;j<q;j++)
          REAL(VECTOR_ELT(total2deriv,j))[i+k*p]=REAL(VECTOR_ELT(total2deriv,j))[i+k*p]+(4*REAL(VECTOR_ELT(VECTOR_ELT(covariates,1),t))[i+k*p]*REAL(VECTOR_ELT(VECTOR_ELT(covariates,1),t))[i+j*p]*REAL(nu0)[i+t*p]*REAL(nu0)[i+t*p])*(REAL(VECTOR_ELT(covariates,0))[t+1]-REAL(VECTOR_ELT(covariates,0))[t]);
      }
    }
  }
  
  // Prepare Outout
  out=PROTECT(allocVector(VECSXP,2));
  SET_VECTOR_ELT(out,0,dbetaV );
  SET_VECTOR_ELT(out,1,total2deriv);
  
  // Clean Up
  UNPROTECT(3+q);
  
  return(out);
}

// Computes first derivative of the criterion function
SEXP compute_derivatives(SEXP V_diag,SEXP alpha,SEXP vector_v,SEXP C,SEXP G,SEXP Gamma,SEXP A,SEXP qR)
{
  // Variables
  SEXP grad;
  int p,q;
  int i,j,k;
  
  // Read Information
  p=LENGTH(alpha);
  q=asInteger(qR);
  
  // Allocate Memory
  grad=PROTECT(allocVector(REALSXP,q+1+p+p*p));
  
  // Compute Gradient
  // First entries are equal to zero
  for(i=0;i<q+1;i++)
    REAL(grad)[i]=0;
  
  // Derivative with respect to alpha
  for(i=0;i<p;i++)
  {
    REAL(grad)[1+q+i]=2*REAL(V_diag)[i]*REAL(alpha)[i]-2*REAL(vector_v)[i];
    for(j=0;j<p;j++)
      REAL(grad)[1+q+i]=REAL(grad)[1+q+i]+2*REAL(C)[i+j*p]*REAL(G)[i+j*p];
  }
  
  // Derivative with respect to C
  for(i=0;i<p;i++)
  {
    for(j=0;j<p;j++)
    {
      REAL(grad)[1+q+p+i+j*p]=2*REAL(alpha)[i]*REAL(G)[i+j*p]-2*REAL(A)[i+j*p];
      for(k=0;k<p;k++)
        REAL(grad)[1+q+p+i+j*p]=REAL(grad)[1+q+p+i+j*p]+REAL(C)[i+k*p]*REAL(Gamma)[k+j*p]+REAL(C)[i+k*p]*REAL(Gamma)[j+k*p];
    }
  }
  
  // Clean Up
  UNPROTECT(1);

  return(grad);
}

// Compute Second Derivative of criterion function with respect to alpha and C
SEXP compute_alphaC_deriv(SEXP G)
{
  // Variables
  int p;
  int i,j,k;
  SEXP out;
  
  // Read Information
  p=(int)sqrt((double)LENGTH(G));
  
  // Allocate Memory
  out=PROTECT(allocVector(REALSXP,p*p*p));
  
  // Compute Derivative
  for(k=0;k<p;k++)
  {
    for(i=0;i<p;i++)
    {
      for(j=0;j<p;j++)
      {
        if(i==k)
          REAL(out)[k+(i+j*p)*p]=REAL(G)[k+j*p];
        else
          REAL(out)[k+(i+j*p)*p]=0;
      }
    }
  }
  
  // Clean Up
  UNPROTECT(1);
  
  return(out);
}


// Compute Second Derivative of criterion function with respect to C
SEXP compute_d2C_deriv(SEXP Gamma)
{
  // Variables
  int p;
  int i,j,k,l;
  SEXP out;
  
  // Read Information
  p=(int)sqrt((double)LENGTH(Gamma));
  
  // Allocate Memory
  out=PROTECT(allocVector(REALSXP,p*p*p*p));
  
  // Compute Derivative
  for(i=0;i<p;i++)
  {
    for(j=0;j<p;j++)
    {
      for(k=0;k<p;k++)
      {
        for(l=0;l<p;l++)
          if(i==k)
            REAL(out)[i+j*p+(k+l*p)*p*p]=REAL(Gamma)[j+l*p]+REAL(Gamma)[l+j*p];
          else
            REAL(out)[i+j*p+(k+l*p)*p*p]=0;
      }
    }
  }
  
  // Clean Up
  UNPROTECT(1);
  
  return(out);
}