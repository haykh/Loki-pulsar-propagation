#pragma PC_DENS

double gFunc (double R);

#ifdef INTBACK
  void rpFromR (double Rmax);
  namespace pcdens {
    extern double *rps;
    extern double *Rs;
    extern int N;
  }
#endif