#include <stdarg.h>
#include <gkyl_app_priv.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#endif

// private function to handle variable argument list for printing
static void
v_app_cout(FILE *fp, const char *fmt, va_list argp)
{
  int rank = 0;
#ifdef GKYL_HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif  
  if ((rank == 0) && fp) {
    vfprintf(fp, fmt, argp);
    fflush(fp);
  }
}

void
gkyl_app_cout(FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_app_cout(fp, fmt, argp);
  va_end(argp);  
}
      
