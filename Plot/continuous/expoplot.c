#include <stdio.h>
#include "bmds_plot.h"
int main(int argc, char *argv[])
{
  int i;

#ifdef DO_LOG
  FILE *fplog;
  fplog = fopen("expoplot-debug.log", "w+");
  for (i=0; i < argc; i++) fprintf(fplog, "argv[%d] = {%s}\n", i, argv[i]);
  fclose(fplog);
#endif

  plot_continuous(evExponential, argc, argv);
  return 0;
}
