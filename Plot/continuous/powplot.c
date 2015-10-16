#include <stdio.h>
#include "bmds_plot.h"
int main(int argc, char *argv[])
{
  int i;

#ifdef DO_LOG
  FILE *fplog;
  fplog = fopen("powplot-debug.log", "w+");
  for (i=0; i < argc; i++) fprintf(fplog, "argv[%d] = {%s}\n", i, argv[i]);
  fclose(fplog);
#endif

  plot_continuous(evPower, argc, argv);
  return 0;
}
