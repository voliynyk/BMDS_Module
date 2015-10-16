#include "sys/types.h"
#include "sys/timeb.h"

 void
gettimeofday_(long *tp, long *tzp)
{
	struct timeb tb;

	ftime(&tb);
	tp[0] = tb.time;
	tp[1] = tb.millitm;
	tzp[0] = tb.timezone;
	tzp[1] = tb.dstflag;
}

int times_(int *b)
{ 
  b[0] = 0;
  return 0;
}

