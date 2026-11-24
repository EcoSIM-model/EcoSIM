#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>

void timec_(double* t1)
{

  struct timeval t;
  gettimeofday(&t,0);
  *t1=(double)t.tv_usec/1e6+(double)t.tv_sec;
}
