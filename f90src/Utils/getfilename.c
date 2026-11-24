#include <stdio.h>
#include <time.h>

// A function that terminates when enter key is pressed

void getfilename_(char *s1, char *s2, char *s3);

void getfilename_(char *s1, char *s2, char *s3){

  while( (*s3++ = *s1++) != ' ')
   ;

  if(*s1 !='\0')s3--;
  while( (*s3++ = *s2++) != ' ')
   ;
  if(*s2 !='\0')s3--;
  *s3='\0';

}
