#include <R.h>

void F77_NAME(fexitc)(char *msg, int *nchar)
{
    int nc = *nchar;
    char buf[256];
    if(nc > 255) {
        warning("invalid character length in fexitc");
        nc = 255;
    }
    strncpy(buf, msg, nc);
    buf[nc] = '\0';
    error(buf);
}
