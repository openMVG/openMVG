#include "stdio.h"
#include "dlfcn.h"
int main() { if (!dlopen("libpba.so", RTLD_LAZY)) printf("error"); }
