// http://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows

//  Windows
#ifdef _WIN32
#include <Windows.h>

inline double get_wall_time()
{
    LARGE_INTEGER time, freq;
    if (!QueryPerformanceFrequency(&freq))
    {
        //  Handle error
        return 0;
    }
    if (!QueryPerformanceCounter(&time))
    {
        //  Handle error
        return 0;
    }
    return (double) time.QuadPart / freq.QuadPart;
}

inline double get_cpu_time()
{
    FILETIME a, b, c, d;
    if (GetProcessTimes(GetCurrentProcess(), &a, &b, &c, &d) != 0)
    {
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return (double) (d.dwLowDateTime |
                         ((unsigned long long) d.dwHighDateTime << 32)) *
            0.0000001;
    }
    else
    {
        //  Handle error
        return 0;
    }
}

//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>

inline double get_wall_time()
{
    struct timeval time;
    if (gettimeofday(&time, NULL))
    {
        //  Handle error
        return 0;
    }
    return (double) time.tv_sec + (double) time.tv_usec * .000001;
}

inline double get_cpu_time()
{
    return (double) clock() / CLOCKS_PER_SEC;
}

#endif
