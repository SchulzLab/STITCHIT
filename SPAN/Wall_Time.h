#ifndef WALL_TIME
#define WALL_TIME

#ifdef __APPLE__
#include <mach/mach.h>
#include <mach/mach_time.h>

double getTime(){
	static double timeConvert = 0.0;
	if ( timeConvert == 0.0 ){
		mach_timebase_info_data_t timeBase;
		(void)mach_timebase_info( &timeBase );
		timeConvert = (double)timeBase.numer / (double)timeBase.denom / 1000000000.0;
	}
	return (double)mach_absolute_time( ) * timeConvert;
}
#elif __linux
#include <sys/time.h>

double getTime(){
	struct timeval tim;
	gettimeofday(&tim, NULL);
	return tim.tv_sec + (tim.tv_usec / 1000000.0);
}
#endif

#endif //WALL_TIME