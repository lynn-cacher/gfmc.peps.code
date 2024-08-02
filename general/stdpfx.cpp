/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#include "stdpfx.h"

#ifdef _MSC_VER
Str present_time_international()
// date and time
{
	struct tm newtime;
	__time64_t seconds;
	char timebuf[26];
	errno_t err;

	_time64(&seconds);
	err = _localtime64_s(&newtime, &seconds);
	if (err) { printf("Invalid argument to _localtime64_s."); exit(1); }
	err = asctime_s(timebuf, 26, &newtime);
	if (err) { printf("Invalid argument to asctime_s."     ); exit(1); }
	timebuf[24] = '\0';
	return Str(timebuf);
}
#else
Str present_time_international()
// date and time
{
	struct tm *newTime;
	time_t szClock;
	time(&szClock);
	newTime = localtime(&szClock);
	char *t = asctime(newTime);
	t[24] = '\0';
	return Str(t);
}
#endif
