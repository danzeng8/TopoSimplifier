/**
 * \file   timer.cpp
 * \brief  thread-safe timer implementation
 *
 * \author Max Resch
 * \date   2012-12-21
 */

#include "timer.h"
#include <boost/timer/timer.hpp>

using namespace std;
using namespace boost::timer;

const Timer Timer::total(true);

void Timer::start()
{
	unique_lock<mutex>(atomic);
	timer.start();
}

const CPUTime Timer::stop()
{
	unique_lock<mutex>(atomic);
	timer.stop();
	return CPUTime(timer.elapsed());
}

const CPUTime Timer::elapsed() const
{
	unique_lock<mutex>(atomic);
	return CPUTime(timer.elapsed());
}

void Timer::pause()
{
	unique_lock<mutex>(atomic);
	timer.stop();
}

void Timer::resume()
{
	unique_lock<mutex>(atomic);
	timer.resume();
}

Timer::Timer(bool autostart)
{
	if (autostart)
		timer.start();
}

Timer::~Timer()
{
}
