/**
 * \file   time.h
 * \brief  data structure for storing wall, user and system time
 *
 * \author Max Resch
 * \date   2013-04-21
 */

#ifndef TIME_H_
#define TIME_H_

#include <boost/timer/timer.hpp>
#include <chrono>
#include <ostream>

typedef std::chrono::duration<double, std::ratio<1>> seconds;

struct CPUTime
{
	typedef std::chrono::nanoseconds _time;

	CPUTime();
	CPUTime(const boost::timer::cpu_times&);

	_time wall;
	_time user;
	_time sys;

	CPUTime& operator+=(const CPUTime&);
	CPUTime operator+(const CPUTime&) const;
	CPUTime operator-(const CPUTime&) const;

	double getSeconds() const;
};

std::ostream& operator<<(std::ostream& os, const CPUTime t);

#endif /* TIME_H_ */
