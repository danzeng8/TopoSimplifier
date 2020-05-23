/**
 * \file   timer.cpp
 * \brief  thread-safe timer
 *
 * \author Max Resch
 * \date   2012-12-21
 */

#ifndef TIMER_H_
#define TIMER_H_

#include "cputime.h"

#include <mutex>
#include <boost/timer/timer.hpp>
#include <boost/chrono.hpp>

class Timer
{
public:
	Timer(bool autostart = false);

	virtual ~Timer();

	void start();

	const CPUTime elapsed() const;

	const CPUTime stop();

	void pause();

	void resume();

	const static Timer total;

private:
	std::mutex atomic;

	boost::timer::cpu_timer timer;
};

#endif /* TIMER_H_ */
