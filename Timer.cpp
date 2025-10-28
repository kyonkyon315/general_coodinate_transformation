#include "Timer.h"
#include <stdexcept>
#include <iomanip>

Timer::Timer()
	:
	running(true), m_stop(0)
{
	m_start = clock();
}

void Timer::start() {
	running = true;
	m_start = clock();
}

void Timer::stop() {
	m_stop = clock();
	if (!running) {
		throw std::runtime_error(
			" In void Timer::stop() \n"
			" timer is not started.\n");
	}
	running = false;
}

std::ostream& operator<<(std::ostream& ostm, const Timer& timer) {
	if (timer.running) {
		throw std::runtime_error(
			" In std::ostream& operator<<(std::ostream& ostm, const Timer& timer) \n"
			" timer is running.\n");
	}
	const double time = (static_cast<double>(timer.m_stop - timer.m_start)) / CLOCKS_PER_SEC * 1000.0;
	ostm << "time " << std::setw(10) << time << "[ms]";
	return ostm;
}






