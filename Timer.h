#pragma once
#include <iostream>
#include <time.h>
class Timer
{
private:
	long long m_start;
	long long m_stop;
	bool running;
public:
	void start();
	void stop();
	Timer();
	friend std::ostream& operator<<(std::ostream& ostm, const Timer& timer);

};


