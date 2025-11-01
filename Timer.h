#pragma once
#include <iostream>
#include <chrono> // <time.h> の代わりに追加

class Timer
{
private:
    // 時間点を保持する型を変更
	std::chrono::high_resolution_clock::time_point m_start;
	std::chrono::high_resolution_clock::time_point m_stop;
	bool running;

public:
	void start();
	void stop();
	Timer(); // コンストラクタ
	friend std::ostream& operator<<(std::ostream& ostm, const Timer& timer);
};