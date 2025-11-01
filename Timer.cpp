#include "Timer.h"
#include <stdexcept>
#include <iomanip>
#include <chrono> // 念のためこちらもインクルード

// コンストラクタ：タイマーを初期化し、即時スタートします。
Timer::Timer()
	:
	running(true) // m_stop は stop() で設定されるため初期化不要
{
	m_start = std::chrono::high_resolution_clock::now(); // clock() から変更
}

// タイマーを（再）スタートします。
void Timer::start() {
	running = true;
	m_start = std::chrono::high_resolution_clock::now(); // clock() から変更
}

// タイマーをストップします。
void Timer::stop() {
	m_stop = std::chrono::high_resolution_clock::now(); // clock() から変更
	if (!running) {
		throw std::runtime_error(
			" In void Timer::stop() \n"
			" timer is not started.\n");
	}
	running = false;
}

// 経過時間を出力ストリームに書き出します。
std::ostream& operator<<(std::ostream& ostm, const Timer& timer) {
	if (timer.running) {
		throw std::runtime_error(
			" In std::ostream& operator<<(std::ostream& ostm, const Timer& timer) \n"
			" timer is running.\n");
	}

    // 経過時間をミリ秒(double)で計算します。
	std::chrono::duration<double, std::milli> time_ms = timer.m_stop - timer.m_start;
	
	ostm << "time " << std::setw(10) << time_ms.count() << "[ms]";
	return ostm;
}