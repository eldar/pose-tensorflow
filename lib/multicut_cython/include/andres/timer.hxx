#pragma once
#ifndef ANDRES_TIMER_HXX
#define ANDRES_TIMER_HXX

#include <chrono>

namespace andres {

class Timer {
public:
    Timer();
    
    double elapsedSeconds() const;

    void reset();
    void start();
    void stop();

private:
    double seconds_;

    decltype(std::chrono::high_resolution_clock::now()) time_;
};

inline
Timer::Timer()
:   seconds_()
{}

inline void
Timer::reset() {
    seconds_ = .0;
}

inline void
Timer::start() {
    time_ = std::chrono::high_resolution_clock::now();
}

inline void
Timer::stop() {
    seconds_ += std::chrono::duration_cast<std::chrono::duration<double>>(
        std::chrono::high_resolution_clock::now() - time_
    ).count();
}

inline double
Timer::elapsedSeconds() const {
    return seconds_;
}

} // namespace andres

#endif
