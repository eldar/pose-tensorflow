#pragma once

#include <sstream>
#include <chrono>

namespace levinkov
{

class Timer
{
public:
	double get_elapsed_seconds() const
	{
		return m_seconds;
	}

	unsigned long long hours() const
	{
		return static_cast<unsigned long long>(m_seconds) / 3600ULL;
	}

	unsigned long long minutes() const
	{
		return (static_cast<unsigned long long>(m_seconds) - hours()*3600ULL) / 60ULL;
	}
	
	void reset()
	{
		m_seconds = .0;
	}

	double seconds() const
	{
		return m_seconds - 3600.0*hours() - 60.0*minutes();
	}

	void start()
	{
		m_timeObject = std::chrono::high_resolution_clock::now();
	}
	
	void stop()
	{
		m_seconds += std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - m_timeObject).count();
	}
	
	std::string to_string() const
	{
		std::ostringstream s(std::ostringstream::out);
		
		s << hours() << "h " << minutes() << "m " << seconds() << "s";
		
		return s.str();
	}
	
private:
	double m_seconds { .0 };

	decltype(std::chrono::high_resolution_clock::now()) m_timeObject;
};

}