// This is a simple implementation for logging, which was modeled after 
//
// http://stackoverflow.com/questions/6168107/how-to-implement-a-good-debug-logging-feature-in-a-project
//
//
#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <iostream>
#include <sstream>
enum loglevel_e
    {logERROR=0, logWARNING=1, logINFO=2, logDEBUG=3, logDEBUG1=4, logDEBUG2=5, logDEBUG3=6, logDEBUG4=7};


class logIt
{
	public:
		logIt(loglevel_e _loglevel = logERROR) {
			_buffer << _loglevel <<" :"
				<< std::string(
					_loglevel > logDEBUG 
					? (_loglevel - logDEBUG) * 4 
					: 1
					, ' ');
		}

		template <typename T>
			logIt & operator<<(T const & value)
			        {
					        _buffer << value;
						        return *this;
							    }
		~logIt()
		{
		        _buffer << std::endl;
			std::cerr << _buffer.str();
		}
	private:
		std::ostringstream _buffer;
};

// global loglevel as reference
extern loglevel_e loglevel;

// actual logging function
#define Log(level) \
	if(level > loglevel); \
	else logIt(level)

#endif

