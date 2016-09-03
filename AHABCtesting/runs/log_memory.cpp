#include "log_memory.hpp"

#include <unistd.h> // sysconf(_SC_CLK_TCK), posix
#include <ios>

//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

void process_mem_usage(double& vm_usage,
						double& resident_set,
						double& process_time)
{
	using std::ios_base;
	using std::ifstream;
	using std::string;

	vm_usage     = 0.0;
	resident_set = 0.0;

	// 'file' stat seems to give the most reliable results
	//
	ifstream stat_stream("/proc/self/stat",ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	//
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want
	//
	unsigned long vsize, utime, stime;
	long rss, cutime, cstime;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
			   >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
			   >> utime >> stime >> cutime >> cstime >> priority >> nice
			   >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

	stat_stream.close();

	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage     = vsize / 1024.0;
	resident_set = rss * page_size_kb;

	double clock_ticks = static_cast<double>(sysconf(_SC_CLK_TCK)); 
	process_time = (utime + stime + cutime + cstime)/clock_ticks;
}

void logUsage(std::ostream& os)
{
	double vm, rss, pt;
	process_mem_usage(vm, rss, pt);
	os << "VM:\t" << vm << "\tRSS:\t" << rss << "\ttime:\t" << pt << std::endl;
}

void logUsage(std::ostream& os, const std::string& logline)
{
	double vm, rss, pt;
	process_mem_usage(vm, rss, pt);
	os << logline << std::endl;
	os << "\t" << "VM:\t" << vm << "\tRSS:\t" << rss << "\ttime:\t" << pt << std::endl;
}

void logUsage(const std::string& filename, 
				const std::string& logline,
				bool append)
{
	std::ofstream os;
	
	if (append) os.open(filename.c_str(), std::ios::app);         // append
	else os.open(filename.c_str()); // don't append
		
	if (os.is_open()) {
		logUsage(os, logline);
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
		<< filename << std::endl << std::endl;
	}
		
   double vm, rss, pt;
   process_mem_usage(vm, rss, pt);
   os << logline << std::endl;
   os << "\t" << "VM:\t" << vm << "\tRSS:\t" << rss << "\ttime:\t" << pt << std::endl;
}
