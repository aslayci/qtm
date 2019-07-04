#include "utils.h"

using namespace std;

long double getElapsedWallTime(std::chrono::steady_clock::time_point begin, std::chrono::steady_clock::time_point end) {
    
    long double timeDiff;
    timeDiff = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    //    timeDiff /= 1000.0; // milliseconds
    timeDiff /= 1000000.0; // seconds
    return timeDiff;
}

long double getElapsedCpuTime(std::clock_t c_start, std::clock_t c_end) {
    
//    long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    long double time_elapsed_sec  = (1000.0 * (c_end-c_start) / CLOCKS_PER_SEC) / 1000.0;

    return time_elapsed_sec;
}



/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t getPeakRSS() {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakWorkingSetSize;
    
#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
        return (size_t)0L;        /* Can't open? */
    if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
    {
        close( fd );
        return (size_t)0L;        /* Can't read? */
    }
    close( fd );
    return (size_t)(psinfo.pr_rssize * 1024L);
    
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif
    
#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;            /* Unsupported. */
#endif
}

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t getCurrentRSS() {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.WorkingSetSize;
    
#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
                   (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
        return (size_t)0L;        /* Can't access? */
    return (size_t)info.resident_size;
    
#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE* fp = NULL;
    if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
        return (size_t)0L;        /* Can't open? */
    if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    {
        fclose( fp );
        return (size_t)0L;        /* Can't read? */
    }
    fclose( fp );
    return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);
    
#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;            /* Unsupported. */
#endif
}



//// for now only for linux -- make it OS generic later
float getCurrentMemoryUsage() { //megabyte verii
    string pid = intToStr(int(getpid()));
//    string command1 = string("mkdir -p forMem");
//    system(command1.c_str());
    string outfile = std::string("forMem_") + pid + std::string(".txt");
    string command = "pmap " + pid + " | grep -i Total | awk '{print $2}' > " + outfile;
    system(command.c_str());
    string mem_str;
    ifstream ifs(outfile.c_str());
    std::getline(ifs, mem_str);
    ifs.close();
    string command1 = string("rm -f ") + outfile;
    system(command1.c_str());

    mem_str = mem_str.substr(0, mem_str.size()-1);
    float mem = (float)strToInt(mem_str);

    return mem/1024; // in MB
}


void split(string str, string delim, vector<string> &result){
//    result.clear();
    while(str.size()){
        int index = str.find(delim);
        if(index!=string::npos){
            string s = str.substr(0,index);
            result.push_back(trim(s));
            str = str.substr(index+delim.size());
            if(str.size()==0)result.push_back(trim(str));
        }else{
            result.push_back(trim(str));
            str = "";
        }
    }
//    return result;
}

void splitDouble(string str, string delim, vector<double> &result){
    //    result.clear();
    while(str.size()){
        int index = str.find(delim);
        if(index!=string::npos){
            string s = str.substr(0,index);
            result.push_back(strToDouble(s));
            str = str.substr(index+delim.size());
            if(str.size()==0)result.push_back(strToDouble(str));
        }else{
            result.push_back(strToDouble(str));
            str = "";
        }
    }
    //    return result;
}

void splitint(string str, string delim, vector<int> &result){

    while(str.size()){
        int index = str.find(delim);
        if(index!=string::npos){
            string s = str.substr(0,index);
            result.push_back(strtoul(s.c_str(),NULL,0));
            str = str.substr(index+delim.size());
            if(str.size()==0) result.push_back(strtoul(str.c_str(),NULL,0));
        }else{
            result.push_back(strtoul(str.c_str(),NULL,0));
            str = "";
        }
    }
}


bool read_next_token(ifstream &input_file, string &token) {
    while (input_file) {
        input_file >> token;
        if (token[0] != '#') return true;
        getline(input_file, token);  // ignore rest of line
    }
    return false;
}

bool read_next_double(ifstream &input_file, double &d) {
    string token;
    if (!read_next_token(input_file, token)) return false;
    d = stod(token);
    return true;
}

void stringTokenizer(string& str, float *tokens, int size, string& delimiters) {
	//Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	//Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);
	
	int i = 0;
	while (string::npos != pos || string::npos != lastPos)
	{
		if(i >= size) {
			cout << "something is wrong" << endl;
			exit(1);
		}
		
		tokens[i++] = strToFloat(str.substr(lastPos, pos - lastPos));
		lastPos = str.find_first_not_of(delimiters, pos);		
		pos = str.find_first_of(delimiters, lastPos);
	}
}
void stringTokenizer(string& str, double *tokens, int size, string& delimiters) {
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);
	
	int i = 0;
	while (string::npos != pos || string::npos != lastPos)
	{
		if(i >= size) {
			cout << "something is wrong" << endl;
			exit(1);
		}
			
		tokens[i++] = strToDouble(str.substr(lastPos, pos - lastPos));
		lastPos = str.find_first_not_of(delimiters, pos);		
		pos = str.find_first_of(delimiters, lastPos);
	}
}


//void rtrim(char *str) {
//    size_t n;
//    n = strlen(str);
//    while (n > 0 && isspace((int char)str[n - 1])) {
//        n--;
//    }
//    str[n] = '\0';
//}
//
//void ltrim(char *str) {
//    size_t n;
//    n = 0;
//    while (str[n] != '\0' && isspace((int char)str[n])) {
//        n++;
//    }
//    memmove(str, str + n, strlen(str) - n + 1);
//}
//
//void trim(char *str) {
//    rtrim(str);
//    ltrim(str);
//}
