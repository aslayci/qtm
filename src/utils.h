#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <iterator>
#include <sys/types.h>
#include <ctime>
#include <sstream>
#include <climits>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <chrono>

#if defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>
#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

typedef long long int int64;
//typedef long long int uint64;
typedef long long unsigned uint64;

size_t getPeakRSS();
size_t getCurrentRSS();


#ifdef _WIN32
#define OS_SEP "\\"
#else
#define OS_SEP "/"
#endif

namespace qtm {
    
    struct VectorHasher {
        std::size_t operator()(std::vector<int> const& vec) const {
            std::size_t seed = vec.size();
            for(auto& i : vec) {
                seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
    
    typedef std::map<std::vector<int>, double> CPT;
    
    struct Query {
        std::vector<int> queryVars;
        std::vector<int> evidenceStates;
    };
    
    typedef std::vector<int> Vi;
    typedef std::vector<Vi> Vvi;

    inline std::ostream& operator<<(std::ostream& os, const Vi& vi) {
        os << "(";
        std::copy(vi.begin(), vi.end(), std::ostream_iterator<int>(os, ", "));
        os << ")";
        return os;
    }
    inline std::ostream& operator<<(std::ostream& os, const Vvi& vvi) {
        os << "(\n";
        for(Vvi::const_iterator it = vvi.begin();
            it != vvi.end();
            it++) {
            os << "  " << *it << "\n";
        }
        os << ")";
        return os;
    }
    
    struct Digits {
        Vi::const_iterator begin;
        Vi::const_iterator end;
        Vi::const_iterator me;
    };
    
    typedef std::vector<Digits> Vd;
    
};

using namespace std;

//struct timespec {
//    time_t tv_sec; /* seconds */
//    long tv_nsec; /* nanoseconds */
//};

//double getRunningTime(std::chrono::steady_clock::time_point begin, std::chrono::steady_clock::time_point end);
long double getElapsedWallTime(std::chrono::steady_clock::time_point begin, std::chrono::steady_clock::time_point end);
long double getElapsedCpuTime(std::clock_t c_start, std::clock_t c_end);

float getCurrentMemoryUsage();

void split(string str, string delim, vector<string> &result);
void splitDouble(string str, string delim, vector<double> &result);
void splitint(string str, string delim, vector<int> &result);

void stringTokenizer(string& str, float *tokens, int size, string& delimiters);
void stringTokenizer(string& str, double *tokens, int size, string& delimiters);
bool read_next_token(ifstream &input_file, string &token);
bool read_next_double(ifstream &input_file, double &d);

// trim from start
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}


inline int strToInt(string s) {
  int i;
  istringstream myStream(s);

  if(myStream >> i) 
    return i;
  
  else {
    std::cout << "String " << s << " is not a number." << endl;
    return atoi(s.c_str());
  }
  
  return i;
}

inline int strToInt(const char* s) {
  int i;
  istringstream myStream(s);

  if (myStream >> i)
    return i;
  
  else {
    cout << "String " << s << " is not a number." << endl;
    return atoi(s);
    }
  
  return i;
}

inline float strToFloat(const char* s) {
  return atof(s);  
}

inline float strToFloat(string s) {
  return atof(s.c_str());
}

inline string floatToStr(float f) {
  stringstream ss;
  ss << f;
  return ss.str();
}

inline int64_t strToInt64(string s) {
  int64_t i;
  istringstream myStream(s);

  if (myStream >> i)
    return i;
  
  else {
    cout << "String " << s << " is not a number." << endl;
    exit(1);    
  }
  
  return i;
}

inline string intToStr(int i) {
  stringstream ss;
  ss << i;
  return ss.str();  
}

inline double strToDouble(string s) {	
// 	return std::stod(s);
	double a = 0.0; 
	stringstream ss;
	ss << s;
	ss >> a;
	return a;
}

inline std::string lowercase(std::string str) {
    std::transform(str.begin(),
                   str.end(),
                   str.begin(),
                   [](char const &c){
                       return ::tolower(c);
                   });
    return str;
}



#endif
