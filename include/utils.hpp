#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <fmt/format.h>
#include <cstdio>
#include <string>
#include <algorithm>

/**
 * Following functions and macros thanks to Dr. Lionel Barnett.
 * Modified by Joshua Brown
 */
#define ERRPT       fmt::print(stderr,"ERROR in '{}' [{}:{}]: ",__FUNCTION__,__FILE__,__LINE__)
#define EEXIT(...)  {ERRPT; fmt::print(stderr,__VA_ARGS__); fputc('\n',stderr); exit(EXIT_FAILURE);}
#define PEEXIT(...) {ERRPT; fmt::print(stderr,__VA_ARGS__); fputc('\n',stderr); perror(NULL); exit(EXIT_FAILURE);}

void magic_header(FILE* fp, const char* const hmessage);
bool check_magic_header(FILE* fp);
void progrep(const char* const msg, const size_t u, const size_t U);
void progrep(FILE* fp, const char* const msg, const size_t u, const size_t U);
void print_debug_info();

 /**
  * End Barnett functions
  */

#ifdef WIN32
inline FILE* popen(const char* command, const char* mode) {
	return _popen(command, mode);
}

inline int pclose(FILE* stream) {
	return _pclose(stream);
}
#endif

void make_dir(std::string path);

// Sets the OMP thread limit. Uses maximum threads if threadLimit = 0
int set_num_threads(int threadLimit, bool verbose = false);

#endif