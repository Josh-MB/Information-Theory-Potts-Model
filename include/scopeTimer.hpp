#include <fmt/format.h>
#include <string>
#include <chrono>

/**
 Helper class for timing how long a scope takes to run
*/
class ScopeTimer {
private:
	std::string name;
	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
public:
	ScopeTimer(std::string name) : name(name) {}
	~ScopeTimer() {
		auto end = std::chrono::high_resolution_clock::now();
		int64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		fmt::print("{} took {} ms\n", name, ms);
	}
	ScopeTimer(ScopeTimer const& other) = delete;
	ScopeTimer(ScopeTimer && other) = delete;
	ScopeTimer& operator=(ScopeTimer const& rhs) = delete;
	ScopeTimer& operator=(ScopeTimer && rhs) = delete;
};