#include "../include/defs.hpp"
#include "../include/utils.hpp"
#include "../include/version.hpp"
#include <fmt/format.h>
#include <cstdint>
#include <omp.h>

bool isLittleEndian() {
	uint32_t num = 1;
	return (*(uint8_t*)& num == 1);
}
void magic_header(FILE * fp, const char* const hmessage)
{
	fmt::print(fp, "{}\n\n", hmessage);

	std::string order = "BIG ENDIAN";
	if (isLittleEndian()) {
		order = "LITTLE ENDIAN";
	}
	//Assumes float word order is same as integer
	fmt::print(fp, "byte  order : {}\n", order);
	fmt::print(fp, "float order : {}\n\n", order);

	fmt::print(fp, "IEEE 754 compliance : {}\n\n", std::numeric_limits<float>::is_iec559 ? "YES" : "NO");

	fmt::print(fp, "sizeof(char)        ={:3}\n", sizeof(char));
	fmt::print(fp, "sizeof(size_t)      ={:3}\n", sizeof(size_t));
	fmt::print(fp, "sizeof(int)         ={:3}\n", sizeof(int));
	fmt::print(fp, "sizeof(long)        ={:3}\n", sizeof(long));
	fmt::print(fp, "sizeof(double)      ={:3}\n\n", sizeof(double));

	fmt::print(fp, "BEGIN MAGIC\n");
	char xc[3] = { 'b','a','h' };
	if (fwrite(&xc, sizeof(char), 3, fp) < 3) PEEXIT("write failed");
	size_t xs = 12345678987654321L;
	if (fwrite(&xs, sizeof(size_t), 1, fp) < 1) PEEXIT("write failed");
	int x2 = 12345678, y2 = -x2;
	if (fwrite(&x2, sizeof(int), 1, fp) < 1) PEEXIT("write failed");
	if (fwrite(&y2, sizeof(int), 1, fp) < 1) PEEXIT("write failed");
	int64_t x3 = 1234567898765432L, y3 = -x3;
	if (fwrite(&x3, sizeof(int64_t), 1, fp) < 1) PEEXIT("write failed");
	if (fwrite(&y3, sizeof(int64_t), 1, fp) < 1) PEEXIT("write failed");
	double x4 = 12345678.987654321, y4 = -x4;
	if (fwrite(&x4, sizeof(double), 1, fp) < 1) PEEXIT("write failed");
	if (fwrite(&y4, sizeof(double), 1, fp) < 1) PEEXIT("write failed");
	fmt::print(fp, "\nEND MAGIC\n");
}

bool check_magic_header(FILE * fp)
{
	bool begin_magic = false;
	constexpr size_t len = 256;
	char line[len];
	while (fgets(line, len, fp) != NULL) {
		if (strcmp(line, "BEGIN MAGIC\n") == 0)
		{
			begin_magic = true;
			fmt::print("found magic\n");
			break;
		}
	}
	if (!begin_magic) PEEXIT("Magic did not begin");

	char xc[3] = { 'z', 'z', 'z' };
	size_t xs = 0;
	int x2 = 0;
	int64_t x3 = 0;
	double x4 = 0;
	if (fread(xc, sizeof(char), 3, fp) != 3) PEEXIT("WARNING - bad magic: char");
	if (xc[0] != 'b' || xc[1] != 'a' || xc[2] != 'h') PEEXIT("WARNING - bad magic: char");
	if (fread(&xs, sizeof(size_t), 1, fp) != 1 || xs != 12345678987654321L) PEEXIT("WARNING - bad magic: size_t");
	if (fread(&x2, sizeof(int), 1, fp) != 1 || x2 != 12345678) PEEXIT("WARNING - bad magic: int");
	if (fread(&x2, sizeof(int), 1, fp) != 1 || x2 != -12345678) PEEXIT("WARNING - bad magic: -int");
	if (fread(&x3, sizeof(int64_t), 1, fp) != 1 || x3 != 1234567898765432L) PEEXIT("WARNING - bad magic: long");
	if (fread(&x3, sizeof(int64_t), 1, fp) != 1 || x3 != -1234567898765432L) PEEXIT("WARNING - bad magic: -long");
	if (fread(&x4, sizeof(double), 1, fp) != 1 || x4 != 12345678.987654321) PEEXIT("WARNING - bad magic: double");
	if (fread(&x4, sizeof(double), 1, fp) != 1 || x4 != -12345678.987654321) PEEXIT("WARNING - bad magic: -double");

	if (fgets(line, len, fp) == NULL) PEEXIT("magic does not end"); //consume extra line
	if (fgets(line, len, fp) == NULL) PEEXIT("magic does not end");
	if (strcmp(line, "END MAGIC\n") != 0) PEEXIT("magic ends badly: {}", line);
	return true;
}

void progrep(const char* const msg, const size_t u, const size_t U)
{
	if ((u + 1) % (U / 10) == 0) fmt::print("{} : {}% complete\n", msg, 10 * (u + 1) / std::max(1ull, (U / 10))); // report progress of iterative sim from 10% - 100%
}

void progrep(FILE * fp, const char* const msg, const size_t u, const size_t U)
{
	if ((u + 1) % (U / 10) == 0) fmt::print(fp, "{} : {}% complete\n", msg, 10 * (u + 1) / std::max(1ull, (U / 10))); // report progress of iterative sim from 10% - 100%
}

void print_debug_info()
{
	fmt::print("Num States: {:d}\n", numStates);

	fmt::print("Git Head: {}\n", GIT_HEAD);
	fmt::print("Git Log: {}\n", GIT_LOG);
	fmt::print("Git Diff: {}\n", GIT_DIFF);
}

void make_dir(std::string path) {
#ifdef WIN32
	// Convert delimiters and just use mkdir <dir> (handles subdirs by default)
	std::replace(path.begin(), path.end(), '/', '\\');
	std::string cmd = fmt::format("mkdir {}", path);
#else
	// Otherwise need to use -p to make intermediate directories
	std::string cmd = fmt::format("mkdir -p {}", path);
#endif
	if (system(cmd.c_str()) == -1) PEEXIT("Output dir \"%s\" creation failed", path.c_str());
}

int set_num_threads(int threadLimit, bool verbose)
{
	int num_threads = omp_get_max_threads();
	if (threadLimit > 0) {
		omp_set_num_threads(threadLimit);
		num_threads = threadLimit;//omp_get_num_threads();
	}
	else {
		omp_set_num_threads(num_threads);
	}
	if (verbose) {
		fmt::print("Using {} threads (Max available: {})\n", num_threads, omp_get_max_threads());
	}
	return num_threads;
}
