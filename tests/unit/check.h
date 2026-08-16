#ifndef PSI_TESTS_UNIT_CHECK_H_
#define PSI_TESTS_UNIT_CHECK_H_

/*
 * A test is a program: it runs checks and returns psi_test_result(). No
 * framework, because the artifact has to configure and build offline.
 */

#include <cstdio>
#include <string>

inline int psi_test_failures = 0;
inline int psi_test_checks = 0;
inline char const *psi_test_case = "";

/* Reports the count as well as the verdict: a suite that silently stopped
 * checking anything would otherwise still print ok. */
inline int psi_test_result()
{
	if (psi_test_checks == 0) {
		std::printf("FAILED: no checks ran\n");
		return 1;
	}
	if (psi_test_failures == 0) {
		std::printf("ok (%d checks)\n", psi_test_checks);
		return 0;
	}
	std::printf("FAILED: %d of %d checks\n", psi_test_failures,
		    psi_test_checks);
	return 1;
}

inline void psi_test_fail(char const *file, int line, char const *expr,
			  std::string const &detail)
{
	std::printf("%s:%d: [%s] %s%s%s\n", file, line, psi_test_case, expr,
		    detail.empty() ? "" : " -- ", detail.c_str());
	psi_test_failures++;
}

/* Names the group a failure came from, so the output says where to look. */
#define CASE(name) psi_test_case = (name)

#define CHECK(expr)                                                   \
	do {                                                          \
		psi_test_checks++;                                    \
		if (!(expr))                                          \
			psi_test_fail(__FILE__, __LINE__, #expr, ""); \
	} while (0)

#define CHECK_EQ(a, b)                                                  \
	do {                                                            \
		psi_test_checks++;                                      \
		auto const psi_a = (a);                                 \
		auto const psi_b = (b);                                 \
		if (!(psi_a == psi_b))                                  \
			psi_test_fail(__FILE__, __LINE__, #a " == " #b, \
				      std::to_string(psi_a) + " vs " +  \
					      std::to_string(psi_b));   \
	} while (0)

#endif /* PSI_TESTS_UNIT_CHECK_H_ */
