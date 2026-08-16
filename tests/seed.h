#ifndef PSI_TESTS_SEED_H_
#define PSI_TESTS_SEED_H_

#include <cstdlib>
#include <iostream>
#include <string>

/*
 * Generated inputs used to be seeded from std::random_device, so a run that
 * failed could never be replayed. Fixed by default; override with PSI_SEED=<n>
 * or data_generator's -seed. Printed once so the value lands in the run's own
 * output.
 */
inline unsigned int generator_seed()
{
	static unsigned int base = [] {
		char const *s = std::getenv("PSI_SEED");
		unsigned int v = s ? static_cast<unsigned int>(
					     std::strtoul(s, nullptr, 10))
				   : 20260816u;
		std::cout << "generator seed " << v << std::endl;
		return v;
	}();
	/* a distinct stream per call site, still reproducible */
	static unsigned int nth = 0;
	return base + nth++;
}

#endif /* PSI_TESTS_SEED_H_ */
