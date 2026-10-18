// Minimal test: returns 0 on success, non-zero on failure (that's all CTest needs).

#include <iostream>

#include "nacawing/version.hpp"

int main()
{
    if (nacawing::version().empty()) {
        std::cerr << "version() returned an empty string\n";
        return 1;
    }
    std::cout << "version test passed\n";
    return 0;
}
