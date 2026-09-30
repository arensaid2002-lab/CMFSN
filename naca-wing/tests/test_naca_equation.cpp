#include <iostream>
#include "nacawing/naca.hpp"
#include <cmath>

int main()
{
    double result = nacawing::y_t(0.12, 0.5);
    if ( std::abs(result - 0.05294) > 0.00001 ) { // The published NACA 0012 table gives 0.05294 at x = 0.5;  NACA012 from page 71  : https://ntrs.nasa.gov/citations/19930090976
        std::cerr << "y_t() returned incorrect value\n";
        return 1;
    }
    std::cout << "NACA equation test passed\n";
    return 0;
}