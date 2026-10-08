#include <iostream>
#include "nacawing/NACA.hpp"
#include <cmath>

int main()
{

    int failures = 0;

    double thickness_0012 = nacawing::y_t(0.12, 0.5);
    if ( std::abs(thickness_0012 - 0.05294) > 0.00001 ) { // The published NACA 0012 table gives 0.05294 at x = 0.5;  NACA0012 from page 71  : https://ntrs.nasa.gov/citations/19930090976
        std::cerr << "thickness_0012 y_t() returned incorrect value\n";
        failures++;
    } 
    
    else {
        std::cout << "thickness_0012 y_t() returned correct value\n";
    }


    double camber_0005 = nacawing::y_c(0.0, 0.0, 0.5);
    if ( std::abs(camber_0005 - 0.0) > 0.00001 ) { // The result should be 0.00
        std::cerr << "camber_0005 y_c() = 0 test returned incorrect value\n";
        failures++;
    }
    else {
        std::cout << "camber_0005 y_c() = 0 test returned correct value\n";
    }


    double camber_2404 = nacawing::y_c(0.02, 0.4, 0.4);
    if ( std::abs(camber_2404 - 0.02) > 0.00001 ) { // The result should be 0.02
        std::cerr << "camber_2404 y_c() returned incorrect value\n";
        failures++;
    }
    else {
        std::cout << "camber y_c() returned correct value\n";
    }


    return failures == 0 ? 0 : 1;

   

}