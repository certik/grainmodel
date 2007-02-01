%module libmeshpy
%{
#include <iostream>
int pr()
{
    std::cout << "heja!";
    return 5;
}
%}

int pr();
