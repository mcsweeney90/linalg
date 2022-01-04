#include "linalg.hpp"
#include <iostream>

int main(int argc, char const *argv[]) {
  MyMatrix A = Poisson(8);
  std::cout << A << std::endl;
  return 0;
}
