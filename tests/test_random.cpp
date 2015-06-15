#include <iostream>
#include <ctime>
using namespace std;

int main()
{
  std::srand(time(NULL));

  std::cout << rand() << std::endl;

}
