#include <string>
#include <iostream>
#include <vector>
#include <set>

typedef struct teststruct{
  int test1;
  int test2[3];
  int *test3;
  std::vector<int> test4;
};

int main(int argc, char *argv[])
{
  teststruct abcd;
  abcd.test1 = 2;
  abcd.test2[0] = 2;
  abcd.test2[1] = 3;
  abcd.test2[2] = 4;
  std::cout << abcd.test2[2]<<"\n";

  abcd.test3 = (int*) malloc(5*sizeof(int));
  abcd.test3[0] = 23;
  abcd.test3[1] = 24;
  abcd.test3[2] = 25;
  abcd.test3[3] = 26;
  abcd.test3[4] = 27;

  for (int i = 0 ; i  < 5; i++)
    std::cout << abcd.test3[i]<<" ";
  std::cout << "\n";
  return 0;
}
