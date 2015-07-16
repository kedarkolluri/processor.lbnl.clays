#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <tuple>
#include <map>
#include <typeinfo>

typedef struct teststruct{
  int test1;
  int test2[3];
  int *test3;
  std::vector<int> test4;
};

typedef struct test{
  int test1;
  int test2;
};

typedef struct toptest{
  int test0;
  std::map<int, test> testmap;
};

toptest filltoptest()
{
  int n=100;
  toptest tt1;
  tt1.test0=100;
  for (int i=0;i<100;i++)
  {
    tt1.testmap[i].test1 = i*i;
    tt1.testmap[i].test2 = i*i*i;
  }
  return tt1;
}

void execute_one()
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

  std::map<std::tuple<int, int>, double> map_double;
  auto tpl = std::make_tuple(1,1);
  map_double[tpl] = 23.2;
  map_double[std::make_tuple(1,2)] = 12.434;
}

void test_tuplesinmap()
{
  std::map< std::tuple<int, int>, int> testmap1;
  auto abcd = std::make_tuple(1, 3);
  testmap1[abcd] = 1;
  std::cout << testmap1[abcd]<<"\n";
  abcd = std::make_tuple(3, 1);
  testmap1[abcd]= 10;
  std::cout << testmap1[abcd]<<"\n";
}

void test_loop_map()
{
  //create a map and fill it
  std::map<int, test> looptestmap;

  for(int i =0;i<10;i++)
  {
    test testtest;
    testtest.test1 = i*i;
    testtest.test2  = i*i*i;
    looptestmap[i] = testtest;
  }
  // retrieve it back again

  for(auto iter = looptestmap.begin(); iter != looptestmap.end(); ++iter)
  {
    test localtest = iter->second;
    std::cout << iter->first<<" "<<localtest.test1 << " "<< localtest.test2<<" \n";
  }

}
void test_tuple_vector()
{
  std::vector< std::tuple<int, int> > tuplevector;
  tuplevector.push_back(std::make_tuple(1, 1));
  tuplevector.push_back(std::make_tuple(2, 2));
  tuplevector.push_back(std::make_tuple(3, 2));
  tuplevector.push_back(std::make_tuple(4, 2));

  for(auto iter = tuplevector.begin(); iter !=tuplevector.end(); ++iter)
  {
    std::cout << typeid(iter).name()<<"\t"<< typeid(*iter).name()<<"\n";
    std::cout << std::get<0>(*iter)<<"\n";
  }
}
int main(int argc, char *argv[])
{
  //test_loop_map();
  test_tuple_vector();

  return 0;
}
