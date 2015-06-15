#include <string>
#include <iostream>

#include <vector>
#include <set>

std::vector<int> return_array(int size)
{
  std::vector<int> return_array;

  for(auto a=0;a<size;a++)
  {
    return_array.push_back(a);
  }
  return return_array;
}

class Bonds {
  std::vector<int> type;
  std::vector<int> id1;
  std::vector<int> id2;
public:
  Bonds(int, int, int);

};

int main(int argc, char *argv[])
{
  /* this does not work!
  int value = 5;
  int *another_value;
  *another_value = value;
  std::cout <<"\t"<< &another_value<<"\n";
  */
  /* checking arrays of string
  std::string *string_names;
  const int NUM = 20;
  string_names = new std::string[NUM];
  for (auto i=0; i<NUM;i++ )
    string_names[i]="abcd";
  for(auto ii=0;ii<NUM;ii++)
    std::cout << string_names[ii]<<"\n";
  */

  /* some string comparison tests
  char test[80]="";
  std::strcat(test, "abcd");
  std::string test1 ="efgh";

  std::string test2;
  test2.assign(test);
  std::cout << test<< " " << test1<< " " <<std::strcmp(test, test1.c_str()) << " is the answer\n";
  if((test1==test2))
    std::cout << " howdy there\n";
  std::cout << test1<< " " << test2<< " " << " is the answer\n";
  */

  /* vector stuff
  std::vector<int> test_vectors;
  test_vectors.push_back(10);
  std::cout << test_vectors.size()<<"\n";
  test_vectors.push_back(11);
  std::cout << test_vectors.size()<<"\n";
  //test_vectors[2]=12;
  std::cout << test_vectors.size()<<"\n";
  std::cout << test_vectors.size()<<"\t"<<test_vectors[0]<<"\t"<<test_vectors[1]<<"\t"<<test_vectors[2]<<"\n";
  for (auto it_vector = test_vectors.begin(); it_vector <=test_vectors.end(); it_vector++)
  {
    std::cout << (*it_vector)<<"\n";
  }
  auto return_vec = return_array(5);
  for(auto r_vec = return_vec.begin(); r_vec < return_vec.end(); r_vec++ )
    std::cout << (*r_vec)<<"\t";
  std::cout <<"\n";
  */
/* tests about set
  std::set<int> test_set;
  test_set.insert(10);
  test_set.insert(20);
  test_set.insert(30);
  for(auto it_set = test_set.begin(); it_set != test_set.end(); it_set++)
  {
    std::cout << (*it_set)<<"\t";
  }
  std::cout <<"\n";
  auto test_find = test_set.find(20);
  if(test_find == test_set.end())
    std::cout << "did not find it\n";
  else
    std::cout <<"found the value to be " << *test_find<<"\n";
  std::vector<int> test_vector(test_set.begin(), test_set.end());
  auto arr = test_vector.data();
  std::cout << sizeof(arr[0])<<"\n";
  std::cout <<sizeof(arr)/sizeof(arr[0])<<"\n";
  for (auto i=0; i<3;i++)
    std::cout << arr[i]<<" from an array\n";
tests about set */

/* an pointer to integer pointers */

typedef int bond_data [2];
using bond_data2 = int [2];
using ab = int;
std::vector<bond_data2*> abcd;
std::vector<int*> abcd2;
for (int i = 0; i <10;i++)
{
  int efgh[2];
  efgh[0] = i;
  efgh[1] = i+20;
  abcd2.push_back(efgh);
}
bond_data2 efgh2;
efgh2[0] = 20;
efgh2[1] = 30;
//abcd.push_back(&efgh2);
for(auto it=abcd2.begin(); it !=abcd2.end(); it++)
{
  std::cout << *it[0]<<"\n";
}
//std::cout << (*(abcd[0]))[1]<<'\n';
if(abcd.begin() != abcd.end()) abcd.clear();
std::cout << abcd.size()<<"\n";

}
