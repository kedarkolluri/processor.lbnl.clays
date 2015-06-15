#include <iostream>
#include <vector>
using namespace std;

int main()
{
  std::vector<int> myvector;
  for(auto i=0;i<10; i++) myvector.push_back(i);

  cout << myvector[5]<<" is value of 5 position\n";

  for (auto i = myvector.begin(); i!=myvector.end(); i++ )
  {
    cout << *i <<" ";
  }
  cout <<"\n";
  myvector.erase(myvector.begin()+5);
  cout << "after erasing \n";
  for (auto i = myvector.begin(); i!=myvector.end(); i++ )
  {
    cout << *i << " ";
  }
  cout <<"\n";

  cout << myvector[5]<<" is value of 5 position\n";

  myvector.erase(myvector.begin()+5);
  cout << "after erasing \n";
  for (auto i = myvector.begin(); i!=myvector.end(); i++ )
  {
    cout << *i << " ";
  }
  cout <<"\n";
  


}
