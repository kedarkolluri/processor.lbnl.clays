#include <iostream>
using namespace std;

class Circle {
  double radius;
  public:
    Circle(double r) : radius(r) {};
    double area() {return 3.14159265*radius*radius;}
};

class Cylinder {
  Circle base;
  double height;

  public:
    Cylinder(double r, double h) : base (r), height(h) {};
    double volume() {return base.area()*height;}
};

int main () {
/*
  Cylinder foo(10.0, 20.0);
  cout << "foo's volume: "<< foo.volume() << '\n';
  return 0;
*/

}
