#ifndef RANDOM_H
#define RANDOM_H

class Random {
 public:
  int seed;

  Random();
  ~Random() {}
  void set_seed(int);
  void init() {}
  double uniform();
  double gaussian(double, double);
  double exponential(double);
};

#endif

