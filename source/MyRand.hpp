//
//  MyRand.hpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 7/28/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#ifndef MyRand_hpp
#define MyRand_hpp

#include <random>

class MyRand
{
private:
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution;
  
public:
  
  // e.g. std::random_device r; seed = r()
  MyRand(unsigned int seed) : distribution(0.0,1.0), generator(seed) {}
  MyRand() : distribution(0.0,1.0) {}
  
  double u() {return distribution(generator);}
  
};

#endif /* MyRand_hpp */
