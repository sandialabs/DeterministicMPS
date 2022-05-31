//
//  TestChockSample.hpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 7/28/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#ifndef TestChockSample_hpp
#define TestChockSample_hpp

#include <vector>
#include <string>
#include <cmath>

void test_phi_sample();
void test_r_sample();

// create n samples
void sample_chock(int n, double phi, std::vector<double> &phis, std::vector<double> &rs);

// make a picture of sampling a chock, to visually see if it looks uniform by area
void test_uniform(double phi = M_PI_4, int num_samples = 4000, std::string fname = "test_uniform");

void chock_sample_figure();

#endif /* TestChockSample_hpp */
