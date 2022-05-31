//
//  TestSimpleMPS.hpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 12/8/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#ifndef TestSimpleMPS_hpp
#define TestSimpleMPS_hpp

void SimpleMPS_sampling(int gridm, int gridn, bool periodic=true, unsigned int rand_seed_p=345038235);

void MPS0();
void MPStime();
void MPStimelong();
void MPSscaling(bool periodic, int maxt = 12);

void SimpleMPS_runtime_table();

void SimpleMPS_stats();


#endif /* TestSimpleMPS_hpp */
