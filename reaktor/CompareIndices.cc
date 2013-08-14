//
//  ReorderIndices.cc
//  
//  Small utility to compare two SQindices
//
//  Created by Masaaki Saitow on 13/05/27.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <SQreaktor.hpp>

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor {

    // *********************************************************
    // Compares two SQindices
    // *********************************************************
    bool compareIndices(SharedIndices const & inds1, SharedIndices const & inds2)
    {
      SQindices temp1, temp2;
      for(auto ip = inds1.cbegin();ip != inds1.cend();++ip) temp1 <= **ip;
      for(auto ip = inds2.cbegin();ip != inds2.cend();++ip) temp2 <= **ip;
      stable_sort(temp1.begin(), temp1.end());
      stable_sort(temp2.begin(), temp2.end());
      return (temp1 == temp2);
    }

}} // Femto::Reaktor
