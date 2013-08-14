//
//  ReorderIndices.cc
//  
//  Small utility to reorder the SharedIndices
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
    // Returns a sorted SharedIndices object
    // *********************************************************
    void reorderIndices(SharedIndices &inds)
    {
      SQindices temp;
      SharedIndices outs;
      for(auto ip = inds.cbegin();ip != inds.cend();++ip) temp <= **ip;
      sort(temp.begin(), temp.end());
      for(auto i = temp.begin();i != temp.end();++i){
	for(auto j = inds.cbegin();j != inds.cend();++j){
          if(*i==**j) outs <= *j;
	} // End j
      } // End i
      inds = outs;
    }

}} // Femto::Reaktor
