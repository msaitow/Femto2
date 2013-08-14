//
//  CalcOrder.cc
//  
//  Small utility to estimate the order of contraction, or size of the intermediate
//
//  Created by Masaaki Saitow on 13/05/27.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <SQreaktor.hpp>

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor {

    long unsigned int CalcOrder(const SharedIndices &Inds)
    {
      int nCore(0), nAct(0), nVirt(0);
      for(size_t num_i = 0;num_i < Inds.size();++num_i){
	if     (Inds[num_i]->get_char() == Femto::core) ++nCore;
	else if(Inds[num_i]->get_char() == Femto::act ) ++nAct ;
	else if(Inds[num_i]->get_char() == Femto::virt) ++nVirt;
      } // End num_i
      return pow(numCore(),nCore) * pow(numAct(), nAct) * pow(numVirt(), nVirt);							
    }

}} // Femto::Reaktor
