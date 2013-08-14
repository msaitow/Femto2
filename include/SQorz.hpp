//
//  SQorz.hpp
//  
//  Utilities to generate the tensor contraction code, which is capable of working with orz suite
//
//  Created by Masaaki Saitow on 13/06/28.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#pragma once

#include <SQreaktor.hpp>
#include <SQnetwork.hpp>

namespace Femto { namespace Reaktor { namespace Orz {

      // *********************************************************
      // Rotate the indices of the tensor, and set the loading indices as external (optionally)
      // *********************************************************
      void SetIndicesNet(ContraNet &InNet, const bool do_rotate=true);    

}}} // Femto::Reaktor::Orz
