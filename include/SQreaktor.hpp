//
//  SQreaktor.hpp
//  
//  Utilities to transform the tensor contractions into the efficient computer code
//
//  Created by Masaaki Saitow on 13/05/27.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#pragma once

#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <Femto.hpp>
#include <SQtensor.hpp>
#include <SQindex.hpp>
#include <SQtensor.hpp>
#include <SQterm.hpp>
#include <SQcontract.hpp>
#include <SQbinary.hpp> 
#include <SQints.hpp>

namespace Femto { namespace Reaktor { 

    // *********************************************************
    // Returns name of the intermediate tensor
    // *********************************************************
    inline const std::string name_Int() 
    { return "X"; }
    
    // *********************************************************
    // Returns a sorted SharedIndices object
    // *********************************************************
    void reorderIndices(SharedIndices &inds);

    // *********************************************************
    // Compares two SQindices
    // *********************************************************
    bool compareIndices(SharedIndices const & inds1, SharedIndices const & inds2);

    inline long unsigned int bigNum () { return 100000000000000000; }

    //////////////////////////////////////////
    // Size of the target system
    inline long unsigned int numCore() { return 30;  }
    inline long unsigned int numAct()  { return 30;  }
    inline long unsigned int numVirt() { return 500; }
    //////////////////////////////////////////

    //////////////////////////////////////
    // Priority in binary decomposition
    enum Priority {Flops=0, Interm=1, ERI=2};
    //////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Divide the given tensor contractions into a series of the binary contractions
    void BinaryDecomposition(const SQcontracts &inContras, SQcont<SQbinaries> &outBins, Priority prior=ERI);
    // Small utility to estimate the size of the intermediate, or the flop count
    long unsigned int CalcOrder(const SQcont<SharedIndex> &Inds);
    //////////////////////////////////////////////////////////////////////////////////////////////////

}} // Femto::Reaktor
