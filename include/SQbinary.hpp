//
//  SQbinary.hpp
//  
//  Class that represents the binary contraction, which is composed of a tensor on LHS 
//  and at least one tensor on RHS. Unlikely to SQcontract, SQbinary contains up to only
//  two tensors on RHS 
//
//  Created by Masaaki Saitow on 13/04/10.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//


#pragma once

#include <SQcontract.hpp>


namespace Femto { namespace Core {

  // Class that represents the binary contraction
  class SQbinary : public SQcontract{
  public:
    SQbinary();
    SQbinary(const SQbinary &obj);
    SQbinary(const double numConst, const SQcont<std::string> Consts,
             const SQtensor &Ltensor, const SQcont<SQtensor> &Rtensors);
    SQbinary(const SQtensor &Ltensor, const SQterm &Rterm);
    
    //SQbinary operator=(const SQbinary &obj);  
  private:
    // Some stuffs for serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int ver){
      ar & boost::serialization::base_object<SQtensor>(*this);
    }     
  };

}} // Femto::SQbinary

////////////////////////// Abbreviations ///////////////////////
typedef Femto::SQcont<Femto::Core::SQbinary> SQbinaries;
////////////////////////////////////////////////////////////////
