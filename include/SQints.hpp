//
//  SQints.hpp
//  
//  Names of the molecular integrals are defined here
//
//  Created by Masaaki Saitow on 13/05/25.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#pragma once

namespace Femto{

  // *********************************************************
  // Returns name of the one-body integrals 
  // *********************************************************
  inline const std::string name_h1() 
  { return "h"; }

  // *********************************************************
  // Returns name of the two-body integrals
  // *********************************************************
  inline const std::string name_h2() 
  { return "V2"; }

  // *********************************************************
  // Returns name of the core Fock matrix of rank 0
  // *********************************************************
  inline const std::string name_cFock0() 
  { return "Fc0"; }

  // *********************************************************
  // Returns name of the core Fock matrix of rank 1
  // *********************************************************
  inline const std::string name_cFock1() 
  { return "Fc1"; }

  // *********************************************************
  // Returns name of the Fock matrix
  // *********************************************************
  inline const std::string name_Fock() 
  { return "F"; }

} // Femto
