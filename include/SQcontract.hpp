//
//  SQcontract.hpp
//  
//  Class that represents the tensor contraction, which is composed of a tensor on LHS 
//  and at least one tensor on RHS. And the associated indices of tensors on both LHS and
//  RHS are linked and stored as summedIndices_
//
//  Created by Masaaki Saitow on 14/4/10.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//


#pragma once

#include <map>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <SQindex.hpp>
#include <SQtensor.hpp>
#include <SQterm.hpp>
#include <Femto.hpp>


namespace Femto { namespace Core {

  // Class that represents the tensor contraction
  class SQcontract{
  
  public:
    SQcontract();
    SQcontract(const SQcontract &obj);
    //SQcontract(const SQbinary &obj);
    SQcontract(const double numConst, const SQcont<std::string> Consts,
               const SQtensor &Ltensor, const SQcont<SQtensor> &Rtensors);
    SQcontract(const SQtensor &Ltensor, const SQterm &Rterm);

    SQcontract operator=(const SQcontract &obj);  
    bool operator==(const SQcontract &obj) const;
    bool operator<(const SQcontract &obj) const;
    friend std::ostream &operator<<(std::ostream &os, SQcontract const &t);

    SQtensor get_Ltensor() const;
    SharedTensor get_Ltensor_ptr();
    SQcont<SQtensor> get_Rtensors() const;
    SQcont<SharedTensor> get_Rtensors_ptr();
    SQcont<SharedIndex> get_summedBody();
    double get_numConst() const;
    SQcont<std::string> get_Consts() const;

    // These three functions are for specifications on the internal indices
    // that defines the body of the Ltensor as a tensor quantities at C++/FORTRAN levels
    void set_Lindices(bool val);          
    void set_Rindices(int num, bool val);
    SQcont<SharedIndex> get_Lindices() const;       
    std::map<SQtensor, SQcont<SharedIndex> > get_Rindices() const; // <SQtensor, its indices>

    void set_Ltensor(const SQtensor &Tensor);
    void set_Rtensors(const SQcont<SQtensor> &Tensors);
    void set_numConst(const double num);
    void set_Consts(const SQcont<std::string> consts);
    void masquerade();
    void print_summedBody();
    void contractkDeltas();

  protected:
    void set_summedBody();

    double   numConst_;               // Numerical factor
    SQcont<std::string> Consts_; // Constants of this contraction
    //*-- --*// SQcont<std::string> InternalNames_;   // Names of indices of the FORTRAN level representation of Ltensor_,
    //*-- --*//                                  // if nothing is set to this, internal indices are used

    bool Lindices_;             // If this is set true, all the indices of Ltensor_ becomes the indices used at the FORTRAN level. 
                                // By default this is set to false, which means only internal indices are used to allocate Ltensor_ at the FORTRAN level
    SQcont<bool> Rindices_;     // If this is set true, all the indices of Rtensors_ becomes the indices used at the FORTRAN level. 
                                // By default this is set to false, which means only internal indices are used to allocate Rtensors_ at the FORTRAN level

    SQtensor Ltensor_;               // Tensor on the left-hand side
    SQcont<SQtensor> Rtensors_;      // Tensors on the right-hand side
    SQcont<SharedIndex> summedIndices_;  // Body of summed indices, of which shared_ptrs are shared
                                         // among the Tensors_ as vector of the SQindex(not as pointer)
  private:
    // Some stuffs for serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int ver){
      ar & boost::serialization::make_nvp("numConst",       numConst_);
      ar & boost::serialization::make_nvp("Consts",         Consts_);
      ar & boost::serialization::make_nvp("Ltensor",        Ltensor_);
      ar & boost::serialization::make_nvp("Rtensors",       Rtensors_);
      ar & boost::serialization::make_nvp("Lindices_",      Lindices_);
      ar & boost::serialization::make_nvp("Rindices_",      Rindices_);
      ar & boost::serialization::make_nvp("Rtensors",       Rtensors_);
      ar & boost::serialization::make_nvp("summedIndices_", summedIndices_);
    }

  };

//old//   // Class that represents the binary contraction
//old//   class SQbinary : public SQcontract{
//old//   public:
//old//     SQbinary();
//old//     //SQbinary(const SQbinary &obj);
//old//     SQbinary(const double numConst, const SQcont<std::string> Consts,
//old//              const SQtensor &Ltensor, const SQcont<SQtensor> &Rtensors);
//old//     
//old//     //SQbinary operator=(const SQbinary &obj);  
//old//   private:
//old//     // Some stuffs for serialization
//old//     friend class boost::serialization::access;
//old//     template<class Archive>
//old//     void serialize(Archive & ar, const unsigned int ver){
//old//       ar & boost::serialization::base_object<SQtensor>(*this);
//old//     }     
//old//   };
  
  // Whether two tensor contractions are factorizable or not
  bool isFactorizable(SQcontract &a, SQcontract &b);

}} // Femto::SQcontract

////////////////////////// Abbreviations ///////////////////////
typedef Femto::SQcont<Femto::Core::SQcontract> SQcontracts;
//typedef Femto::SQcont<Femto::Core::SQbinary> SQbinaries;
////////////////////////////////////////////////////////////////

