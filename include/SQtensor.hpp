//
//  SQtensor.hpp
//  
//  Class that represents the tensorial object, such as spin-free generator,
//  the amplitudes and so on. 
//
//  Created by Masaaki Saitow on 13/04/10.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//


#pragma once

#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <boost/serialization/vector.hpp>
#include <Femto.hpp>
#include <SQindex.hpp>
#include <SQcont.hpp>

//#define _ABORT

namespace Femto { namespace Core {

  // Class that represents the tensor object
  class SQtensor{
  
  public:
    SQtensor();
    SQtensor(const SQtensor &obj);
    SQtensor(const std::string &name, const SQcont<SharedIndex> indices, 
             const Symmetry &symm, const bool isCommutable=true, const notation nota=Dirac);

//*OLD*     SQtensor(const std::string &name, const SQcont<SharedIndex> indices, 
//*OLD*              const Symmetry &symm, const bool isCommutable=true);
  
    std::string get_name() const;
    std::string get_notation_str() const;
    notation get_notation() const;
    SQcont<SharedIndex> get_indices() const;
    //  bool isAdditive(const SQtensor &obj) const;
    bool isCommutable() const;
    bool hasIndex(const SQindex *i) const;
    SQtensor &operator=(const SQtensor &obj); 
    bool operator<(const SQtensor &obj) const;
    bool operator>(const SQtensor &obj) const;
    bool operator==(const SQtensor &obj) const; // Should be improved considering all the permutational symetry
    IIvector get_perms() const;
    //    void makeAllPerms();
    int sortIndices();
    void rotateIndices(const size_t i);
    friend std::ostream &operator<<(std::ostream &os, const SQtensor &s);
    void put_indices(size_t i, SharedIndex v);
    void print_symm() const; // Print permutational symmetry
    void convertD2M();       // Convert all the indices and permutational symmetry from Dirac into Mulliken order
    void convertM2D();       // Convert all the indices and permutational symmetry from Mulliken into Dirac order
    std::string convert2LaTeX() const; // Convert to LaTeX
  
  protected:
    void makeAllPerms();

    //int sign_;                 // Sign of this tensor
    std::string name_;         // Name of tensor
    notation not_;             // Notation of tensor
    bool isCommutable_;        // Flag that shows this tensor is commutable
    SQcont<SharedIndex> indices_; // Vector of the associated indices
    Symmetry symm_;            // Permutational symmetry object
    IIvector permutations_;    // All patterns of permutation
    Ivector  factors_;         // All patterns of factors

  private:
    // Some stuffs for serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int ver){
      ar & boost::serialization::make_nvp("name",         name_);
      ar & boost::serialization::make_nvp("not",          not_);
      ar & boost::serialization::make_nvp("isCommutable", isCommutable_);
      ar & boost::serialization::make_nvp("indices",      indices_);
      ar & boost::serialization::make_nvp("symm",         symm_);
      ar & boost::serialization::make_nvp("permutations", permutations_);
      ar & boost::serialization::make_nvp("factors",      factors_);
    }

//*TEST*     // Some stuffs for serialization
//*TEST*     friend class boost::serialization::access;
//*TEST*     template<class Archive>
//*TEST*     void serialize(Archive & ar, const unsigned int ver){
//*TEST*       ar & name_;
//*TEST*       ar & not_;
//*TEST*       ar & isCommutable_;
//*TEST*       ar & indices_;
//*TEST*       ar & symm_;
//*TEST*       ar & permutations_;
//*TEST*       ar & factors_;
//*TEST*     }

  };
  
  // Class that represents the Kronecker's delta object
  class kDelta : public SQtensor{
  public:
    kDelta();
    kDelta(const SQcont<SharedIndex> indices);
  private:
    // Some stuffs for serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int ver){
      ar & boost::serialization::base_object<SQtensor>(*this);
    }    
  };

  // Class that represents the Spin-free unitary group generator object
  class sfGen : public SQtensor{
  public:
    //    sfGen();
    sfGen(const SQcont<SharedIndex> indices);
    //private:
    //int order_; // Order of sfGen
  private:
    // Some stuffs for serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int ver){
      ar & boost::serialization::base_object<SQtensor>(*this);
    }    
  };

  // Class that represents the reduced-density matrix (RDM) object
  class RDM : public SQtensor{
  public:
    //    RDM();
    RDM(const SQcont<SharedIndex> indices);
    //private:
    //int order_; // Order of RDM
  private:
    // Some stuffs for serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int ver){
      ar & boost::serialization::base_object<SQtensor>(*this);
    }    
  };

  // Class that represents the creation operator object
  class aCre : public SQtensor{
  public:
    aCre(SharedIndex index);
  private:
    // Some stuffs for serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int ver){
      ar & boost::serialization::base_object<SQtensor>(*this);
    }    
  };

  // Class that represents the destruction operator object
  class aDes : public SQtensor{
  public:
    aDes(SharedIndex index);
  private:
    // Some stuffs for serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int ver){
      ar & boost::serialization::base_object<SQtensor>(*this);
    }    
  };

  // Judge just the two tensors are same without rotating the indices
  bool are_sameforms(const SQtensor &ta, const SQtensor &tb);

}} // Femto::core

////////////////////////// Abbreviations ///////////////////////
// Acronym for SQtensor* and this would be changed in future
typedef Femto::Core::SQtensor* SharedTensor;
typedef Femto::SQcont<Femto::Core::SQtensor> SQtensors;
typedef Femto::SQcont<Femto::Core::SQtensor*> SharedTensors;
////////////////////////////////////////////////////////////////
