//
//  SQcontract.cc
//  
//
//  Created by Masaaki Saitow on 12/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <cmath>
#include <vector>
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <Femto.hpp>
#include <SQcontract.hpp>
#include <SQbinary.hpp>
#include <SQtensor.hpp>
#include <SQterm.hpp>

#define _NEW_COMBINATION

using namespace std;

namespace Femto { namespace Core {

  // *********************************************************
  // 
  // *********************************************************
  SQcontract::SQcontract()
  {}

  // *********************************************************
  // 
  // *********************************************************
  SQcontract::SQcontract(const double numConst, const SQcont<string> Consts, 
		         const SQtensor &Ltensor, const SQcont<SQtensor> &Rtensors)
    : numConst_(numConst),
      Consts_(Consts),
      Ltensor_(Ltensor),
      Rtensors_(Rtensors),
      Lindices_(false)
  {   
    if(is_RDM(Ltensor_.get_name()) || is_sfGen(Ltensor_.get_name()) || Ltensor_.get_name() == kDelta_name())
      { cout << "SQcontract: Ltensor cannot be either of sfGen, RDM or kDelta" << endl; abort(); }
    for(auto t = Rtensors_.begin();t != Rtensors_.end();++t){
      Rindices_ <= false;
      if(is_sfGen(t->get_name())){
        cout << "Spin-free unitary group generator is detected: " << endl;
        abort();
      } // End if
    } // End t
    set_summedBody();
  }

  // *********************************************************
  // 
  // *********************************************************
  SQcontract::SQcontract(const SQtensor &Ltensor, const SQterm &Rterm)
    : numConst_(Rterm.get_numConst()),
      Consts_(Rterm.get_Consts()),
      Ltensor_(Ltensor),
      Rtensors_(Rterm.get_tensors()),
      Lindices_(false)
  {
    if(is_RDM(Ltensor_.get_name()) || is_sfGen(Ltensor_.get_name()) || Ltensor_.get_name() == kDelta_name())
      { cout << "SQcontract: Ltensor cannot be either of sfGen, RDM or kDelta" << endl; abort(); }
    for(auto t = Rtensors_.begin();t != Rtensors_.end();++t){
      Rindices_ <= false;
      if(is_sfGen(t->get_name())){
        cout << "Spin-free unitary group generator is detected: " << endl;
        abort();
      } // End if
    } // End t
    set_summedBody();    
  }

  // *********************************************************
  // 
  // *********************************************************
  SQcontract::SQcontract(const SQcontract &obj)
    : Ltensor_(obj.Ltensor_),
      numConst_(obj.numConst_),
      Consts_(obj.Consts_),
      Lindices_(obj.Lindices_),
      Rindices_(obj.Rindices_),
      Rtensors_(obj.Rtensors_)
  { set_summedBody(); }

//*   // *********************************************************
//*   // 
//*   // *********************************************************
//*   SQcontract::SQcontract(const SQbinary &obj)
//*     : Ltensor_(obj.Ltensor_),
//*       numConst_(obj.numConst_),
//*       Consts_(obj.Consts_),
//*       Lindices_(obj.Lindices_),
//*       Rindices_(obj.Rindices_),
//*       Rtensors_(obj.Rtensors_)
//*   { set_summedBody(); }

  // *********************************************************
  // 
  // *********************************************************
  SQcontract SQcontract::operator=(const SQcontract &obj){
    numConst_      = obj.numConst_;
    Consts_        = obj.Consts_;
    Ltensor_       = obj.Ltensor_;
    Rtensors_      = obj.Rtensors_;
    Lindices_      = obj.Lindices_;    
    Rindices_      = obj.Rindices_;    

    set_summedBody();
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::set_summedBody()
  {

    // Copy all the indices found in Rtensors_ into allIndices
    SQcont<SharedIndex> allIndices;
    for(auto t = Rtensors_.begin();t != Rtensors_.end();++t){
      for(size_t inum = 0;inum < t->get_indices().size();++inum){
	bool is_OK(true);
        for(auto I = allIndices.begin();I != allIndices.end();++I){
          if(**I==*t->get_indices()[inum]) is_OK = false;
	} // End I
        if(is_OK){
          allIndices <= SharedIndex( new SQindex(t->get_indices()[inum]->get_index(),
                                                 t->get_indices()[inum]->get_char(),
                                                 t->get_indices()[inum]->get_isSummed(),
						 t->get_indices()[inum]->get_isExt()));
	} // End if
      } // End inum 
    } // End t

    // Then, copy all the indices of the LTensor
    for(size_t inum = 0;inum < Ltensor_.get_indices().size();++inum){
      bool is_OK(true);
      for(auto I = allIndices.begin();I != allIndices.end();++I){
	if(**I==*Ltensor_.get_indices()[inum]) is_OK = false;
      } // End I
      if(is_OK){
	allIndices <= SharedIndex( new SQindex(Ltensor_.get_indices()[inum]->get_index(),
					       Ltensor_.get_indices()[inum]->get_char(),
					       Ltensor_.get_indices()[inum]->get_isSummed(),
					       Ltensor_.get_indices()[inum]->get_isExt()));
      } // End if      
    } // End inum

    // Link indices in Rtensors_ to allIndices
    for(auto t = Rtensors_.begin();t != Rtensors_.end();++t){
      for(size_t inum = 0;inum < t->get_indices().size();++inum){
	for(auto j = allIndices.begin();j != allIndices.end();++j)
	  if(*t->get_indices()[inum] == **j) t->put_indices(inum, *j);
      } // End inum
    } // End t

    // Link indices in Ltensor_ to allIndices
    for(size_t inum = 0;inum < Ltensor_.get_indices().size();++inum){
      for(auto j = allIndices.begin();j != allIndices.end();++j)
	if(*Ltensor_.get_indices()[inum] == **j) Ltensor_.put_indices(inum, *j);
    } // End inum

    summedIndices_.clear();
    summedIndices_.reserve(allIndices.size());

    // Copy the content of allIndices to summedIndices_
    for(auto i = allIndices.begin();i != allIndices.end();++i)
      summedIndices_ <= SharedIndex( new SQindex((*i)->get_index(),
						 (*i)->get_char(),
						 (*i)->get_isSummed(),
						 (*i)->get_isExt()));
    // Then re-link to summedIndices_
    for(auto t = Rtensors_.begin();t != Rtensors_.end();++t){
      for(size_t inum = 0;inum < t->get_indices().size();++inum){
	for(auto j = summedIndices_.begin();j != summedIndices_.end();++j)
	  if(*t->get_indices()[inum] == **j) t->put_indices(inum, *j);
      } // End inum
    } // End t

    for(size_t inum = 0;inum < Ltensor_.get_indices().size();++inum){
      for(auto j = summedIndices_.begin();j != summedIndices_.end();++j)
	if(*Ltensor_.get_indices()[inum] == **j) Ltensor_.put_indices(inum, *j);
    } // End inum

    ////////////////////////////////////////////////////////////////////////////
    // // Firstly, assign all the &summedIndices_[:] to &tempIndices[:]
    // vector<SQindex> tempIndices(summedIndices_);
    // for(size_t i = 0;i < tempIndices.size();++i){
    //   for(size_t I = 0;I < Rtensors_.size();++I){
    //     vector<SQindex*> temp1(Rtensors_[I].get_indices());
    //     for(size_t j = 0;j < temp1.size();++j){
    //       if(*temp1[j]==tempIndices[i]) Rtensors_[I].put_indices(j, &tempIndices[i]);
    // 	} // End j
    //   } // End I
    //   vector<SQindex*> temp1(Ltensor_.get_indices());
    //   for(size_t j = 0;j < temp1.size();++j){
    //     if(*temp1[j]==tempIndices[i]) Ltensor_.put_indices(j, &tempIndices[i]);
    //   } // End j
    // } // End i
    // summedIndices_.clear();
    // 
    // // Search for the dummy indices .... 
    // // Maybe it's nice for the sigma construction, but not for the diagonal preconditionor (2012/10/27)
    // for(size_t I = 0;I < Rtensors_.size();++I){
    //   vector<SQindex*> temp1(Rtensors_[I].get_indices());
    //   vector<SQindex> indices;
    //   for(size_t i = 0;i < temp1.size();++i) indices.push_back(*temp1[i]);
    //   for(size_t j = 0;j < indices.size();++j){
    //     if(find(summedIndices_.begin(), summedIndices_.end(), indices[j])==summedIndices_.end()) {
    //       summedIndices_.push_back(*temp1[j]);
    //     } // End if
    //   } // End j
    // } // End I
    // vector<SQindex*> temp1(Ltensor_.get_indices());
    // vector<SQindex> indices;
    // for(size_t i = 0;i < temp1.size();++i) indices.push_back(*temp1[i]);
    // for(size_t j = 0;j < indices.size();++j){
    //   if(find(summedIndices_.begin(), summedIndices_.end(), indices[j])==summedIndices_.end()) {
    //     summedIndices_.push_back(*temp1[j]);
    //   } // End if
    // } // End j
    // 
    // // Replace all the dummy indices, which are shared among all the terms, with those
    // // copied in the summedIndices
    // for(size_t i = 0;i < summedIndices_.size();++i){
    //   SQindex temp_i(summedIndices_[i]);
    //   for(size_t I = 0;I < Rtensors_.size();++I){
    //     vector<SQindex*> temp1 = Rtensors_[I].get_indices();
    //     vector<SQindex> indices;        
    //     for(size_t j = 0;j < temp1.size();++j) indices.push_back(*temp1[j]);
    //     for(size_t j = 0;j < indices.size();++j){
    //       if(indices[j] == temp_i) Rtensors_[I].put_indices(j, &summedIndices_[i]);
    // 	}
    //   } // End I
    //   for(size_t num_i = 0;num_i < Ltensor_.get_indices().size();++num_i)
    //     if(*(Ltensor_.get_indices()[num_i]) == temp_i) Ltensor_.put_indices(num_i, &summedIndices_[i]);
    // } // End i
    // 
    // //masquerade();
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::masquerade()
  {
     // Set the dummy names for all the dummy indices
     int Ccount(0);
     int Ocount(0);
     int Vcount(0);
     for(size_t i = 0;i < summedIndices_.size();++i){
       // In case of core
       if(summedIndices_[i]->get_char() == core && summedIndices_[i]->get_isSummed()){
         ++Ccount;
         ostringstream stm;
         stm << Ccount;
         summedIndices_[i]->put_index("c" + stm.str());
       }
       else
       // In case of active
 	if(summedIndices_[i]->get_char() == act && summedIndices_[i]->get_isSummed()){
         ++Ocount;
         ostringstream stm;
         stm << Ocount;
         summedIndices_[i]->put_index("o" + stm.str());
       }
       else
       // In case of virtual
 	if(summedIndices_[i]->get_char() == virt && summedIndices_[i]->get_isSummed()){
         ++Vcount;
         ostringstream stm;
         stm << Vcount;
         summedIndices_[i]->put_index("v" + stm.str());
       }     
     } // End i

  }

  // *********************************************************
  // 
  // *********************************************************
  SQcont<SQtensor> SQcontract::get_Rtensors() const
  { return Rtensors_; }

  // *********************************************************
  // 
  // *********************************************************
  SQtensor SQcontract::get_Ltensor() const
  { return Ltensor_; }

  // *********************************************************
  // 
  // *********************************************************
  SQcont<SharedTensor> SQcontract::get_Rtensors_ptr()
  {
    SQcont<SharedTensor> retval; retval.reserve(Rtensors_.size());
    for(size_t i = 0;i < Rtensors_.size();++i) retval.push_back(&(Rtensors_[i])); 
    return retval; 
  }

  // *********************************************************
  // 
  // *********************************************************
  SharedTensor SQcontract::get_Ltensor_ptr()
  { return &Ltensor_; }

  // *********************************************************
  // 
  // *********************************************************
  double SQcontract::get_numConst() const
  { return numConst_; }

  // *********************************************************
  // 
  // *********************************************************
  SQcont<string> SQcontract::get_Consts() const
  { return Consts_; }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::set_Rtensors(const SQcont<SQtensor> &Tensors)
  {
    Rtensors_ = Tensors;
    for(size_t i = 0;i < Rtensors_.size();++i) Rindices_.p()[i] = false;
    set_summedBody();    
  }

  //*-- --*// // *********************************************************
  //*-- --*// // 
  //*-- --*// // *********************************************************
  //*-- --*// vector<SQindex*> SQcontract::get_Lindices()
  //*-- --*// {
  //*-- --*//   vector<SQindex*> retval;
  //*-- --*//   // If InternalNames are not set, just return the internal indices.
  //*-- --*//   if(!InternalNames_.size()){
  //*-- --*//     for(size_t num_i = 0;num_i < summedIndices_.size();++num_i)
  //*-- --*//       if(!summedIndices_[num_i].get_isExt()) retval.push_back(&summedIndices_[num_i]);
  //*-- --*//   }
  //*-- --*//   // If InternalNames are given specifically, return only these indices.
  //*-- --*//   {
  //*-- --*//     for(size_t num_i = 0;num_i < summedIndices_.size();++num_i)
  //*-- --*//       if(find(InternalNames_.begin(), InternalNames_.end(), summedIndices_[num_i].get_index()) != InternalNames_.end()) 
  //*-- --*// 	  retval.push_back(&summedIndices_[num_i]);      
  //*-- --*//   }
  //*-- --*//   return retval;
  //*-- --*// }

  //*-- --*// // *********************************************************
  //*-- --*// // 
  //*-- --*// // *********************************************************
  //*-- --*// void SQcontract::set_Lindices(vector<string> &name_list)
  //*-- --*// { InternalNames_ = name_list; }
  //*-- --*// 
  //*-- --*// // *********************************************************
  //*-- --*// // 
  //*-- --*// // *********************************************************
  //*-- --*// void SQcontract::clear_Lindices()
  //*-- --*// { InternalNames_.clear(); }

  // *********************************************************
  // 
  // *********************************************************
  // NOTE :: 
  //   Basically, Ltensor_ can be sigma vector, preconditionor or the intermediate tensor.
  //   Lindices_ are intended to be used to represent how the Ltensor_ is treated in 
  //   the code generation step, which means, in case that Ltensor_ is an intermediate type
  //   tensor associated with the loading index of ERI, it's sufficient to use only the internal
  //   indices as the actual indices of DTensor, or symblock. But in the other case, in which, 
  //   Ltensor_ doesn't have the loading index, all the indices_ may have to be treated as the 
  //   actual indices of the DTensor representation of the Ltensor_. Lindices_ can stand for which types
  //   of the treatment is appropriate for LTensor_. In case of false (by default), only the internal indices
  //   are treated as the actual indices of DTensor representation. But if it's true, all the indices are used.  
  void SQcontract::set_Lindices(bool val)
  { Lindices_ = val; }

  // *********************************************************
  // 
  // *********************************************************
  // NOTE ::
  //   Meaning of Rindices_ is same to Lindices_ for LTensor_. For n-th member of Rtensors_, 
  //   Rindices_[n] stands for how the associated indices are treated in the code generating steps.
  void SQcontract::set_Rindices(int num, bool val)
  { Rindices_.p()[num] = val; }

  // *********************************************************
  // 
  // *********************************************************
  SharedIndices SQcontract::get_Lindices() const 
  { 
    SharedIndices inds;
    if(Lindices_) inds = Ltensor_.get_indices();
    else 
      for(size_t num_i = 0;num_i < Ltensor_.get_indices().size();++num_i)
	if(!Ltensor_.get_indices()[num_i]->get_isExt()) inds <= Ltensor_.get_indices()[num_i];  
    return inds;
  } 

  // *********************************************************
  // 
  // *********************************************************
  map<SQtensor, SharedIndices> SQcontract::get_Rindices() const
  { 
    map<SQtensor, SharedIndices> t_inds;
    for(size_t num_t = 0;num_t < Rtensors_.size();++num_t){
      SharedIndices inds(Rtensors_[num_t].get_indices());
      if(Rindices_.p()[num_t]){
        for(size_t num_i = 0;num_i < Rtensors_[num_t].get_indices().size();++num_i)
	  if(!Rtensors_[num_t].get_indices()[num_i]->get_isExt())
	    inds <= Rtensors_[num_t].get_indices()[num_i];
      } // End if
      else inds = Rtensors_[num_t].get_indices();
      t_inds.insert(map<SQtensor, SharedIndices>::value_type(Rtensors_[num_t], inds));
    } // End num_t

    return t_inds; 
  }  

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::set_Ltensor(const SQtensor &Tensor)
  {
    Ltensor_ = Tensor;
    set_summedBody();    
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::set_numConst(const double num)
  { numConst_ = num; }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::set_Consts(const SQcont<string> consts)
  { Consts_ = consts; }

  // *********************************************************
  // 
  // *********************************************************
  ostream& operator <<(std::ostream &os, SQcontract const &t)
  {
    os << t.get_Ltensor() << "+= ";
    os << boost::format("(%14.8f)") % t.get_numConst();
    if(t.get_Consts().size())
      os << " " <<= t.get_Consts();
    if(t.get_Rtensors().size())
      os << " " <<= t.get_Rtensors();       

    //*OLD* os << boost::format("(%14.8f)") % t.get_numConst() << " ";
    //*OLD* SQcont<string> strs(t.get_Consts());
    //*OLD* for(size_t i = 0;i < strs.size();++i) os << strs[i] << " ";
    //*OLD* 
    //*OLD* SQcont<SQtensor> tensors(t.get_Rtensors());    
    //*OLD* for(size_t i = 0;i < tensors.size();++i) os << tensors[i];

    return os;
  }

  // *********************************************************
  // 
  // *********************************************************
  SQcont<SharedIndex> SQcontract::get_summedBody()
  {
    SQcont<SharedIndex> retval;
    for(auto ind = summedIndices_.begin();ind != summedIndices_.end();++ind) 
      retval <= *ind;
    return retval;
  }

  // *********************************************************
  // Contract Kronecker deltas
  // *********************************************************
  void SQcontract::contractkDeltas()
  {

    // *******************************************************************
    // * The case that Kronecker's delta has at least one dummy index is *
    // * not considered because this class is not intended for the deri- *
    // * vation of the many-body equation                                *
    // *******************************************************************
    for(auto t = Rtensors_.begin(); t != Rtensors_.end();){
      // In case of delta_{p}^{p} ....
      if((t->get_name()==kDelta_name()) && (*(t->get_indices()[0]) == *(t->get_indices()[1])))
        t = Rtensors_.erase(t);
//*- -*//      else if(t->get_name()==kDelta_name() && t->get_indices()[0]->get_isSummed() || t->get_indices()[1]->get_isSummed()){
//*- -*//	cout << "SQcontract: kDelta with dummy index is detected. It's not the one to be handled at code-generation process" << endl;
//*- -*//	abort();
//*- -*//      } // End if
      // In case of both indices of kDelta are not dummy ....
      else if((t->get_name()==kDelta_name()) && !(t->get_indices()[0]->get_isSummed()) && !(t->get_indices()[1]->get_isSummed())){
        if(t->get_indices()[0]->get_char() != t->get_indices()[1]->get_char()){
	  numConst_ = 0; // This term shoud be zero
	  return;
	} // End if
	// If kdelta is composed if purely, non-dummy indices (not implemented in SQterm) ... 
	else{
	  SharedIndex killed;
	  SharedIndex killer;
	  if(*(t->get_indices()[0]) > *(t->get_indices()[1])) { killer = t->get_indices()[1]; killed = t->get_indices()[0]; }
	  else                                                { killer = t->get_indices()[0]; killed = t->get_indices()[1]; }
	  
	  *killed = *killer; 
	  t = Rtensors_.erase(t);
	} // End else
	++t;
      } // End else
      // In case of either of indices in kDelta is dummy .... 
      else if((t->get_name()==kDelta_name()) && (t->get_indices()[0]->get_isSummed() || t->get_indices()[1]->get_isSummed())){
        if(t->get_indices()[0]->get_char() != t->get_indices()[1]->get_char()){
          numConst_ = 0;
	  return;
	} // End if
	else{
          SharedIndex killed;
          SharedIndex killer;
          if(t->get_indices()[1]->get_isSummed()) { killed = t->get_indices()[1]; killer = t->get_indices()[0];}
          else                                    { killed = t->get_indices()[0]; killer = t->get_indices()[1];}

	  *killed = *killer; 
 	  t = Rtensors_.erase(t);
	} // End else
      } // End else if
      else ++t;

    } // End for

  }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::print_summedBody()
  {
    cout << ">> summedBody <<" << endl;
    auto i(summedIndices_.begin());
    int count(0);
    for(;i != summedIndices_.end();++i)
      cout << (boost::format("[%10d] ") % count++) << *i << " (" << **i << ") <" <<  ((*i)->get_isSummed() ? "D" : "N") << " : "<<  ((*i)->get_isExt() ? "E" : "I") << ">" << endl; 
  }

  // *********************************************************
  // 
  // *********************************************************
  bool SQcontract::operator<(const SQcontract &obj) const
  {
    bool retval;
    if      (numConst_  < obj.numConst_) retval = true;
    else if (numConst_  > obj.numConst_) retval = false;
    else if (numConst_ == obj.numConst_){
      SQcont<string> myConsts(Consts_);
      SQcont<string> o_Consts(obj.Consts_);
      if      (myConsts  < o_Consts) retval = true;
      else if (myConsts  > o_Consts) retval = false;
      else if (myConsts == o_Consts){
	if      (Ltensor_  < obj.Ltensor_) retval = true;
	else if (Ltensor_  > obj.Ltensor_) retval = false;
	else if (Ltensor_ == obj.Ltensor_){
	  retval = false;
	  SQtensors myRten(Rtensors_);
	  SQtensors o_Rten(obj.Rtensors_);
	  if      (myRten  < o_Rten) retval = true;
	  else if (myRten  > o_Rten) retval = false;
	  else if (myRten == o_Rten) retval = false;
	} // End else if
      } // End else if
    } // End else if
    return retval;
  }

  // *********************************************************
  // 
  // *********************************************************
  bool SQcontract::operator==(const SQcontract &obj) const
  {
    if(numConst_  != obj.numConst_) return false;
    if(!(Ltensor_ == obj.Ltensor_)) return false;
    SQcont<bool> myRinds(Rindices_);
    SQcont<bool> o_Rinds(obj.Rindices_);
    if(myRinds != o_Rinds) return false;
    auto avator1(*this);
    auto avator2(obj);
    return isFactorizable(avator1, avator2);
  }
    
  // *********************************************************
  // Return whether two SQcontract are factorizable or not
  // *********************************************************
  bool isFactorizable(SQcontract &a, SQcontract &b)
  {
    // Compare constants
    SQcont<string> Const_a(a.get_Consts());
    SQcont<string> Const_b(b.get_Consts());
    if(Const_a != Const_b) return false;

//*//     // Compare number of indices and names of (L)tensors
//*//     if(a.get_Ltensor().get_indices().size() != b.get_Ltensor().get_indices().size() || 
//*//        a.get_Ltensor().get_name() != b.get_Ltensor().get_name()) return false;

    // Compare names of (R)tensors and associated indices
    SQtensors ten_a(a.get_Rtensors());
    SQtensors ten_b(b.get_Rtensors());
    if(ten_a.size() != ten_b.size()) return false;
    SQcont<string> name_a;
    SQcont<string> name_b;
//*//     stable_sort(ten_a.begin(), ten_a.end());
//*//     stable_sort(ten_b.begin(), ten_b.end());
    for(size_t i = 0;i < ten_a.size();++i){
      if(!name_a.count(ten_a[i].get_name())) name_a <= ten_a[i].get_name();
      if(!name_b.count(ten_b[i].get_name())) name_b <= ten_b[i].get_name();
    } // End i
    if(name_a != name_b) return false;

    // Count all the dummies in each orbital group
    SharedIndices a_indices(a.get_summedBody());
    SharedIndices c_ptr;
    SharedIndices o_ptr;
    SharedIndices v_ptr;
    for(auto i = a_indices.begin();i != a_indices.end();++i){
      if     ((*i)->get_char()==core && (*i)->get_isSummed()) c_ptr <= *i; 
      else if((*i)->get_char()==act  && (*i)->get_isSummed()) o_ptr <= *i; 
      else if((*i)->get_char()==virt && (*i)->get_isSummed()) v_ptr <= *i;
    } // End i

    // Before comparing the tensors, masquerade the other one!
    b.masquerade();
    SQcont<SQtensor> b_tensors(b.get_Rtensors()); 
    Femto::for_each(b_tensors, boost::bind(&Femto::Core::SQtensor::sortIndices, _1));
    stable_sort(b_tensors.begin(), b_tensors.end());

    // Permute pairs of indices in each orbital group
    IIvector c_perms(makePermutations((int)c_ptr.size()));
    IIvector o_perms(makePermutations((int)o_ptr.size()));
    IIvector v_perms(makePermutations((int)v_ptr.size()));
    int c_max(c_perms.size() ? c_perms.size() : 1);
    int o_max(o_perms.size() ? o_perms.size() : 1);
    int v_max(v_perms.size() ? v_perms.size() : 1);
    for(int ic = 0;ic < c_max;++ic){
      for(int io = 0;io < o_max;++io){
        for(int iv = 0;iv < v_max;++iv){
          // Replace each name of index
          if(c_perms.size())
            //cout << c_perms.size() << endl;
            for(size_t i = 0;i < c_ptr.size();++i){
              stringstream num; num << c_perms[ic][i]+1;
              c_ptr[i]->put_index("c"+num.str());
            }
	  if(o_perms.size())
            for(size_t i = 0;i < o_ptr.size();++i){
              stringstream num; num << o_perms[io][i]+1;
              o_ptr[i]->put_index("o"+num.str());
            }
          if(v_perms.size())
            for(size_t i = 0;i < v_ptr.size();++i){
              stringstream num; num << v_perms[iv][i]+1;
              v_ptr[i]->put_index("v"+num.str()); 
            }
	  SQtensors a_tensors(a.get_Rtensors());
	  Femto::for_each(a_tensors, boost::bind(&Femto::Core::SQtensor::sortIndices, _1));
	  stable_sort(a_tensors.begin(), a_tensors.end());
	  if(a_tensors == b_tensors && a.get_Ltensor() == b.get_Ltensor()) return true;
        } // End iv
      } // End io
    } // End ic

    return false;
  }

  // *********************************************************
  // 
  // *********************************************************
  SQbinary::SQbinary()
  {}

  // *********************************************************
  // 
  // *********************************************************
  SQbinary::SQbinary(const double numConst, const SQcont<string> Consts, 
		     const SQtensor &Ltensor, const SQcont<SQtensor> &Rtensors)
  {

    numConst_ = numConst;
    Consts_   = Consts;
    Ltensor_  = Ltensor;
    Rtensors_ = Rtensors;
    Lindices_ = false;
    Rindices_;
    Rindices_.push_back(false);
    Rindices_.push_back(false);

    if(is_RDM(Ltensor_.get_name()) || is_sfGen(Ltensor_.get_name()) || Ltensor_.get_name() == kDelta_name())
      { cout << "SQcontract: Ltensor cannot be either of sfGen, RDM or kDelta" << endl; abort(); }
   
    int num_non_kdeltas = 0;
    for(auto t = Rtensors_.begin();t != Rtensors_.end();++t)
      if(t->get_name() != kDelta_name()) ++num_non_kdeltas;
    if(num_non_kdeltas >= 3){
      cout << "Number of non Kronecker's delta type tensors has to be less than 3." << endl;
      for(auto t = Rtensors_.begin();t != Rtensors_.end();++t) cout << *t;
      cout << endl;
      abort();
    } // End if
    for(auto t = Rtensors_.begin();t != Rtensors_.end();++t){
      if(is_sfGen(t->get_name())){
        cout << "Spin-free unitary group generator is detected: " << endl;
        abort();
      } // End if
    } // End t
    set_summedBody();
  }

  // *********************************************************
  // 
  // *********************************************************
  SQbinary::SQbinary(const SQbinary &obj)
  {
    numConst_ = obj.numConst_;
    Consts_   = obj.Consts_;
    Ltensor_  = obj.Ltensor_;
    Rtensors_ = obj.Rtensors_;
    Rindices_ = obj.Rindices_;
    set_summedBody();
  }

  // *********************************************************
  // 
  // *********************************************************
  SQbinary::SQbinary(const SQtensor &Ltensor, const SQterm &Rterm)
  {
    numConst_ = Rterm.get_numConst();
    Consts_   = Rterm.get_Consts();
    Ltensor_  = Ltensor;
    Rtensors_ = Rterm.get_tensors();
    Lindices_ = false;
    if(is_RDM(Ltensor_.get_name()) || is_sfGen(Ltensor_.get_name()) || Ltensor_.get_name() == kDelta_name())
      { cout << "SQbinary: Ltensor cannot be either of sfGen, RDM or kDelta" << endl; abort(); }
    for(auto t = Rtensors_.begin();t != Rtensors_.end();++t){
      Rindices_ <= false;
      if(is_sfGen(t->get_name())){
        cout << "Spin-free unitary group generator is detected: " << endl;
        abort();
      } // End if
    } // End t
    set_summedBody();    
  }

}} // Femto::

