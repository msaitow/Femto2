//
//  SQterm.cc
//  
//
//  Created by Masaaki Saitow on 12/06/30.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <boost/function.hpp>
#include <boost/format.hpp>
#include <boost/bind.hpp>
#include <Femto.hpp>
#include <SQterm.hpp>
#include <SQtensor.hpp>
#include <SQindex.hpp>

//#define _DEBUG_TERM

#define _SORT
//#define _NEW_FACT

#define _NEG_SUMMED  // Don't touch this
  
//#define _GEN_DEBUG   // Debug option for utilization of the generic index

using namespace std;

namespace Femto { namespace Core { 

  // *********************************************************
  // 
  // *********************************************************
  SQterm::SQterm()
  {}

  // *********************************************************
  // 
  // *********************************************************
  SQterm::~SQterm()
  {}

  // *********************************************************
  // 
  // *********************************************************
  SQterm::SQterm(const double numConst, const SQcont<string> Consts, 
                 const SQtensors Tensors, const bool isInCanonical)
    : numConst_(numConst),
      Consts_(Consts),
      isInCanonical_(isInCanonical)
  {
    int num_kDeltas(0);
    SQtensors Tnon_commutes;
    for(size_t I = 0;I < Tensors.size();++I){
      if(Tensors[I].isCommutable()){
        // Kronecker's delta with same indices aren't needed (2012/10/22)
	if(Tensors[I].get_name() == kDelta_name()){
	  if(find(Tensors_.begin(), Tensors_.end(), Tensors[I]) == Tensors_.end()) 
	    Tensors_ <= Tensors[I];
          else ++num_kDeltas;
	}
	else Tensors_ <= Tensors[I];
      }
      else   Tnon_commutes <= Tensors[I];
    }
    if(Tensors_.size()+Tnon_commutes.size()+num_kDeltas != Tensors.size()){
      cout << "Algorithmic Error .... " << endl;
      abort();
    }
    sort(Consts_.begin(), Consts_.end());
    sort(Tensors_.begin(), Tensors_.end());

    Tensors_.insert(Tensors_.end(), Tnon_commutes.begin(), Tnon_commutes.end());
#ifndef _NEG_SUMMED
    set_summedBody();
#endif
    // Count number of sfGen
    int count(0);
    for(size_t I = 0;I < Tensors_.size();++I){
      if(is_sfGen(Tensors_[I].get_name())) ++count;
    } // End I
    // Already normal ordered
    if(count == 0 || count == 1) set_isInCanonical(true);

#ifdef _DEBUG_TERM
    cout << "Content of the Tensors_ .... " << endl;               //*TEST* 
    for(size_t i = 0;i < Tensors_.size();++i) cout << Tensors_[i]; //*TEST*
    cout << endl << endl;                                          //*TEST*
#endif

  }

  // *********************************************************
  // 
  // *********************************************************
  SQterm::SQterm(const SQterm &obj)
    : isInCanonical_(obj.isInCanonical_),
      numConst_(obj.numConst_),
      Consts_(obj.Consts_),
      Tensors_(obj.Tensors_)
  {
    //#ifndef _NEG_SUMMED 
    set_summedBody(); 
    //#endif
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQterm::set_summedBody()
  {

    SharedIndices allIndices;
    // Copy all the indices found in Tensors_ into allIndices (is summedBody_ is already set)
    if(summedIndices_.size()){
      for(auto t = Tensors_.begin();t != Tensors_.end();++t){
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
      
      // Link indices in Tensors_ to allIndices
      for(auto t = Tensors_.begin();t != Tensors_.end();++t){
	for(size_t inum = 0;inum < t->get_indices().size();++inum){
	  for(auto j = allIndices.begin();j != allIndices.end();++j)
	    if(*t->get_indices()[inum] == **j) t->put_indices(inum, *j);
	} // End inum
      } // End t
      summedIndices_.clear();
      summedIndices_.reserve(allIndices.size());
    } // End if
    // if not just copy the indices to the summedBody_
    else{
      for(auto t = Tensors_.begin();t != Tensors_.end();++t){
	for(size_t inum = 0;inum < t->get_indices().size();++inum){
	  bool is_OK(true);
	  for(auto I = allIndices.begin();I != allIndices.end();++I){
	    if(**I==*t->get_indices()[inum]) is_OK = false;
	  } // End I
	  if(is_OK) allIndices <= t->get_indices()[inum];
	} // End inum 
      } // End t
    } // End else

    // Copy the content of allIndices to summedIndices_
    for(auto i = allIndices.begin();i != allIndices.end();++i)
      summedIndices_ <= SharedIndex( new SQindex((*i)->get_index(),
						 (*i)->get_char(),
						 (*i)->get_isSummed(),
						 (*i)->get_isExt()));

    // Then re-link to summedIndices_
    for(auto t = Tensors_.begin();t != Tensors_.end();++t){
      for(size_t inum = 0;inum < t->get_indices().size();++inum){
	for(auto j = summedIndices_.begin();j != summedIndices_.end();++j)
	  if(*t->get_indices()[inum] == **j) t->put_indices(inum, *j);
      } // End inum
    } // End t

    ///////////////////////// Outdated implementation ////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    // Copy all the indices found in Tensors_ into allIndices
    // SQcont<SharedIndex> allIndices;
    // for(auto t = Tensors_.begin();t != Tensors_.end();++t){
    //   for(size_t inum = 0;inum < t->get_indices().size();++inum){
    // 	bool is_OK(true);
    //     for(auto I = allIndices.begin();I != allIndices.end();++I){
    //       if(**I==*t->get_indices()[inum]) is_OK = false;
    // 	} // End I
    //     if(is_OK){
    //       allIndices <= SharedIndex( new SQindex(t->get_indices()[inum]->get_index(),
    //                                              t->get_indices()[inum]->get_char(),
    //                                              t->get_indices()[inum]->get_isSummed(),
    // 						 t->get_indices()[inum]->get_isExt()));
    // 	} // End if
    //   } // End inum 
    // } // End t
    // 
    // // Link indices in Tensors_ to allIndices
    // for(auto t = Tensors_.begin();t != Tensors_.end();++t){
    //   for(size_t inum = 0;inum < t->get_indices().size();++inum){
    // 	for(auto j = allIndices.begin();j != allIndices.end();++j)
    // 	  if(*t->get_indices()[inum] == **j) t->put_indices(inum, *j);
    //   } // End inum
    // } // End t
    // 
    // summedIndices_.clear();
    // summedIndices_.reserve(allIndices.size());
    // 
    // // Copy the content of allIndices to summedIndices_
    // for(auto i = allIndices.begin();i != allIndices.end();++i)
    //   summedIndices_ <= SharedIndex( new SQindex((*i)->get_index(),
    // 						 (*i)->get_char(),
    // 						 (*i)->get_isSummed(),
    // 						 (*i)->get_isExt()));
    // 
    // // Then re-link to summedIndices_
    // for(auto t = Tensors_.begin();t != Tensors_.end();++t){
    //   for(size_t inum = 0;inum < t->get_indices().size();++inum){
    // 	for(auto j = summedIndices_.begin();j != summedIndices_.end();++j)
    // 	  if(*t->get_indices()[inum] == **j) t->put_indices(inum, *j);
    //   } // End inum
    // } // End t
    //////////////////////////////////////////////////////////////////////////////////////

    //////////////////// From the previous version ////////////////////////////////
    // ////////////////////////////////////////////////////////////////////////////
    // // Firstly, assign all the &summedIndices_[:] to &tempIndices[:]
    // vector<SQindex> tempIndices(summedIndices_);
    // for(size_t i = 0;i < tempIndices.size();++i){
    //   for(size_t I = 0;I < Tensors_.size();++I){
    //     vector<SQindex*> temp1(Tensors_[I].get_indices());
    //     for(size_t j = 0;j < temp1.size();++j){
    //       if(*temp1[j]==tempIndices[i]) Tensors_[I].put_indices(j, &tempIndices[i]);
    // 	} // End j
    //   } // End I
    // } // End i
    // summedIndices_.clear();
    // 
    // // Search for the dummy indices .... 
    // for(size_t I = 0;I < Tensors_.size();++I){
    //   vector<SQindex*> temp1(Tensors_[I].get_indices());
    //   vector<SQindex> indices;
    //   for(size_t i = 0;i < temp1.size();++i) indices.push_back(*temp1[i]);
    //   for(size_t j = 0;j < indices.size();++j){
    //     if(find(summedIndices_.begin(), summedIndices_.end(), indices[j])==summedIndices_.end()) {
    //       summedIndices_.push_back(*temp1[j]);
    //     } // End if
    //   } // End j
    // } // End I
    // 
    // // Replace all the dummy indices, which are shared among all the terms, with those
    // // copied in the summedIndices
    // for(size_t i = 0;i < summedIndices_.size();++i){
    //   SQindex temp_i(summedIndices_[i]);
    //   for(size_t I = 0;I < Tensors_.size();++I){
    //     vector<SQindex*> temp1(Tensors_[I].get_indices());
    //     vector<SQindex> indices;        
    //     for(size_t j = 0;j < temp1.size();++j) indices.push_back(*temp1[j]);
    //     for(size_t j = 0;j < indices.size();++j){
    //       if(indices[j] == temp_i) Tensors_[I].put_indices(j, &summedIndices_[i]);
    // 	}
    //   } // End I
    // } // End i
    //////////////////////////////////////////////////////////////////////////////////////////

  }


    // *********************************************************
    // 
    // *********************************************************
    void SQterm::masquerade()
    {
      // Set the dummy names for all the dummy indices
      int Ccount(0);
      int Ocount(0);
      int Vcount(0);
      int Gcount(0);
      for(size_t i = 0;i < summedIndices_.size();++i){
        // In case of core
        if(summedIndices_[i]->get_char() == core && summedIndices_[i]->get_isSummed()){
          ++Ccount;
          ostringstream stm;
          stm << Ccount;
          summedIndices_[i]->put_index("c" + stm.str());
        } // End if
	// In case of active
        else if(summedIndices_[i]->get_char() == act && summedIndices_[i]->get_isSummed()){
	  ++Ocount;
	  ostringstream stm;
	  stm << Ocount;
	  summedIndices_[i]->put_index("o" + stm.str());
	} // End if
	// In case of virtual
	else if(summedIndices_[i]->get_char() == virt && summedIndices_[i]->get_isSummed()){
	  ++Vcount;
	  ostringstream stm;
	  stm << Vcount;
	  summedIndices_[i]->put_index("v" + stm.str());
	} // End if    
	// In case of generic
	else if(summedIndices_[i]->get_char() == gen && summedIndices_[i]->get_isSummed()){
	  ++Gcount;
	  ostringstream stm;
	  stm << Gcount;
	  summedIndices_[i]->put_index("g" + stm.str());
	} // End if    
      } // End i
      for(auto t = Tensors_.begin();t != Tensors_.end();++t){
        t->sortIndices();
      } // End t    
      
    }
  
//*NOT_YET* 
//*NOT_YET* //*INCOMPLETE*   // *********************************************************
//*NOT_YET* //*INCOMPLETE*   // 
//*NOT_YET* //*INCOMPLETE*   // *********************************************************
//*NOT_YET* //*INCOMPLETE*   void SQterm::erase_duplication()
//*NOT_YET* //*INCOMPLETE*   {
//*NOT_YET* //*INCOMPLETE*     vector<int> count;   count.resize((int)summedIndices_.size());
//*NOT_YET* //*INCOMPLETE*     vector<int> location;location.resize((int)summedIndices_.size());
//*NOT_YET* //*INCOMPLETE*     for(size_t i = 0;i < summedIndices_.size()-1;++i){
//*NOT_YET* //*INCOMPLETE*       for(size_t j = i+1;j < summedIndices_.size();++j){
//*NOT_YET* //*INCOMPLETE*         if(summedBody_[i]==summedBody_[j]) { ++count[i]; location[i] = j; }
//*NOT_YET* //*INCOMPLETE*       } // End j    
//*NOT_YET* //*INCOMPLETE*     } // End i
//*NOT_YET* //*INCOMPLETE*     for(size_t i = 0;i < summedIndices.size();++i){
//*NOT_YET* //*INCOMPLETE*       if(count[i]){
//*NOT_YET* //*INCOMPLETE*         for(vector<SQtensor>::iterator t = Tensors_.begin();t != Tensors.end();++t){
//*NOT_YET* //*INCOMPLETE* 	  if(*(t->get_indices()) == summedBody_[i]) t->put_indices(i, &summedBody[i]);
//*NOT_YET* //*INCOMPLETE* 	} // End t
//*NOT_YET* //*INCOMPLETE*         summedIndices
//*NOT_YET* //*INCOMPLETE*       } // End if
//*NOT_YET* //*INCOMPLETE*     } // End for  
//*NOT_YET* //*INCOMPLETE*   }

    // *********************************************************
    // 
    // *********************************************************
    SQtensors SQterm::get_tensors() const
    { return Tensors_; }

    // *********************************************************
    // 
    // *********************************************************
    SQcont<SharedTensor> SQterm::get_tensors_ptr()
    {
      SQcont<SharedTensor> retval;
      for(size_t i = 0;i < Tensors_.size();++i) retval <= &(Tensors_[i]); 
      return retval; 
    }

    // *********************************************************
    // 
    // *********************************************************
    double SQterm::get_numConst() const
    { return numConst_; }
 
    // *********************************************************
    // 
    // *********************************************************
    SQcont<string> SQterm::get_Consts() const
    { return Consts_; }
 
    // *********************************************************
    // 
    // *********************************************************
    bool SQterm::get_isInCanonical() const
    { return isInCanonical_; }
 
    // *********************************************************
    // 
    // *********************************************************
    void SQterm::set_isInCanonical(const bool isInCanonical)
    { isInCanonical_ = isInCanonical; }

    // *********************************************************
    // 
    // *********************************************************
    void SQterm::set_tensors(const SQcont<SQtensor> &Tensors)
    {
      SQtensors Tcommutes, Tnon_commutes;
      Tcommutes.reserve(Tensors.size());
      Tnon_commutes.reserve(Tensors.size());
      for(size_t I = 0;I < Tensors.size();++I){
	if(Tensors[I].isCommutable()) Tcommutes <= Tensors[I];
	else                          Tnon_commutes <= Tensors[I];
      }
      if(Tcommutes.size()+Tnon_commutes.size() != Tensors.size()){
	cout << "Algorithmic Error .... " << endl;
	abort();
      }
      sort(Tcommutes.begin(), Tcommutes.end());
      
      Tensors_.clear();
      Tensors_.reserve(Tensors.size());
      
      // Modified 2012/10/23
      Tensors_.insert(Tensors_.begin(), Tcommutes.begin(),     Tcommutes.end());
      Tensors_.insert(Tensors_.end(),   Tnon_commutes.begin(), Tnon_commutes.end());
#ifndef _NEG_SUMMED            
      set_summedBody();
#endif      
    }
 
    // *********************************************************
    // 
    // *********************************************************
    SQterm SQterm::operator=(const SQterm &obj){
      isInCanonical_ = obj.isInCanonical_;
      numConst_      = obj.numConst_;     
      Consts_        = obj.Consts_;       
      Tensors_       = obj.Tensors_;      
      //#ifndef _NEG_SUMMED     
      set_summedBody();
      //#endif      
      return (*this);
    }

    // *********************************************************
    // 
    // *********************************************************
    bool SQterm::operator==(const SQterm &obj)
    {
      if(numConst_ != obj.numConst_) return false;
      if(Consts_   != obj.Consts_)   return false;
      if(Tensors_  != obj.Tensors_)  return false;
      
      return true;
    }
    
    // *********************************************************
    // 
    // *********************************************************
    void SQterm::set_numConst(const double num)
    { numConst_ = num; }
    
    // *********************************************************
    // 
    // *********************************************************
    void SQterm::set_Consts(const SQcont<string> consts)
    { Consts_ = consts; }
 
    // *********************************************************
    // 
    // *********************************************************
    ostream& operator <<(std::ostream &os, const SQterm &t)
    {
      os << boost::format("(%14.8f)") % t.get_numConst();
      if(t.get_Consts().size())
	os << " " <<= t.get_Consts();
      if(t.get_tensors().size())
	os << " " <<= t.get_tensors();       

      return os;
    }

    // *********************************************************
    // 
    // *********************************************************
    SQcont<SharedIndex> SQterm::get_summedBody()
    {
      SharedIndices retval;
      auto ind(summedIndices_.begin());
      for(;ind != summedIndices_.end();++ind) retval <= *ind;
      return retval;
    }

    // *********************************************************
    // If the term is in canonical order, transform sfGen to RDM 
    // *********************************************************
    void SQterm::transform2RDM()
    {
      //set_summedBody(); //*TEST* 
      auto ten(Tensors_.begin());
      if(!get_isInCanonical()) return;
      bool set_flag(false); 
      for(;ten != Tensors_.end();){
	if(is_sfGen(ten->get_name())){
	  RDM temp(ten->get_indices());
	  Tensors_.erase(ten);
	  Tensors_ <= temp;
	  set_flag = true;        
	} // End if
	else ++ten;
      } // End ten
      if(set_flag) {
	SQtensors commT;
	SQtensors ncommT;
	for(size_t i = 0;i < Tensors_.size();++i)
	  if(Tensors_[i].isCommutable())  commT <= Tensors_[i];
	  else                           ncommT <= Tensors_[i];
	
	sort(commT.begin(), commT.end());
	commT.insert(commT.end(), ncommT.begin(), ncommT.end());
	Tensors_ = commT;
	//set_summedBody(); //*TEST*
      }
    }
 
    // *********************************************************
    // Contract Kronecker deltas
    // *********************************************************
    void SQterm::contractkDeltas()
    {
#ifdef _NEG_SUMMED
      set_summedBody();
#endif
      for(auto t = Tensors_.begin(); t != Tensors_.end();){
	// In case of delta_{p}^{p} ....
	if((t->get_name()==kDelta_name()) && (*(t->get_indices()[0]) == *(t->get_indices()[1])))
	  t = Tensors_.erase(t);
	// In case of both indices of kDelta are dummy () ....
	else if((t->get_name()==kDelta_name()) && !(t->get_indices()[0]->get_isSummed()) && !(t->get_indices()[1]->get_isSummed())){
	  if(t->get_indices()[0]->get_char() != t->get_indices()[1]->get_char() /*&& (t->get_indices()[0]->get_char() != gen || t->get_indices()[1]->get_char() != gen)*/){
	    numConst_ = 0; // This term shoud be zero
	    return;
	  } // End if
	  ++t;
	} // End else if
	// In case of either of indices in kDelta is dummy ....
	// If both indices are different each other and not of generic indices 
	else if((t->get_name()==kDelta_name()) && (t->get_indices()[0]->get_isSummed() || t->get_indices()[1]->get_isSummed())){
	  if(t->get_indices()[0]->get_char() != t->get_indices()[1]->get_char() && t->get_indices()[0]->get_char()!=gen && t->get_indices()[1]->get_char()!=gen){
	    numConst_ = 0;
	    return;
	  } // End if
	  else{
//*SLEEP* 	    cout<<"LLLL" << endl;
//*SLEEP* 	    print_summedBody(); //*TEST*
//*SLEEP* 	    cout<<"MMMM" << endl;
	    SharedIndex killed;
	    SharedIndex killer;
	    if(t->get_indices()[0]->get_isSummed() && t->get_indices()[1]->get_isSummed())
	      if(t->get_indices()[0]->get_char() == gen)
		{ killed =  t->get_indices()[0]; killer = t->get_indices()[0]; }
	      else if(t->get_indices()[1]->get_char() == gen)
		{ killed =  t->get_indices()[1]; killer = t->get_indices()[0]; }
	      else                       
		{ killed = t->get_indices()[0]; killer = t->get_indices()[1]; }
	    else if(t->get_indices()[1]->get_isSummed()) { killed = t->get_indices()[1]; killer = t->get_indices()[0]; }
	    else                                         { killed = t->get_indices()[0]; killer = t->get_indices()[1]; }
//*SLEEP* 	    cout << " * Killed " << *killed << endl;
//*SLEEP* 	    cout << " * Killer " << *killer << endl;
	    *killed = *killer; 
//*SLEEP* 	    cout<<"NNNN" << endl;
//*SLEEP* 	    print_summedBody(); //*TEST*
//*SLEEP* 	    cout<<"OOOO" << endl;
	    t = Tensors_.erase(t);
	  } // End else
	} // End else if
	else ++t;
	
      } // End for
      //set_summedBody(); //*SLOW* ???
      // This should be improved like  
    }
    
    // *********************************************************
    // Decompose sfGen, or RDM (only if isInCanonical = true)
    // *********************************************************
    void SQterm::decomposeRDMVirt()
    {
      if(!isInCanonical_) return;
      auto t(Tensors_.begin());
      for(;t != Tensors_.end();++t){
	if(is_sfGen(t->get_name()) || is_RDM(t->get_name())){
	  SharedIndices indices(t->get_indices());
	  auto i(indices.begin());
	  for(;i != indices.end();++i){
	    // In case of core .... 
	    if((*i)->get_char() == core) {
	      /* Do for core indices */
	    } // End if
	    // In case of virtual ....
	    else if((*i)->get_char() == virt) {
	      numConst_ = 0;
	      return;
	    } // End else if
	  } // End for
	} // End if
      } // End for
    }
 
    // *********************************************************
    // 
    // *********************************************************
    void SQterm::print_summedBody()
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
    string SQterm::convert2LaTeX() const
    {
      string term("");
      if     (numConst_ == 0.0) return "";
      
      ostringstream stm;
      stm << fabs(numConst_);
      term += (numConst_ > 0.0 ? "+ " : "- ");
      if(fabs(numConst_) != 1.0) term += stm.str() + " ";
      for(size_t i = 0;i < Consts_.size();++i)
	if(Consts_[i] != "") term += Consts_[i] + " ";
      
      for(size_t i = 0;i < Tensors_.size();++i)
	term += Tensors_[i].convert2LaTeX();
      return term;
    }
 
    // *********************************************************
    // Return whether two SQterms are additive or not
    // *********************************************************
    // *NOTE* :: If there are common tensor type in terms, the current algorithm 
    //           may not work perfectly (2012/10/22)
    bool isAdditive(SQterm &a, SQterm &b)
    {
      // Compare constants
      SQcont<string> Const_a(a.get_Consts());
      SQcont<string> Const_b(b.get_Consts());
      if(Const_a != Const_b) return false;
      
#ifndef _NEW_FACT
      // Modified 2012/10/24
      return isFactorizable(a, b);
      //return isFactorizable_new(a, b);
#else
      return isFactorizable2(a, b);
#endif

//* --- *     // Compare names of tensors and associated indices
//* --- *     vector<SQtensor> ten_a(a.get_tensors());
//* --- *     vector<SQtensor> ten_b(b.get_tensors());
//* --- *     if(ten_a.size() != ten_b.size()) return false;
//* --- *     for(size_t i = 0;i < ten_a.size();++i){
//* --- *       if(ten_a[i].get_name()           != ten_b[i].get_name())           
//* --- *         return false;
//* --- *       if(ten_a[i].get_indices().size() != ten_b[i].get_indices().size())
//* --- *         return false;
//* --- *     } // End i
//* --- * 
//* --- *     // Count all the dummies in each orbital group
//* --- *     vector<SQindex*> a_indices(a.get_summedBody());
//* --- *     vector<SQindex*> c_ptr;
//* --- *     vector<SQindex*> o_ptr;
//* --- *     vector<SQindex*> v_ptr;
//* --- *     for(vector<SQindex*>::iterator i = a_indices.begin();i != a_indices.end();++i){
//* --- *       if     ((*i)->get_char()==(char_state)0 && (*i)->get_isSummed()) c_ptr.push_back(*i); 
//* --- *       else if((*i)->get_char()==(char_state)1 && (*i)->get_isSummed()) o_ptr.push_back(*i); 
//* --- *       else if((*i)->get_char()==(char_state)2 && (*i)->get_isSummed()) v_ptr.push_back(*i);
//* --- *     } // End i
//* --- * 
//* --- *     // Before comparing the tensors, masquerade the other one!
//* --- *     b.masquerade();
//* --- * 
//* --- *     // Permute pairs of indices in each orbital group
//* --- *     IIvector c_perms = makePermutations((int)c_ptr.size());
//* --- *     IIvector o_perms = makePermutations((int)o_ptr.size());
//* --- *     IIvector v_perms = makePermutations((int)v_ptr.size());
//* --- *     int c_max = (c_perms.size() ? c_perms.size() : 1);
//* --- *     int o_max = (o_perms.size() ? o_perms.size() : 1);
//* --- *     int v_max = (v_perms.size() ? v_perms.size() : 1);
//* --- *     for(int ic = 0;ic < c_max;++ic){
//* --- *       for(int io = 0;io < o_max;++io){
//* --- *         for(int iv = 0;iv < v_max;++iv){
//* --- *           // Replace each name of index
//* --- *           if(c_perms.size())
//* --- *             //cout << c_perms.size() << endl;
//* --- *             for(size_t i = 0;i < c_ptr.size();++i){
//* --- *               stringstream num; num << c_perms[ic][i]+1;
//* --- *               c_ptr[i]->put_index("c"+num.str());
//* --- *             }
//* --- * 	  if(o_perms.size())
//* --- *             for(size_t i = 0;i < o_ptr.size();++i){
//* --- *               stringstream num; num << o_perms[io][i]+1;
//* --- *               o_ptr[i]->put_index("o"+num.str());
//* --- *             }
//* --- *           if(v_perms.size())
//* --- *             for(size_t i = 0;i < v_ptr.size();++i){
//* --- *               stringstream num; num << v_perms[iv][i]+1;
//* --- *               v_ptr[i]->put_index("v"+num.str()); 
//* --- *             }
//* --- * 	    if(a.get_tensors() == b.get_tensors()) return true;
//* --- *         } // End iv
//* --- *       } // End io
//* --- *     } // End ic
//* --- * 
//* --- *     return false;
    }


    // *Created (2012/10/24)*
    // *********************************************************
    // Return whether two SQterms are factorizable or not
    // *********************************************************
    // *NOTE* :: If there are common tensor type in terms, the current algorithm 
    //           may not work perfectly (2012/10/22)
    bool isFactorizable(SQterm &a, SQterm &b)
    {
      
      // Compare names of tensors and associated indices
      SQtensors ten_a(a.get_tensors());
      SQtensors ten_b(b.get_tensors());
      if(ten_a.size() != ten_b.size()) return false;
      for(size_t i = 0;i < ten_a.size();++i){
	if(ten_a[i].get_name()           != ten_b[i].get_name())           
	  return false;
	if(ten_a[i].get_indices().size() != ten_b[i].get_indices().size())
	  return false;
      } // End i
      
      //*UNCOMPLETE*     // Consider the which tensor holds the non-dummy indices
      //*UNCOMPLETE*     map<SQindex*, vector<string> > a_ind2ten, b_ind2ten;
      //*UNCOMPLETE*     for(size_t a = 0;a < ten_a.size();++a){
      //*UNCOMPLETE*       for(size_t i = 0;i < ten_a[a].get_indices().size();++i){
      //*UNCOMPLETE*         if()
      //*UNCOMPLETE*       } // End i 
      //*UNCOMPLETE*     } // End a
      
      // Count all the dummies in each orbital group
      SharedIndices a_indices(a.get_summedBody());
      SharedIndices c_ptr;
      SharedIndices o_ptr;
      SharedIndices v_ptr;
      SharedIndices g_ptr;
      for(auto i = a_indices.begin();i != a_indices.end();++i){
	if     ((*i)->get_char()==core && (*i)->get_isSummed()) c_ptr <= *i; 
	else if((*i)->get_char()==act  && (*i)->get_isSummed()) o_ptr <= *i; 
	else if((*i)->get_char()==virt && (*i)->get_isSummed()) v_ptr <= *i;
	else if((*i)->get_char()==gen  && (*i)->get_isSummed()) g_ptr <= *i;
      } // End i
      
      // Before comparing the tensors, masquerade the other one!
      b.masquerade();
      SQcont<SQtensor> b_tensors(b.get_tensors());
      Femto::for_each(b_tensors, boost::bind(&Femto::Core::SQtensor::sortIndices, _1));  
      //OLD//for(auto t = b_tensors.begin();t != b_tensors.end();++t) t->sortIndices();
      sort(b_tensors.begin(), b_tensors.end()); //*TEST*

      // Permute pairs of indices in each orbital group
      IIvector c_perms(makePermutations((int)c_ptr.size()));
      IIvector o_perms(makePermutations((int)o_ptr.size()));
      IIvector v_perms(makePermutations((int)v_ptr.size()));
      IIvector g_perms(makePermutations((int)g_ptr.size()));
      int c_max(c_perms.size() ? c_perms.size() : 1);
      int o_max(o_perms.size() ? o_perms.size() : 1);
      int v_max(v_perms.size() ? v_perms.size() : 1);
      int g_max(g_perms.size() ? g_perms.size() : 1);
      for(int ic = 0;ic < c_max;++ic){
	for(int io = 0;io < o_max;++io){
	  for(int iv = 0;iv < v_max;++iv){
	    for(int ig = 0;ig < g_max;++ig){
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
	    if(g_perms.size())
	      for(size_t i = 0;i < g_ptr.size();++i){
		stringstream num; num << g_perms[ig][i]+1;
		g_ptr[i]->put_index("g"+num.str()); 
	      }
#ifdef _SORT
	    SQcont<SQtensor> a_tensors(a.get_tensors());
	    //OLD//for(auto t = a_tensors.begin();t != a_tensors.end();++t) t->sortIndices();
	    Femto::for_each(a_tensors, boost::bind(&Femto::Core::SQtensor::sortIndices, _1));
 	    sort(a_tensors.begin(), a_tensors.end()); //*TEST* 

	    bool all_ok(true);
	    for(size_t num = 0;num < a_tensors.size();++num)
	      if(!are_sameforms(a_tensors[num], b_tensors[num])) all_ok = false;
 	    if(all_ok) return true;

 	    //if(a_tensors == b_tensors) return true; //*TEST
#else
 	    if(a.get_tensors() == b.get_tensors()) return true;
#endif
	    } // End ig
	  } // End iv
	} // End io
      } // End ic
      
      return false;
    }
    
 
    // *Created (2012/11/14)*
    // *********************************************************
    // Return whether two SQterms are factorizable or not
    // *********************************************************
    // *NOTE* :: This will work perfectly, even if there are common tensors
    //           But it seems to be somewhat heavier
    bool isFactorizable2(SQterm &a, SQterm &b)
    {
      // Compare names of tensors and associated indices
      SQtensors ten_a(a.get_tensors());
      SQtensors ten_b(b.get_tensors());
      if(ten_a.size() != ten_b.size()) return false;
      
      //Divide the tensors by the names   
      typedef map<string, SQcont<SharedTensor> > ten_ptr; // name of tensor and the positions 
      ten_ptr a_ten, b_ten;
      for(size_t i = 0;i < a.get_tensors().size();++i) 
	a_ten[a.get_tensors()[i].get_name()].push_back(a.get_tensors_ptr()[i]);
      for(size_t i = 0;i < b.get_tensors().size();++i) 
	b_ten[b.get_tensors()[i].get_name()].push_back(b.get_tensors_ptr()[i]);
      for(auto i = a_ten.begin();i != a_ten.end();++i){
	string a_key(i->first);
	if(b_ten.find(a_key) == b_ten.end())        return false;
	if(b_ten[a_key].size() != i->second.size()) return false;
      } // End i
      
      int count(0);
      vector<IIvector> i_perms(a_ten.size());
      for(auto i = a_ten.begin();i != a_ten.end();++i,++count) i_perms[count] = makePermutations((int)i->second.size());
      
      // Count all the dummies in each orbital group
      SharedIndices a_indices(a.get_summedBody());
      SharedIndices c_ptr;
      SharedIndices o_ptr;
      SharedIndices v_ptr;
      SharedIndices g_ptr;
      for(auto i = a_indices.begin();i != a_indices.end();++i){
	if     ((*i)->get_char()==core && (*i)->get_isSummed()) c_ptr.push_back(*i); 
	else if((*i)->get_char()==act  && (*i)->get_isSummed()) o_ptr.push_back(*i); 
	else if((*i)->get_char()==virt && (*i)->get_isSummed()) v_ptr.push_back(*i);
	else if((*i)->get_char()==gen  && (*i)->get_isSummed()) g_ptr.push_back(*i);
      } // End i
      
      // Before comparing the tensors, masquerade the other one!
      b.masquerade();    
      
      // Permute pairs of indices in each orbital group
      IIvector c_perms(makePermutations((int)c_ptr.size()));
      IIvector o_perms(makePermutations((int)o_ptr.size()));
      IIvector v_perms(makePermutations((int)v_ptr.size()));
      IIvector g_perms(makePermutations((int)g_ptr.size()));
      int c_max(c_perms.size() ? c_perms.size() : 1);
      int o_max(o_perms.size() ? o_perms.size() : 1);
      int v_max(v_perms.size() ? v_perms.size() : 1);
      int g_max(g_perms.size() ? g_perms.size() : 1);
      for(int ic = 0;ic < c_max;++ic){
	for(int io = 0;io < o_max;++io){
	  for(int iv = 0;iv < v_max;++iv){
	    for(int ig = 0;ig < g_max;++ig){
	    // Replace each name of index
	    if(c_perms.size())
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
	    if(g_perms.size())
	      for(size_t i = 0;i < g_ptr.size();++i){
		stringstream num; num << g_perms[ig][i]+1;
		g_ptr[i]->put_index("g"+num.str()); 
	      }
	    // Compare all the members of the corresponding subgroups
	    // in every possible orders
	    {  
	      bool all_ok(true);
	      int count(0);
	      for(auto i = a_ten.begin();i != a_ten.end();++i,++count){
		bool ok;
		//IIvector i_perms(makePermutations((int)i->second.size()));
		for(int num_p = 0;num_p < i_perms[count].size();++num_p){
		  ok = true;
		  for(size_t num = 0;num < i->second.size();++num){
		    if(!(*(b_ten[i->first].at(num)) == *(i->second.at(i_perms[count][num_p][num])))) ok = false;
		  } // End num
		  if(ok) break;
		} // End p
		if(!ok) { all_ok = false; break; }
	      } // End i
	      if(all_ok) return true;
	    } // End scope
	    } // End ig
	  } // End iv
	} // End io
      } // End ic
      
      return false;    
    }
     

    // *Created (2012/10/24)*
    // *********************************************************
    // Return whether two SQterms are factorizable or not
    // *********************************************************
    // *NOTE* :: If there are common tensor type in terms, the current algorithm 
    //           may not work perfectly (2012/10/22)
    bool isFactorizable_new(SQterm &a, SQterm &b)
    {
      
      // Compare names of tensors and associated indices
      SQtensors ten_a(a.get_tensors());
      SQtensors ten_b(b.get_tensors());
      if(ten_a.size() != ten_b.size()) return false;
      for(size_t i = 0;i < ten_a.size();++i){
	if(ten_a[i].get_name()           != ten_b[i].get_name())           
	  return false;
	if(ten_a[i].get_indices().size() != ten_b[i].get_indices().size())
	  return false;
      } // End i
      
      //*UNCOMPLETE*     // Consider the which tensor holds the non-dummy indices
      //*UNCOMPLETE*     map<SQindex*, vector<string> > a_ind2ten, b_ind2ten;
      //*UNCOMPLETE*     for(size_t a = 0;a < ten_a.size();++a){
      //*UNCOMPLETE*       for(size_t i = 0;i < ten_a[a].get_indices().size();++i){
      //*UNCOMPLETE*         if()
      //*UNCOMPLETE*       } // End i 
      //*UNCOMPLETE*     } // End a
      
      // Count all the dummies in each orbital group
      SharedIndices a_indices(a.get_summedBody());
      SharedIndices c_ptr;
      SharedIndices o_ptr;
      SharedIndices v_ptr;
      SharedIndices g_ptr;
      for(auto i = a_indices.begin();i != a_indices.end();++i){
	if     ((*i)->get_char()==core && (*i)->get_isSummed()) c_ptr.push_back(*i); 
	else if((*i)->get_char()==act  && (*i)->get_isSummed()) o_ptr.push_back(*i); 
	else if((*i)->get_char()==virt && (*i)->get_isSummed()) v_ptr.push_back(*i);
	else if((*i)->get_char()==gen  && (*i)->get_isSummed()) g_ptr.push_back(*i);
      } // End i
      
      // Before comparing the tensors, masquerade the other one!
      b.masquerade();
      
      // Permute pairs of indices in each orbital group
      IIvector c_perms(makePermutations((int)c_ptr.size()));
      IIvector o_perms(makePermutations((int)o_ptr.size()));
      IIvector v_perms(makePermutations((int)v_ptr.size()));
      IIvector g_perms(makePermutations((int)g_ptr.size()));
      int c_max(c_perms.size() ? c_perms.size() : 1);
      int o_max(o_perms.size() ? o_perms.size() : 1);
      int v_max(v_perms.size() ? v_perms.size() : 1);
      int g_max(g_perms.size() ? g_perms.size() : 1);
      for(int ic = 0;ic < c_max;++ic){
	for(int io = 0;io < o_max;++io){
	  for(int iv = 0;iv < v_max;++iv){
	    for(int ig = 0;ig < g_max;++ig){
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
	    if(g_perms.size())
	      for(size_t i = 0;i < g_ptr.size();++i){
		stringstream num; num << g_perms[ig][i]+1;
		g_ptr[i]->put_index("g"+num.str()); 
	      }
#ifdef _SORT
 	    SQcont<SQtensor> a_tensors(a.get_tensors()); sort(a_tensors.begin(), a_tensors.end()); //*TEST* 
 	    //vector<SQtensor> b_tensors(b.get_tensors()); sort(b_tensors.begin(), b_tensors.end()); //*TEST*
 	    if(a_tensors == b.get_tensors()) return true; //*TEST
#else
	    bool all_ok(true);
	    for(size_t num = 0;num < a.get_tensors().size();++num) {
	      a.get_tensors_ptr()[num]->sortIndices();
	      b.get_tensors_ptr()[num]->sortIndices();
	      if(!are_sameforms(a.get_tensors()[num], b.get_tensors()[num])) all_ok = false;
	    } // End if  
 	    if(all_ok) return true;
#endif
	    } // End ig
	  } // End iv
	} // End io
      } // End ic
      
      return false;
    }

 
    // *********************************************************
    // Delete all the terms in inTerms with zero factors
    // *********************************************************
    void screenTerms(const SQcont<SQterm> &inTerms, SQcont<SQterm> *outTerms, const double crit)
    {
      outTerms->reserve(inTerms.size());
      for(size_t i = 0;i < inTerms.size();++i){
	if(fabs(inTerms[i].get_numConst()) > crit) outTerms->push_back(inTerms[i]);
      } // End i
    }
 
    // *********************************************************
    // Delete all the terms in inTerms with zero factors
    // *********************************************************
    SQterms screenTerms(const SQterms &inTerms, const double crit)
    {
      SQterms retval;
      for(size_t i = 0;i < inTerms.size();++i){
	if(fabs(inTerms[i].get_numConst()) > crit) retval.push_back(inTerms[i]);
      } // End for 
      return retval;
    }
 

    // *********************************************************
    // Combine terms with isAdditive()==true
    // *********************************************************
    SQterms combineTerms(SQterms &inTerms)
    {
      SQterms retval;
      int count(0);
      for(auto t1 = inTerms.begin();t1 != inTerms.end();++t1, ++count){
	//retval.push_back(*t1);
	//vector<SQterm>::iterator t_ret = retval.end(); --t_ret;
	for(auto t2 = inTerms.begin();t2 != inTerms.end();){
	  if(isAdditive(*t1, *t2) && t1!=t2) {
	    t1->set_numConst(t1->get_numConst()+t2->get_numConst()); 
	    t2 = inTerms.erase(t2);
	    //t1->masquerade();
	  } // End if
	  else ++t2;
	} // End for
	if(fabs(t1->get_numConst())) retval.push_back(*t1);
      } // End for
      if((int)inTerms.size() != count){
	cout << "Something is wrong in combining terms .... " << endl;
	abort();
      } // End if
      
      return retval;
    }
 
    // *********************************************************
    // Combine terms with isAdditive()==true
    // *********************************************************
    void combineTerms(SQterms &inTerms, SQterms *outTerms)
    {
      int count(0);
      outTerms->reserve(inTerms.size());
      for(auto t1 = inTerms.begin();t1 != inTerms.end();++t1, ++count){
	for(auto t2 = inTerms.begin();t2 != inTerms.end();){
	  if(isAdditive(*t1, *t2) && t1!=t2) { 
	    t1->set_numConst(t1->get_numConst()+t2->get_numConst()); 
	    t2 = inTerms.erase(t2);
	    //t1->masquerade();
	  } // Enf if
	  else ++t2;
	} // End t2
	if(fabs(t1->get_numConst())) outTerms->push_back(*t1);
      } // End for
      //cout << " : " << inTerms.size() << ", " << count << endl;
      if((int)inTerms.size() != count){
	cout << "Something is wrong in combining terms .... " << endl;
	abort();
      } // End if
    }
 
    // *********************************************************
    // Combine terms with isAdditive()==true
    // *********************************************************
    void combineTerms2(SQterms &inTerms, SQterms *outTerms)
    {
      int count(0);
      outTerms->reserve(inTerms.size());
      SQcont<SQterm*> p_inTerms; p_inTerms.reserve(inTerms.size());
      for(auto t1 = inTerms.begin();t1 != inTerms.end();++t1) 
	p_inTerms.push_back(&(*t1));    
      for(auto t1 = p_inTerms.begin();t1 != p_inTerms.end();++t1, ++count){
	for(auto t2 = p_inTerms.begin();t2 != p_inTerms.end();){
	  if(isAdditive(**t1, **t2) && t1!=t2) { 
	    (*t1)->set_numConst((*t1)->get_numConst()+(*t2)->get_numConst()); 
	    t2 = p_inTerms.erase(t2);
	    //t1->masquerade();
	  } // Enf if
	  else ++t2;
	} // End t2
	if(fabs((*t1)->get_numConst())) outTerms->push_back(**t1);
      } // End for
      //cout << " : " << inTerms.size() << ", " << count << endl;
      if((int)p_inTerms.size() != count){
	cout << "Something is wrong in combining terms .... " << endl;
	abort();
      } // End if
    }
 
    // *********************************************************
    // Combine terms with isAdditive()==true
    // This implementation may be the most elegant one, but no warrentry!
    // *********************************************************
    void combineTerms_new(SQterms &inTerms)
    {
      for(auto t1 = inTerms.begin();t1 != inTerms.end();){
	for(auto t2 = inTerms.begin();t2 != inTerms.end();){
	  if(isAdditive(*t1, *t2) && t1!=t2) { 
	    t1->set_numConst(t1->get_numConst()+t2->get_numConst()); 
	    t2 = inTerms.erase(t2);
	    //t1->masquerade();
	  } // Enf if
	  else ++t2;
	} // End t2
	if(!fabs(t1->get_numConst())) t1 = inTerms.erase(t1);
	else ++t1;
      } // End for
    }
 

    // *********************************************************
    // Decompose the RDM to contract core indices
    // *********************************************************
    void decomposeRDMCore(SQterm &inTerm, SQterms *outTerms)
    {
      if(!inTerm.get_isInCanonical()) return;
      SQtensors sfGen_list;
      SQtensors other_list;
      
      SQtensors tensors(inTerm.get_tensors());
      //cout << inTerm << endl;
      for(size_t t = 0;t < tensors.size();++t){
	if(is_sfGen(tensors[t].get_name()) || is_RDM(tensors[t].get_name())) 
	  sfGen_list <= tensors[t];
	else other_list <= tensors[t];
	SharedIndices inds(tensors[t].get_indices());
	// Check whether there is any generic, or non-MO indices 
	bool OK(true);
        for(auto i = inds.begin();i != inds.end();++i)
	  if(!isSF(**i) || (*i)->get_char()==gen) OK = false;
	if(!OK && inds.size()){
	  cout << " >>>> Femto::Core::decomposeRDMCore produces error <<<< " << endl;
	  cout << " >>>> " << inTerm << endl;
	  abort();
	} // End if
      } // End t
      if(!sfGen_list.size()) return;
      
      SQtensor e1(*sfGen_list.begin());
      int o1((int)(e1.get_indices().size())/2);
      
      SQcont<int> CindCre;
      SQcont<int> CindDes;
      for(size_t i = 0;i < (size_t)o1;++i){
	if(e1.get_indices()[i   ]->get_char()==core) CindCre <= i;
	if(e1.get_indices()[i+o1]->get_char()==core) CindDes <= i + o1;      
      } // End i
      
      int nc((int)CindCre.size());
      outTerms->reserve(Nterms());
      if(nc != CindDes.size()) return;
      if(nc == 0) { *outTerms <= inTerm; return; }
      
      // New rank of sfGen
      int newOrder(o1 - nc);
      IIvector perms(makePermutations(nc));
      for(auto perm = perms.begin();perm != perms.end();++perm){
	// Compute all the pairs of contraction ....
	// Example: (conPairs[0][p], conPairs[1][p]) is the p-th contraction pair
	IIvector conPairs;
	Ivector p;
	conPairs.push_back(p);
	conPairs.push_back(p);
	for(size_t i = 0;i < nc;++i){
	  conPairs[0].push_back(CindCre[(*perm)[i]]);
	  conPairs[1].push_back(CindDes[i]);
	} // End i
	
	// Evaluate constant prefactor
	int prefactor = 1;
	Ivector ind;
	for(int i = 0;i < nc;++i){
	  auto i_index(find(ind.begin(), ind.end(), i));
	  
	  if(i_index == ind.end()){
	    ind.push_back(i);
	    int i0(conPairs[0][i]);
	    int i1(conPairs[1][i] - o1);
	    if(i1 == i0){
	      prefactor = prefactor * 2;
	      continue;
	    } // End if
	    // While i1 in conPair[0] ////////////
	    auto i_index(find(conPairs[0].begin(), conPairs[0].end(), i1));
	    for(;i_index!=conPairs[0].end();){
	      ind.push_back((int)(i_index-conPairs[0].begin()));
	      i1 = conPairs[1][(size_t)(i_index-conPairs[0].begin())] - o1;
	      if(i1==i0){ prefactor = prefactor * 2; break; }
	      i_index = find(conPairs[0].begin(), conPairs[0].end(), i1);
	    } // End i_index
	  } // End if
	} // End i
	
	// Initialize the term's tensor list
	SQtensors     tensorList(other_list);
	SharedIndices indexList;
	SharedIndex ptr(NULL);
	for(int i = 0;i < 2*newOrder;++i) indexList.push_back(ptr);
	
	// Polulate the index list for the new excitation operator
	// Also, create a kronecker delta for each contraction pair
	int count(0);
	SQtensor* kDelta_ptr(NULL);
	for(int i1 = 0;i1 < o1;++i1){
	  auto i_index(find(conPairs[0].begin(), conPairs[0].end(), i1));  
	  if(i_index != conPairs[0].end()){
	    int i2 = conPairs[1][(size_t)(i_index-conPairs[0].begin())];
	    SharedIndex ind1(e1.get_indices()[(size_t)i1]);
	    SharedIndex ind2(e1.get_indices()[(size_t)i2]);
	    SharedIndices ind_delta;
	    ind_delta <= ind1;
	    ind_delta <= ind2;
	    tensorList <= kDelta(ind_delta);
	    kDelta_ptr = &(tensorList[tensorList.size()-1]);
	  } // End if
	  else{
	    indexList[count] = e1.get_indices()[i1];
	    int i2(i1 + o1);
	    auto i_index(find(conPairs[1].begin(), conPairs[1].end(), i2));
	    for(;i_index!=conPairs[1].end();) {
	      i2 = conPairs[0][(size_t)(i_index-conPairs[1].begin())] + o1;
	      i_index = find(conPairs[1].begin(), conPairs[1].end(), i2);
	    } // End i_index
          indexList[count+newOrder] = e1.get_indices()[i2];
          ++count;
	  } // End else
	} // End i1
	
	// Ensure that all slots in the index list have been filled
	for(auto ind = indexList.begin();ind != indexList.end();++ind)
	  if(*ind == NULL){
	    cout << "There is at least one unassigned index in the new spin-free unitary group generator" << endl;
	    abort(); 
	  } // End if
	
	// If kDelta has two indices belong to different orbital group, eliminate this term.
	if(kDelta_ptr != NULL){
	  if(kDelta_ptr->get_indices()[0]->get_char() != kDelta_ptr->get_indices()[1]->get_char())
	    continue;
	} // End if
	
	// Add the new excitation operator to the tensor list
	if(indexList.size()!=0) tensorList.push_back(sfGen(indexList));
	
	// Determine the sign 
	IIvector indPairs;
	Ivector q;
	indPairs.push_back(q);      
	indPairs.push_back(q);      
	for(int i = 0;i < o1;++i){
	  indPairs[0].push_back(i);
	  auto i_index(find(conPairs[0].begin(), conPairs[0].end(), i));
	  if(i_index != conPairs[0].end())
	    indPairs[1].push_back(conPairs[1][(size_t)(i_index-conPairs[0].begin())]-o1);
	  else{
	    SharedIndex a(e1.get_indices()[i]);
	    auto a_index(find(indexList.begin(), indexList.end(), a));
	    SharedIndex b(indexList[(size_t)(a_index-indexList.begin())+newOrder]);
	    SharedIndices temp;
	    for(size_t j = o1;j < 2*o1;++j) temp <= e1.get_indices()[j];
	    auto b_index(find(temp.begin(), temp.end(), b));
	    indPairs[1].push_back((int)(b_index-temp.begin()));
	  } // End else
	} // End i
	int nperm(get_num_perms(indPairs[0], indPairs[1]));
	double sign(nperm%2 ? -1 : 1);
	
	*outTerms <= SQterm(sign*prefactor*inTerm.get_numConst(), inTerm.get_Consts(), tensorList); 
      } // End perm
      
    }
    

    // *********************************************************
    // Decompose the RDM to extract the generic indices
    // *********************************************************
    void makeInteractions(SQterm &inTerm, SQterms *outTerms, SQcont<char_state> states)
    {

      int Nstates(states.size());

      SQtensors tensors(inTerm.get_tensors());      
      for(size_t t = 0;t < tensors.size();++t){
	SharedIndices inds(tensors[t].get_indices());
	// Check whether there is any generic, or non-MO indices 
	bool OK(true);
        for(auto i = inds.begin();i != inds.end();++i) if(!isSF(**i)) OK = false;
	if(!OK){
	  cout << " >>>> Femto::Core::makeInteractions produces error [1] <<<< " << endl;
	  abort();
	} // End if
      } // End t
      
      SharedIndices inds(inTerm.get_summedBody());
      if(!inds.size()) return;

      // Generate all the index patterns (convert integers up to pow(Nstates,numgen) into notation system of based Nstates)
      int numgen(0);
      for(auto I = inds.begin();I != inds.end();++I) if((*I)->get_char() == gen) ++numgen;
      if(!numgen){ *outTerms <= inTerm; return; }

      SQcont<SQcont<char_state> > allInteractions; // <core, occ, virt> .... 
      int nummax((int)pow(Nstates,numgen));
      for(int I = 0;I < nummax;++I){
        int max(I);
	int m(numgen - 1);
	SQcont<char_state> thispattern;
	while(m >= 0){
	  int thisnum((int)pow(Nstates,m));
	  int n(max/thisnum);
	  max -= (n*thisnum);
	  for(int ns = 0;ns < Nstates;++ns) if(n == ns) thispattern <= states[ns];
	  //*SVD* if     (n == 0) thispattern <= core;
	  //*SVD* else if(n == 1) thispattern <= act;
	  //*SVD* else if(n == 2) thispattern <= virt;

	  //else{ cout << "!!! ERROR !!!" << endl; abort(); } //*DEBUG
          --m;
	} // End m
	allInteractions <= thispattern;
	////////////////////////////
	// Like this ::           //
	// <<0,0,0,0,0, .... ,0>, //
	//  <0,0,0,0,0, .... ,1>, //
	//  <0,0,0,0,0, .... ,2>, //
	//  ....................  //
	//  <1,1,1,1,1, .... ,1>, //
	//  ....................  //
	//  <2,2,2,2,2, .... ,2>> //
	////////////////////////////
	//*DEBUG* cout << " M : " << m << " ::::::::: TEST :::::::::::::::::::::::::::::::" << endl;
	//*DEBUG* cout << "<numgen, nummax> : " << numgen << ", " << nummax << endl;
	//*DEBUG* (cout <<= thispattern) << endl;
      } // End I
      
      // Extract the generic indices ....
      SharedIndices OrigIs(inTerm.get_summedBody());
      Ivector gen_i;
      int count(0);
      for(auto ind = OrigIs.begin();ind != OrigIs.end();++ind){
	if((*ind)->get_char() == Femto::gen) {++count; gen_i.push_back((size_t)(ind-OrigIs.begin()));}
      } // End count
      if(count != numgen) { cout << " >>>> Femto::Core::decomposeRDMGen produces error [1] <<<< " << endl; abort(); }
      for(auto P = allInteractions.begin();P != allInteractions.end();++P){
	*outTerms <= inTerm;
        SharedIndices CpyIs(outTerms->back().get_summedBody());
	if(OrigIs.size() != CpyIs.size())
	  { cout << " >>>> Femto::Core::decomposeRDMGen produces error [2] <<<< " << endl; abort(); }
	for(size_t num_i = 0;num_i < CpyIs.size();++num_i)
	  if(*(CpyIs[num_i]) != *(OrigIs[num_i])) 
	    { cout << " >>>> Femto::Core::decomposeRDMGen produces error [3] <<<< " << endl; abort(); }
	for(size_t num = 0;num < numgen;++num){
	  // Not yet tested enough!!!!!
	  SharedIndex dum(new SQindex(CpyIs[gen_i[num]]->get_index(), P->at(num), CpyIs[gen_i[num]]->get_isSummed(), CpyIs[gen_i[num]]->get_isExt())); 
#ifdef _GEN_DEBUG
	  cout << " ++ DDDD I(" << num << ") : " << *dum << endl; //*DEBUG* 
	  cout << " -- EEEE I(" << num << ") : " << *(CpyIs[gen_i[num]]) << endl; //*DEBUG* 
#endif
          *(CpyIs[gen_i[num]]) = *dum;
#ifdef _GEN_DEBUG
	  cout << " // AAAA I(" << num << ") : " << *(CpyIs[gen_i[num]]) << endl; //*DEBUG* 
#endif
	} // End num
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Check if the term evidently vanishes or not
	if(outTerms->back().get_isInCanonical()){
	  bool isDead(false);
          for(size_t num_t = 0;num_t < outTerms->back().get_tensors().size();++num_t){
            if(is_RDM(outTerms->back().get_tensors()[num_t].get_name()) || is_sfGen(outTerms->back().get_tensors()[num_t].get_name())){
              SharedIndices inds(outTerms->back().get_tensors()[num_t].get_indices());
	      int order(inds.size()/2);
	      pair<int, int> numCore(0, 0);
	      //numCore.first = 0; numCore.second = 0;
	      for(size_t num_i = 0;num_i < inds.size();++num_i){
		if     (inds[num_i]->get_char() == Femto::virt) isDead = true;
		else if(inds[num_i]->get_char() == Femto::core && num_i < order ) ++numCore.first;
		else if(inds[num_i]->get_char() == Femto::core && num_i >= order) ++numCore.second;
		if(isDead) break;
	      } // End if
	      if(numCore.first != numCore.second) isDead = true;
	    } // End if
	  } // End num_t
	  if(isDead) outTerms->pop_back();
	} // End if
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
      } // End P 

    }


//*SVD*     // *********************************************************
//*SVD*     // Decompose the RDM to extract the generic indices
//*SVD*     // *********************************************************
//*SVD*     void makeInteractions(SQterm &inTerm, SQcont<SQterm> *outTerms)
//*SVD*     {
//*SVD*       if(!inTerm.get_isInCanonical()) return;
//*SVD*       //inTerm.set_summedBody(); //*TEST* 
//*SVD* 
//*SVD*       SQtensors tensors(inTerm.get_tensors());      
//*SVD*       for(size_t t = 0;t < tensors.size();++t){
//*SVD* 	SharedIndices inds(tensors[t].get_indices());
//*SVD* 	// Check whether there is any generic, or non-MO indices 
//*SVD* 	bool OK(true);
//*SVD*         for(auto i = inds.begin();i != inds.end();++i) if(!isSF(**i)) OK = false;
//*SVD* 	if(!OK){
//*SVD* 	  cout << " >>>> Femto::Core::decomposeRDMGen produces error [0] <<<< " << endl;
//*SVD* 	  abort();
//*SVD* 	} // End if
//*SVD*       } // End t
//*SVD*       
//*SVD*       SharedIndices inds(inTerm.get_summedBody());
//*SVD*       if(!inds.size()) return;
//*SVD*       // Generate all the index patterns (convert integers up to pow(3,numgen) into ternary)
//*SVD*       int numgen(0);
//*SVD*       for(auto I = inds.begin();I != inds.end();++I) if((*I)->get_char() == gen) ++numgen;
//*SVD*       SQcont<SQcont<char_state> > core_act_virt; // <core, occ, virt>
//*SVD*       int nummax((int)pow(3,numgen));
//*SVD*       for(int I = 0;I < nummax;++I){
//*SVD*         int max(I);
//*SVD* 	int m(numgen - 1);
//*SVD* 	SQcont<char_state> thispattern;
//*SVD* 	while(m >= 0){
//*SVD* 	  int thisnum((int)pow(3,m));
//*SVD* 	  int n(max/thisnum);
//*SVD* 	  max -= (n*thisnum);
//*SVD* 	  if     (n == 0) thispattern <= core;
//*SVD* 	  else if(n == 1) thispattern <= act;
//*SVD* 	  else if(n == 2) thispattern <= virt;
//*SVD* 	  //else{ cout << "!!! ERROR !!!" << endl; abort(); } //*DEBUG
//*SVD*           --m;
//*SVD* 	} // End m
//*SVD* 	core_act_virt <= thispattern;
//*SVD* 	////////////////////////////
//*SVD* 	// Like this ::           //
//*SVD* 	// <<0,0,0,0,0, .... ,0>, //
//*SVD* 	//  <0,0,0,0,0, .... ,1>, //
//*SVD* 	//  <0,0,0,0,0, .... ,2>, //
//*SVD* 	//  ....................  //
//*SVD* 	//  <1,1,1,1,1, .... ,1>, //
//*SVD* 	//  ....................  //
//*SVD* 	//  <2,2,2,2,2, .... ,2>> //
//*SVD* 	////////////////////////////
//*SVD* 	//*DEBUG* cout << " M : " << m << " ::::::::: TEST :::::::::::::::::::::::::::::::" << endl;
//*SVD* 	//*DEBUG* cout << "<numgen, nummax> : " << numgen << ", " << nummax << endl;
//*SVD* 	//*DEBUG* (cout <<= thispattern) << endl;
//*SVD*       } // End I
//*SVD*       
//*SVD*       // Extract the generic indices ....
//*SVD*       SharedIndices OrigIs(inTerm.get_summedBody());
//*SVD*       Ivector gen_i;
//*SVD*       int count(0);
//*SVD*       for(auto ind = OrigIs.begin();ind != OrigIs.end();++ind){
//*SVD* 	if((*ind)->get_char() == gen) {++count; gen_i.push_back((size_t)(ind-OrigIs.begin()));}
//*SVD*       } // End count
//*SVD*       if(count != numgen) { cout << " >>>> Femto::Core::decomposeRDMGen produces error [1] <<<< " << endl; abort(); }
//*SVD*       for(auto P = core_act_virt.begin();P != core_act_virt.end();++P){
//*SVD* 	*outTerms <= inTerm;
//*SVD*         SharedIndices CpyIs(outTerms->back().get_summedBody());
//*SVD* 	if(OrigIs.size() != CpyIs.size())
//*SVD* 	  { cout << " >>>> Femto::Core::decomposeRDMGen produces error [2] <<<< " << endl; abort(); }
//*SVD* 	for(size_t num_i = 0;num_i < CpyIs.size();++num_i)
//*SVD* 	  if(*(CpyIs[num_i]) != *(OrigIs[num_i])) 
//*SVD* 	    { cout << " >>>> Femto::Core::decomposeRDMGen produces error [3] <<<< " << endl; abort(); }
//*SVD* 	for(size_t num = 0;num < numgen;++num){
//*SVD* 	  // Not yet tested enough!!!!!
//*SVD* 	  SharedIndex dum(new SQindex(CpyIs[gen_i[num]]->get_index(), P->at(num), CpyIs[gen_i[num]]->get_isSummed(), CpyIs[gen_i[num]]->get_isExt())); 
//*SVD* #ifdef _GEN_DEBUG
//*SVD* 	  cout << " ++ DDDD I(" << num << ") : " << *dum << endl; //*DEBUG* 
//*SVD* 	  cout << " -- EEEE I(" << num << ") : " << *(CpyIs[gen_i[num]]) << endl; //*DEBUG* 
//*SVD* #endif
//*SVD*           *(CpyIs[gen_i[num]]) = *dum;
//*SVD* #ifdef _GEN_DEBUG
//*SVD* 	  cout << " // AAAA I(" << num << ") : " << *(CpyIs[gen_i[num]]) << endl; //*DEBUG* 
//*SVD* #endif
//*SVD* 	} // End num
//*SVD* 	///////////////////////////////////////////////////////////////////////////////////////////////////////////
//*SVD* 	// Check if the term evidently vanishes or not
//*SVD* 	if(outTerms->back().get_isInCanonical()){
//*SVD* 	  bool isDead(false);
//*SVD*           for(size_t num_t = 0;num_t < outTerms->back().get_tensors().size();++num_t){
//*SVD*             if(is_RDM(outTerms->back().get_tensors()[num_t].get_name()) || is_sfGen(outTerms->back().get_tensors()[num_t].get_name())){
//*SVD*               SharedIndices inds(outTerms->back().get_tensors()[num_t].get_indices());
//*SVD* 	      int order(inds.size()/2);
//*SVD* 	      pair<int, int> numCore;
//*SVD* 	      numCore.first = 0; numCore.second = 0;
//*SVD* 	      for(size_t num_i = 0;num_i < inds.size();++num_i){
//*SVD* 		if     (inds[num_i]->get_char() == Femto::virt) isDead = true;
//*SVD* 		else if(inds[num_i]->get_char() == Femto::core && num_i < order ) ++numCore.first;
//*SVD* 		else if(inds[num_i]->get_char() == Femto::core && num_i >= order) ++numCore.second;
//*SVD* 		if(isDead) break;
//*SVD* 	      } // End if
//*SVD* 	      if(numCore.first != numCore.second) isDead = true;
//*SVD* 	    } // End if
//*SVD* 	  } // End num_t
//*SVD* 	  if(isDead) outTerms->pop_back();
//*SVD* 	} // End if
//*SVD* 	///////////////////////////////////////////////////////////////////////////////////////////////////////////
//*SVD*       } // End P 
//*SVD* 
//*SVD*     }

 
    // *********************************************************
    // Created 2012/11/26 (not tested yet)
    // *********************************************************
    //SQterm SQterm::operator*(const SQterm &obj)
    SQterm times(SQterm &a, SQterm &b)
    {
      SQterm retval;
      // Merge Tensors_
      SQtensors ten_a(a.get_tensors());
      SQtensors ten_b(b.get_tensors());
      ten_a.insert(ten_a.end(), ten_b.begin(), ten_b.end());
      
      // Rename the dummy indices in b to work with those of a
      Ivector a_dummies(3); // core=0, active=1, virtual=2
      for(size_t num = 0;num < a.get_summedBody().size();++num)
	if(a.get_summedBody()[num]->get_isSummed()){
	  if     (a.get_summedBody()[num]->get_char() == core) ++a_dummies[0];
	  else if(a.get_summedBody()[num]->get_char() == act ) ++a_dummies[1];
	  else if(a.get_summedBody()[num]->get_char() == virt) ++a_dummies[2];
	} // End if
      
      Ivector b_dummies(3);  // core=0, active=1, virtual=2
      for(size_t num = 0;num < b.get_summedBody().size();++num)    
	if(b.get_summedBody()[num]->get_isSummed()){
	  ostringstream stm;
	  if     (b.get_summedBody()[num]->get_char() == core) { 
	    stm << a_dummies[0] + (b_dummies[0]++);
	    b.get_summedBody()[num]->put_index("c"+stm.str());
	  } // End if
	  else if(b.get_summedBody()[num]->get_char() == act ) {
	    stm << a_dummies[1] + (b_dummies[1]++);
	    b.get_summedBody()[num]->put_index("o"+stm.str());
	  } // End if
	  else if(b.get_summedBody()[num]->get_char() == virt) {
	    stm << a_dummies[2] + (b_dummies[2]++);
	    b.get_summedBody()[num]->put_index("v"+stm.str());
	  } // End if
	  else if(b.get_summedBody()[num]->get_char() == gen) {
	    stm << a_dummies[2] + (b_dummies[2]++);
	    b.get_summedBody()[num]->put_index("g"+stm.str());
	  } // End if
	} // End if
      
      retval.set_tensors(ten_a);
      b.masquerade(); // Revert all the names of dummy indices
      
      // Merge Consts_
      SQcont<string> Consts_a(a.get_Consts());
      SQcont<string> Consts_b(b.get_Consts());
      Consts_a.insert(Consts_a.end(), Consts_b.begin(), Consts_b.end());
      retval.set_Consts(Consts_a);
      
      // Merge numConst_
      retval.set_numConst((a.get_numConst())*(b.get_numConst()));
      
      return retval;
    }

    // *********************************************************
    // A convenient function to process terms
    // *********************************************************
    void processTerms(SQterms &inTerms, SQterms *outTerms)
    {
      outTerms->reserve(inTerms.size());
      Femto::for_each(inTerms, boost::bind(&Femto::Core::SQterm::contractkDeltas, _1));
      Femto::Core::screenTerms(inTerms, outTerms); // Screen terms with negligible factor
      Femto::for_each(*outTerms, boost::bind(&Femto::Core::SQterm::transform2RDM, _1)); // Transform sfGen to RDM (only if isInCanonical is true)
    }

    // *********************************************************
    // Exclude all the core indices from the RDMs
    // *********************************************************
    void excludeCore(SQterms &inTerms, SQterms *outTerms, string allSpace)
    {
      cout << endl << " >> Exclude core indices from the RDMs << " << endl;
      cout         << "   >> Check if there is any generic index << " << endl;  
      /////////////////////////// Generalized case //////////////////////////////
      SQcont<char_state> states(Femto::GenerateInteractions(allSpace));
      ///////////////////////////////////////////////////////////////////////////
      SQterms outs;
      outs.reserve(Nterms());

      /////////// Cool way with boost ///////////////
      boost::function<void(SQterm&)> makeInts(boost::bind<void, SQterm&, SQterms*, SQcont<char_state> >(&Femto::Core::makeInteractions, _1, &outs, states));
      Femto::for_each(inTerms, makeInts);
      ///////////////////////////////////////////////
      //*SVD* for(auto t = inTerms.begin();t != inTerms.end();++t){
      //*SVD* 	SQterms batch;
      //*SVD* 	Femto::Core::makeInteractions(*t, &batch, states);
      //*SVD* 	outs += batch; 
      //*SVD* } // End t
      ///////////////////////////////////////////////

      if(outs.size()) {
	inTerms.reserve(Nterms());
	inTerms = outs;
	outs.clear();
	outs.reserve(Nterms());
      } // End if
#ifdef _GEN_DEBUG
      ////////////////// DEBUG ////////////////////////////////
      cout << endl;
      cout << " <-- DDD Content --> " << endl;
      Femto::for_each(inTerms, boost::bind(&Femto::Core::SQterm::masquerade, _1));
      cout << inTerms;
      cout << " <-----------------> " << endl << endl;
      /////////////////////////////////////////////////////////
#endif
      for(auto t = inTerms.begin();t != inTerms.end();++t){
	SQterms batch;
	Femto::Core::decomposeRDMCore(*t, &batch);
	Femto::for_each(batch, boost::bind(&Femto::Core::SQterm::contractkDeltas, _1));
	Femto::for_each(batch, boost::bind(&Femto::Core::SQterm::transform2RDM, _1));
	Femto::for_each(inTerms, boost::bind(&Femto::Core::SQterm::decomposeRDMVirt, _1));   
	outs += batch;
      } // End t
      outTerms->reserve(Nterms());
      processTerms(outs, outTerms);
    }

}} // Femto::core

