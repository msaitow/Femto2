//
//  SQnetwork.hpp
//  
//  Small wraper to call tree container and its utils (tree.hh)
//  Network of tensor contraction can be representable as ContarNet and BinaryNet
//
//  Created by Masaaki Saitow on 14/04/10.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#pragma once

#include <iostream>
#include <SQcontract.hpp>
#include "../tree-2.81/src/tree.hh"
#include "../tree-2.81/src/tree_util.hh"

#define _DEBUG_TREE

namespace Femto { namespace Core {


    // *********************************************************
    // Convert SQcont<...> to tree<...> 
    // *********************************************************
    template<typename T>
    void WeaveNet(SQcont<T> const &incont, tree<T> &outnet)
    {
      outnet.clear();
      auto top(outnet.begin());
      auto now(outnet.insert(top, *incont.cbegin()));
      for(auto c = incont.cbegin()+1;c != incont.cend();++c)
	now = outnet.append_child(now, *c);
    }

    // *********************************************************
    // Convert SQcont<SQcont<...> > to tree<...> (with multiple heads) 
    // *********************************************************
    template<typename T>
    void WeaveNets(SQcont<SQcont<T> > &incont, tree<T> &outnet)
    {
      outnet.clear();
      auto top(outnet.begin());
      for(auto c = incont.begin();c != incont.end();++c){
	if(c->size()){
	  auto now(outnet.insert(top, c->at(0)));
	  for(size_t num = 1;num < c->size();++num)
	    now = outnet.append_child(now, c->at(num));
	} // End if
      } // End c
    }

    // *********************************************************
    // Print the tree structure 
    // *********************************************************
    template<typename T>
    void printTree(tree<T> const &intree, std::ostream& str=std::cout)
    {
      int nhead(0);
      auto now(intree.begin());
      while(intree.is_valid(now)){
	str << " " << *now << boost::format(" <-- [Head](%3d)") % nhead++ << std::endl;
	for(typename tree<T>::iterator sib = intree.begin(now);sib != intree.end(now);++sib){
	  for(int i = 0;i < intree.depth(sib);++i) str << "---";
	  str << "> ";
	  str << *sib << std::endl;
	} // End sib
	now = intree.next_at_same_depth(now);
      } // End while
    }

    // *********************************************************
    // Merge the nodes at same depth of a network
    // *********************************************************
    template<typename T>
    void mergeNodes(tree<T> &intree)
    {
      // First, determine  which top node has the shallowest depth
      int minimal_d(1000);
      for(typename tree<T>::leaf_iterator c = intree.begin_leaf();c != intree.end_leaf();++c){
	int this_depth(intree.depth(c));
	if(this_depth < minimal_d) minimal_d = this_depth;
      } // End c
      if(minimal_d == -1) return;

      typedef std::map<T, typename tree<T>::iterator> Info;
      auto now(intree.begin());
      //const int max_depth(intree.max_depth()-1);
      const int max_depth(minimal_d);
      for(int depth = 0;depth < max_depth;++depth){
        typename tree<T>::fixed_depth_iterator itr_begin(intree.begin_fixed(intree.begin(), depth));
	Info siblings;
	while(intree.is_valid(itr_begin)){
	  int counts(siblings.count(*itr_begin));
	  typename tree<T>::fixed_depth_iterator itr_next(intree.next_at_same_depth(itr_begin));

	  // If current node and the next node are same, combine them
	  if(!counts){
	    typename tree<T>::iterator current(itr_begin);
	    siblings.insert(typename Info::value_type(*itr_begin, current));
	  } // End if
	  else{
	    typename tree<T>::iterator current(itr_begin);
	    typename tree<T>::iterator next(siblings[*itr_begin]);
	    intree.merge(intree.begin(next), intree.end(next), intree.begin(current), intree.end(current));
#ifdef _DEBUG_TREE
	    std::cout << "C: " << *current << " == N: " << *next << std::endl;
#endif
	    typename tree<T>::iterator c_back(itr_begin);
	    const int mydepth(intree.depth(c_back));
	    c_back -= mydepth;
	    intree.erase(c_back); // Erase the redundunt nodes
	  } // End else
	  itr_begin = itr_next;
	} // End while 

      } // End depth
    }

///////////////////////// Not so cool implementation, even though it works ///////////////////////
//*VERIFIED_TO_WORK*//     // *********************************************************
//*VERIFIED_TO_WORK*//     // Merge the nodes at same depth of a network
//*VERIFIED_TO_WORK*//     // *********************************************************
//*VERIFIED_TO_WORK*//     template<typename T>
//*VERIFIED_TO_WORK*//     void mergeNodes(tree<T> &intree)
//*VERIFIED_TO_WORK*//     {
//*VERIFIED_TO_WORK*//       typedef std::map<T, typename tree<T>::iterator> Info;
//*VERIFIED_TO_WORK*//       auto now(intree.begin());
//*VERIFIED_TO_WORK*//       const int max_depth(intree.max_depth()-1);
//*VERIFIED_TO_WORK*//       SQcont<T> toBeErased;
//*VERIFIED_TO_WORK*//       for(int depth = 0;depth < max_depth;++depth){
//*VERIFIED_TO_WORK*//         typename tree<T>::fixed_depth_iterator itr_begin(intree.begin_fixed(intree.begin(), depth));
//*VERIFIED_TO_WORK*// 	Info siblings;
//*VERIFIED_TO_WORK*// 	while(intree.is_valid(itr_begin)){
//*VERIFIED_TO_WORK*// 	  // If current node and the next node are same, combine them
//*VERIFIED_TO_WORK*// 	  if(siblings.count(*itr_begin)){
//*VERIFIED_TO_WORK*// 	    typename tree<T>::iterator current(itr_begin);
//*VERIFIED_TO_WORK*// 	    typename tree<T>::iterator next(siblings[*itr_begin]);
//*VERIFIED_TO_WORK*// 	    //intree.merge(intree.begin(current), intree.end(current), intree.begin(next), intree.end(next));
//*VERIFIED_TO_WORK*// 	    intree.merge(intree.begin(next), intree.end(next), intree.begin(current), intree.end(current));
//*VERIFIED_TO_WORK*// #ifdef _DEBUG_TREE
//*VERIFIED_TO_WORK*// //*TEST*// 	    typename tree<T>::iterator c_next(itr_begin); --c_next;
//*VERIFIED_TO_WORK*// //*TEST*// 	    typename tree<T>::iterator n_next(siblings[*itr_begin]); --n_next;
//*VERIFIED_TO_WORK*// //*TEST*// 	    std::cout << "C: " << *current << "( : " << *c_next << ") == N: " << *next << " ( : "<< *n_next << ")"<< std::endl;
//*VERIFIED_TO_WORK*// 	    std::cout << "C: " << *current << " == N: " << *next << std::endl;
//*VERIFIED_TO_WORK*// #endif
//*VERIFIED_TO_WORK*// 	    //typename tree<T>::iterator Erased(find(intree.begin(), intree.end(), *current));
//*VERIFIED_TO_WORK*// 	    typename tree<T>::iterator c_back(itr_begin); //--c_back;
//*VERIFIED_TO_WORK*// 	    const int mydepth(intree.depth(c_back));
//*VERIFIED_TO_WORK*// 	    c_back -= mydepth;
//*VERIFIED_TO_WORK*// 	    toBeErased <= *c_back;
//*VERIFIED_TO_WORK*// 	    //current = intree.erase(c_back); // Erase the redundunt nodes
//*VERIFIED_TO_WORK*// 	    //typename tree<T>::fixed_depth_iterator temp(current);
//*VERIFIED_TO_WORK*// 	    //itr_begin = temp;
//*VERIFIED_TO_WORK*// 	  } // End if
//*VERIFIED_TO_WORK*// 	  //else {
//*VERIFIED_TO_WORK*// 	    std::cout << *itr_begin << ", " << std::endl;
//*VERIFIED_TO_WORK*// 	    typename tree<T>::iterator current(itr_begin);
//*VERIFIED_TO_WORK*// 	    siblings.insert(typename Info::value_type(*itr_begin, current));
//*VERIFIED_TO_WORK*// 	    //std::cout << " PRE : " <<  << std::endl;
//*VERIFIED_TO_WORK*// 	    itr_begin = intree.next_at_same_depth(itr_begin);
//*VERIFIED_TO_WORK*// 	    if(intree.is_valid(itr_begin)) std::cout << *itr_begin << ":: " << std::endl;
//*VERIFIED_TO_WORK*// 	    else std::cout << "Edge" << std::endl;
//*VERIFIED_TO_WORK*// 	    //}
//*VERIFIED_TO_WORK*// 	} // End while 
//*VERIFIED_TO_WORK*// 	for(auto c = toBeErased.begin();c != toBeErased.end();++c){
//*VERIFIED_TO_WORK*// 	  typename tree<T>::iterator place(std::find(intree.begin(), intree.end(), *c));
//*VERIFIED_TO_WORK*// 	  //std::cout << *place << std::endl;
//*VERIFIED_TO_WORK*// 	  if(intree.is_valid(place)) intree.erase(place);
//*VERIFIED_TO_WORK*// 	} // End c
//*VERIFIED_TO_WORK*// 
//*VERIFIED_TO_WORK*// //*TEST*// /////////////////////////////// TEST/////////////////////////////////////
//*VERIFIED_TO_WORK*// //*TEST*// 	while(intree.is_valid(itr_begin)){
//*VERIFIED_TO_WORK*// //*TEST*// 	  std::cout << *itr_begin << ", ";;
//*VERIFIED_TO_WORK*// //*TEST*// 	  typename tree<T>::iterator current(itr_begin);
//*VERIFIED_TO_WORK*// //*TEST*// 	  ++current;
//*VERIFIED_TO_WORK*// //*TEST*// 	  if(intree.is_valid(current)) std::cout << " I :: " << *current;
//*VERIFIED_TO_WORK*// //*TEST*// 	  else std::cout << "Edge_C "; 
//*VERIFIED_TO_WORK*// //*TEST*// 	  std::cout << ", N :: ";
//*VERIFIED_TO_WORK*// //*TEST*//           itr_begin = intree.next_at_same_depth(itr_begin);
//*VERIFIED_TO_WORK*// //*TEST*// 	  if(intree.is_valid(itr_begin)) std::cout << *itr_begin << ":: " << std::endl;
//*VERIFIED_TO_WORK*// //*TEST*// 	  else std::cout << "Edge" << std::endl;	  
//*VERIFIED_TO_WORK*// //*TEST*// 	}
//*VERIFIED_TO_WORK*// //*TEST*// /////////////////////////////// TEST/////////////////////////////////////
//*VERIFIED_TO_WORK*// 
//*VERIFIED_TO_WORK*//       } // End depth
//*VERIFIED_TO_WORK*//     }
///////////////////////////////////////////////////////////////////////////////////////////////////////

}} // Femto::Core
////////////////////////// Abbreviations ///////////////////////
typedef tree<Femto::Core::SQcontract> ContraNet;
typedef tree<Femto::Core::SQbinary> BinaryNet;
typedef Femto::SQcont<tree<Femto::Core::SQcontract> > ContraNets;
typedef Femto::SQcont<tree<Femto::Core::SQbinary> > BinaryNets;
////////////////////////////////////////////////////////////////
