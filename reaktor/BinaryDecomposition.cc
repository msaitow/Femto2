//
//  BinaryDecomposition.cc
//  
//  A bunch of code that transforms a stream of tensor contractions into the binary contractions
//
//  Created by Masaaki Saitow on 13/05/27.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <SQreaktor.hpp>

using namespace std;
using namespace Femto;
using namespace Femto::Core;

#define _DEBUG1
#define _DEBUG2

namespace Femto { namespace Reaktor {

    void BinaryDecomposition(const SQcontracts &inContra, SQcont<SQbinaries> &outBins, Priority prior)
    {
      outBins.reserve(Nterms());
      for(auto Con = inContra.cbegin();Con != inContra.cend();++Con){

        const int numSteps(Con->get_Rtensors().size()-1); // Number of steps required to perfome right-hand side contraction
	typedef vector<tuple<size_t, long unsigned int, long unsigned int> > flop_data;
	vector<flop_data> all_data;
	SQbinaries tempBins; tempBins.reserve(numSteps);            // Binary contarctions

	// Evaluate all possible orders of tensors in the *sequential* contractions
	// < <0,1,2, ..., n>,
	//   <....         > >
	IIvector nCombinations(makePermutations(Con->get_Rtensors().size()));	//cout << "DIM " << nCombinations.size() << endl; //*TEST* 
	for(auto I = nCombinations.begin();I != nCombinations.end();)
	  if(I->at(0) > I->at(1)) I = nCombinations.erase(I);
	  else ++I;

 	// Evaluate all possibilities to form the intermediates and a series of binary contarctions 
	for(auto pos = nCombinations.begin();pos != nCombinations.end();++pos){
	  flop_data X_flops; X_flops.reserve(numSteps); // <Number of step, Size of intermediate, FLOPS>
          SharedIndices prevInds;
	  Ivector perm(*pos);
	  for(size_t num = 1;num < Con->get_Rtensors().size();++num){
	    if(num == 1){ // if first
	      SQtensor t1(Con->get_Rtensors()[perm[0]]);
	      SQtensor t2(Con->get_Rtensors()[perm[1]]);
	      SharedIndices t1Inds(t1.get_indices());
	      SharedIndices t2Inds(t2.get_indices());
	      cout << "t1, " << t1 << " t2, " << t2 << endl;
	      SharedIndices FlopInds(make_union(t1Inds, t2Inds));   // Indices processed in loops
	      SharedIndices IntInds(make_symmdiff(t1Inds, t2Inds)); // Indices associated with the intermediate

#ifdef _DEBUG2
	      {
		cout << " F >> ";
		for(auto I = FlopInds.cbegin();I != FlopInds.cend();++I) cout << **I << " ";
		cout << endl; 
		cout << " I >> ";
		for(auto I = IntInds.cbegin();I != IntInds.cend();++I) cout << **I << " "; 
		cout << endl; 
		cout << " >> CalcO, " << CalcOrder(IntInds) << ", " << CalcOrder(FlopInds) << endl; //*TEST* 
	      }
#endif
	      X_flops.push_back(make_tuple(num, CalcOrder(IntInds), CalcOrder(FlopInds)));
	      prevInds = IntInds;
	    } // End if
	    else{
	      SQtensor t_n(Con->get_Rtensors()[perm[num]]);
	      SharedIndices t_nInds(t_n.get_indices());
	      cout << "t_n, " << t_n << endl;
	      SharedIndices FlopInds(make_union(prevInds, t_nInds));   // Indices processed in loops
	      SharedIndices IntInds(make_symmdiff(prevInds, t_nInds)); // Indices associated with the intermediate

#ifdef _DEBUG2
	      {
		cout << " F >> ";
		for(auto I = FlopInds.cbegin();I != FlopInds.cend();++I) cout << **I << " ";
		cout << endl; 
		cout << " I >> ";
		for(auto I = IntInds.cbegin();I != IntInds.cend();++I) cout << **I << " "; 
		cout << endl; 
		cout << " >> CalcO, " << CalcOrder(IntInds) << ", " << CalcOrder(FlopInds) << endl; //*TEST* 
	      }
#endif
	      X_flops.push_back(make_tuple(num, CalcOrder(IntInds), CalcOrder(FlopInds)));
	      prevInds = IntInds;
	    } // End else
	  } // End num
	  all_data.push_back(X_flops);
	} // End pos

	///////////////////////////////////////////////////////////////////////////
	// Judge which pattern seems to be the best
	///////////////////////////////////////////////////////////////////////////
	// First, determine the least flop count and size of the intermediate
	vector<pair<long unsigned int, long unsigned int> > MaxScores; // Measure the score (and choose the lowest path) 
	for(auto pattern = all_data.cbegin();pattern != all_data.cend();++pattern){
	  pair<long unsigned int ,long unsigned int> Costs(0, 0); // <Size of intermediate, FLOPS>
	  for(auto step = pattern->cbegin();step != pattern->cend();++step){
	    if(get<1>(*step) > Costs.first ) Costs.first  = get<1>(*step);
	    if(get<2>(*step) > Costs.second) Costs.second = get<2>(*step);
	  } // End step
	  MaxScores.push_back(Costs);
	} // End pattern

	pair<int, long unsigned int> int_score (-1, bigNum()); // <Number, lowest score(int)>
	pair<int, long unsigned int> flop_score(-1, bigNum()); // <Number, lowest score(flops)>
	for(auto S = MaxScores.begin();S != MaxScores.end();++S){
	  if(S->first  < int_score.second) 
	    { int_score.first  = (int)(S-MaxScores.begin()); int_score.second  = S->first; } 
	  if(S->second < flop_score.second) 
	    { flop_score.first = (int)(S-MaxScores.begin()); flop_score.second = S->second; } 
	} // End S 
	
	cout << endl;
	cout << "  << Now, we know which one has the minimal score. Then, let's find the dead-heaters >> " << endl;
	vector<int> draw;
	for(auto S = MaxScores.begin();S != MaxScores.end();++S){
	  if(S->first == int_score.second && S->second == flop_score.second) draw.push_back((int)(S-MaxScores.begin()));
	} // End S

	cout << endl;
	cout << " --------------------------------------------------------------------------- " << endl;
	cout << " ||-> Minimal inetrmediate --> " << int_score.first  << endl;
	cout << " ||-> Minimal flop count   --> " << flop_score.first << endl;
	if(draw.size()){
	  cout << " ||-> ";
	  for(auto d = draw.begin();d != draw.end();++d) cout << *d << " ";
	  cout << endl;
	} // End if
	cout << endl;

#ifdef _DEBUG1
	cout << " >> The optimal contraction pattern(s) << " << endl;
	for(auto pattern = all_data.cbegin();pattern != all_data.cend();++pattern){
	  cout << boost::format(" --[%5d]--  |-> ") % (int)(pattern-all_data.cbegin()); 
	  cout << string(Con->get_Rtensors().size()-1, '(');
	  for(size_t num = 0;num < Con->get_Rtensors().size();++num) {
	    cout << Con->get_Rtensors()[nCombinations[(size_t)(pattern-all_data.cbegin())][num]] << (num ? ") " : "");
	  } // End num
	  cout << endl;
	  cout << "      ";
	  for(auto step = pattern->cbegin();step != pattern->cend();++step)
	    cout << boost::format("< [%5d] Xinterm=%15ld, FLOPS=%15ld > ") % ((int)(step-pattern->cbegin())) % get<1>(*step) % get<2>(*step);
	  cout << endl;
	} // End pattern
	cout << endl;
#endif

	{ //*-- Determine and formation of the optimal binary contraction pattern --*//
	  size_t theBest;
	  size_t thepriority;
	  
	  bool hasERI(false); // true, if Con->Rtensors contains at least one ERI
	  SQtensors tens(Con->get_Rtensors());
	  for(auto t = tens.begin();t != tens.end();++t) if(t->get_name() == name_h2()) hasERI = true;
	  
	  // If priority is set to ``ERI", contraction pattern that has ERI in the first contraction is preferred
	  if(prior == ERI && hasERI){
	    cout << " >> Due to the priority, contraction with ERI is preferred." << endl; 
	    if     (Con->get_Rtensors()[nCombinations[flop_score.first][0]].get_name() == name_h2() || Con->get_Rtensors()[nCombinations[flop_score.first][1]].get_name() == name_h2()) theBest = flop_score.first; 
	    else if(Con->get_Rtensors()[nCombinations[int_score.first ][0]].get_name() == name_h2() || Con->get_Rtensors()[nCombinations[int_score.first ][1]].get_name() == name_h2()) theBest = int_score.first; 
	    else{
	      for(auto d = draw.cbegin();d != draw.cend();++d){
		if(Con->get_Rtensors()[nCombinations[*d][0]].get_name() == name_h2() || Con->get_Rtensors()[nCombinations[*d][1]].get_name() == name_h2()) { theBest = *d; break; } 
	      } // End d
	    } // End else 
	  } // End if
	  else if(prior != ERI){
	    if     (prior == Flops)  theBest = flop_score.first; 
	    else if(prior == Interm) theBest = int_score.first;
	    else { cout << " >>>> BinaryContraction : Can't recognize the priority" << endl; abort(); }
	  }
	  else theBest = int_score.first;

	  SQbinaries bins;
	  SharedIndices prevInds;
	  SQtensor X;  
 	  for(size_t num = 1;num < Con->get_Rtensors().size();++num) {
	    SQtensors out_t;
	    if(num == 1){
	      out_t <= Con->get_Rtensors()[nCombinations[theBest][0]];
	      out_t <= Con->get_Rtensors()[nCombinations[theBest][1]];
	      
	      SharedIndices inds1(out_t.at(0).get_indices());
	      SharedIndices inds2(out_t.at(1).get_indices());
	      prevInds = make_symmdiff(inds1, inds2);
	      reorderIndices(prevInds);  // <<-- sortIndices(prevInds), to sort to correct order
	      ostringstream stm;
	      stm << num-1;
	      string name(name_Int() + stm.str());
	      X = SQtensor(name, prevInds, un_symm(prevInds.size()));

	      if(num == Con->get_Rtensors().size()-1){
                SharedIndices LHindices(Con->get_Ltensor().get_indices());
		reorderIndices(LHindices);
		if(LHindices == prevInds) X =Con->get_Ltensor();
	      } // End if
 
	      bins <= SQbinary(Con->get_numConst(), Con->get_Consts(), X, out_t);
	    } // End if
	    else{
	      out_t <= Con->get_Rtensors()[nCombinations[theBest][num]];
	      out_t <= X;

	      SharedIndices inds1(out_t.at(0).get_indices());
	      SharedIndices inds2(out_t.at(1).get_indices());
	      prevInds = make_symmdiff(inds1, inds2);
	      reorderIndices(prevInds); // <<-- sortIndices(prevInds), to sort to correct order

	      ostringstream stm;
	      stm << num-1;
	      string name(name_Int() + stm.str());
	      X = SQtensor(name, prevInds, un_symm(prevInds.size()));

	      if(num == Con->get_Rtensors().size()-1){
                SharedIndices LHindices(Con->get_Ltensor().get_indices());
		reorderIndices(LHindices);
		if(LHindices == prevInds) X =Con->get_Ltensor();
	      } // End if

	      bins <= SQbinary(1.0, UNITY, X, out_t);
	    }
 	  } // End num
	  if(!(bins.back().get_Ltensor() == Con->get_Ltensor())){
            SQtensors RHtensors(X);
	    SQtensor LHtensor(bins.back().get_Ltensor());
	    bins <= SQbinary(1.0, UNITY, LHtensor, RHtensors);
	  } // End if
	  cout << " -->> |||| <<--" << endl << bins << endl; //*TEST* 
 	  outBins <= bins;
	} // End scope
	///////////////////////////////////////////////////////////////////////////

      } // End Con

    }

}} // Femto::Reaktor
