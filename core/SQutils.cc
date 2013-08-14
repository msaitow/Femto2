//
//  SQutils.cc
//  
//
//  Created by Masaaki Saitow on 12/07/01.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <Femto.hpp>

using namespace std;

namespace Femto {

  // *********************************************************
  // Returns n-tuple
  // *********************************************************
  IIvector makePermutations(const int n)
  {
    if(n == 1){
      IIvector outList;
      Ivector temp;
      temp.push_back(0);
      outList.push_back(temp);
      return outList;
    }
    IIvector outList;
    for(size_t i = 0;i < (size_t)n;++i){
      IIvector temp = makePermutations(n-1);
      for(size_t j = 0;j < temp.size();++j){
        for(size_t k = 0;k < temp[j].size();++k){
          if(temp[j][k] >= i) ++temp[j][k];
	} // End k
        outList.push_back(temp[j]); 
        // Looks a little bia tricky, but at least, it works ....
	auto it = outList.end(); --it;
        it->insert(it->begin(), 1, i);
      } // End j
    } // End i
    return outList;
  }

  // *********************************************************
  // Returns n-tuple, composed of elements of (0,1,,...,p) 
  // with redundancy
  // *********************************************************
  IIvector makeTuples1(int n, Ivector &inList)
  {
    IIvector outList;
    int p = (int)inList.size();

    if(n == 0) {
      Ivector temp;
      outList.push_back(temp);
      return outList;
    }
    if(n == 1) {
      for(size_t i = 0;i < p;++i) {
	Ivector temp;
        temp.push_back(inList[i]);
        outList.push_back(temp);
      } // End i      
      return outList;
    } // End if
    if(n == p) {
      Ivector temp;
      for(size_t i = 0;i < p;++i) temp.push_back((int)i);
      outList.push_back(temp);
      return outList;
    } // End i
    if(p < n){
      cout << "Size of p should be larger than that of n" << endl;
      abort();
    } // End if

    Ivector tempList;
    for(size_t i = 0;i < p;++i) tempList.push_back(inList[i]);
    for(size_t i = 0;i < tempList.size();++i){
      Ivector temp2;
      for(size_t t = 0;t < inList.size();++t){
        if(inList[i]!=inList[t]) temp2.push_back(inList[t]);
      }
      IIvector subList = makeTuples1(n-1, temp2);
      for(size_t j = 0;j < subList.size();++j){
        Ivector temp2;
        temp2.push_back(tempList[i]);
        outList.push_back(temp2);
        for(size_t k = 0;k < subList[j].size();++k){
	  auto it = outList.end(); --it;
          it->push_back(subList[j][k]);
	} // End k
      } // End j
    } // End i

    return outList;
  }

  // *********************************************************
  // Returns n-tuple, composed of elements of (0,1,,...,p)
  // without redundancy
  // *********************************************************
  IIvector makeTuples2(int n, Ivector &inList)
  {
    size_t p = inList.size();
    IIvector outList;
    if(n == 0) {
      Ivector temp;
      outList.push_back(temp);
      return outList;
    }
    if(n == (int)p) {
      Ivector temp;
      for(size_t i = 0;i < p;++i) temp.push_back(inList[i]);
      outList.push_back(temp);
      return outList;
    } // End if
    if(n == 1){
      for(size_t i = 0;i < p;++i) {
	Ivector temp;
        temp.push_back(inList[i]);
        outList.push_back(temp);
      } // End i
      return outList;
    } // End if
    if((int)p < n){
      cout << "Size of p should be larger than that of n" << endl;
      abort();
    } // End if
    
    Ivector tempList;
    for(size_t i = 0;i < p-n+1;++i) tempList.push_back(inList[i]);
    for(size_t i = 0;i < tempList.size();++i){
      Ivector temp2;
      for(size_t t = i+1;t < inList.size();++t) temp2.push_back(inList[t]);
      IIvector subList = makeTuples2(n-1, temp2);
      for(size_t j = 0;j < subList.size();++j){
        Ivector temp2;
        temp2.push_back(tempList[i]);
        outList.push_back(temp2);
        for(size_t k = 0;k < subList[j].size();++k){
	  auto it = outList.end(); --it;
          it->push_back(subList[j][k]);
	} // End k
      } // End j
    } // End i

    return outList;
  }

  // *********************************************************
  // Returns n-tuple, composed of elements of (0,1,,...,p) 
  // with partial redundancy
  // *********************************************************
  IIvector makeTuples3(int n, int order, Ivector &inList)
  {
    IIvector outList;
    int p = (int)inList.size();

    if(p%2 != 0){
      cout << "Function makeTupleRDM is designed only for RDM class. So, inList.size()" << endl;
      cout << "should be an even number." << endl;
      abort();
    }

    if(n == 0) return outList;
    if(n == 1) {
      for(size_t i = 0;i < p;++i) {
	Ivector temp;
        temp.push_back(inList[i]);
        outList.push_back(temp);
      } // End i      
      return outList;
    } // End if
    if(n == p) {
      Ivector temp;
      for(size_t i = 0;i < p;++i) temp.push_back((int)i);
      outList.push_back(temp);
      return outList;
    } // End i
    if(p < n){
      cout << "Size of p should be larger than that of n" << endl;
      abort();
    } // End if

    Ivector tempList;
    for(size_t i = 0;i < p;++i) tempList.push_back(inList[i]);
    for(size_t i = 0;i < tempList.size();++i){
      Ivector temp2;
      for(size_t t = 0;t < inList.size();++t){
        if(inList[i]!=inList[t] and (inList[i]>=order ? inList[i]-order : inList[i]+order)!=inList[t]) temp2.push_back(inList[t]);
      }
      IIvector subList = makeTuples3(n-1, order, temp2);
      for(size_t j = 0;j < subList.size();++j){
        Ivector temp2;
        temp2.push_back(tempList[i]);
        outList.push_back(temp2);
        for(size_t k = 0;k < subList[j].size();++k){
	  auto it = outList.end(); --it;
          it->push_back(subList[j][k]);
	} // End k
      } // End j
    } // End i

    return outList;
  }

  // *********************************************************
  // Returns the number of permutations between two Ivectors
  // *********************************************************
  int get_num_perms(vector<int> &ti, vector<int> &bi)
  {
    if(ti.size()!=bi.size()) abort();
    typedef pair<int, int> p_int;
    vector<p_int> x;
    for(size_t i = 0;i < ti.size();++i){
      p_int temp;
      temp.first = ti[i];
      temp.second = bi[i];
      x.push_back(temp);
    } // End i
    sort(x.begin(), x.end(), SecGreat());
    vector<int> y;
    for(size_t i = 0;i < x.size();++i) y.push_back(x[i].first);

    int n_perms = 0;
    for(size_t i = 0;i < y.size();){
      if(y[i] != (int)i){
        int t = y[i];
        y[i]  = y[t];
        y[t]  = t;
        ++n_perms; 
      } // End if
      else ++i;
    } // End i
    return n_perms;
  }

  // *********************************************************
  // Returns factorial up to n
  // *********************************************************
  int fact(const int n){
    if(n==0) return 1;
    int retval = 1;
    for(int i = 0;i < n;++i) retval *= i + 1;
    return retval;
  }

  // *********************************************************
  // Returns symmetry vector for one-body integrals
  // *********************************************************
  Symmetry h1_symm()
  {
    Symmetry h1_symm;

    Femto::Ivector h0; // {0,1}
    h0.push_back(0);
    h0.push_back(1);

    Femto::Ivector h1; // {1,0}
    h1.push_back(1);
    h1.push_back(0);
 
    h1_symm.first.push_back(h0);
    h1_symm.first.push_back(h1);
    h1_symm.second.push_back(1);
    h1_symm.second.push_back(1);

    return h1_symm;
  }

  // *********************************************************
  // Returns symmetry vector for ERI
  // *********************************************************
  Symmetry h2_symm()
  {
    Symmetry V2_symm;

    Femto::Ivector S0; // {0,1,2,3}
    S0.push_back(0);
    S0.push_back(1); 
    S0.push_back(2); 
    S0.push_back(3); 

    Femto::Ivector S1; // {2,1,0,3};
    S1.push_back(2);
    S1.push_back(1); 
    S1.push_back(0); 
    S1.push_back(3); 

    Femto::Ivector S2; // {0,3,2,1};
    S2.push_back(0);
    S2.push_back(3); 
    S2.push_back(2); 
    S2.push_back(1); 

    Femto::Ivector S3; // {1,0,3,2};
    S3.push_back(1);
    S3.push_back(0); 
    S3.push_back(3); 
    S3.push_back(2); 

    V2_symm.first.push_back(S0);
    V2_symm.first.push_back(S1);
    V2_symm.first.push_back(S2);
    V2_symm.first.push_back(S3);

    V2_symm.second.push_back(1);
    V2_symm.second.push_back(1);
    V2_symm.second.push_back(1);
    V2_symm.second.push_back(1);

    return V2_symm;
  }

  // *********************************************************
  // Returns symmetry vector for ERI in Mulliken notation
  // *********************************************************
  Symmetry h2_symmM()
  {
    Symmetry V2_symm;

    Femto::Ivector S0; // {0,1,2,3}
    S0.push_back(0);
    S0.push_back(1); 
    S0.push_back(2); 
    S0.push_back(3); 

    Femto::Ivector S1; // {1,0,2,3};
    S1.push_back(1);
    S1.push_back(0); 
    S1.push_back(2); 
    S1.push_back(3); 

    Femto::Ivector S2; // {0,1,3,2};
    S2.push_back(0);
    S2.push_back(1); 
    S2.push_back(3); 
    S2.push_back(2); 

    Femto::Ivector S3; // {2,3,0,1};
    S3.push_back(2);
    S3.push_back(3); 
    S3.push_back(0); 
    S3.push_back(1); 

    V2_symm.first.push_back(S0);
    V2_symm.first.push_back(S1);
    V2_symm.first.push_back(S2);
    V2_symm.first.push_back(S3);

    V2_symm.second.push_back(1);
    V2_symm.second.push_back(1);
    V2_symm.second.push_back(1);
    V2_symm.second.push_back(1);

    return V2_symm;
  }

  // *********************************************************
  // Returns unit symmetry for 2-index tensor
  // *********************************************************
  Symmetry u2_symm()
  {
    Symmetry u2_symm;
    Femto::Ivector S2_0;
    S2_0.push_back(0);
    S2_0.push_back(1);
    u2_symm.first.push_back(S2_0);
    u2_symm.second.push_back(1);

    return u2_symm;
  }

  // *********************************************************
  // Returns unit symmetry for 4-index tensor
  // *********************************************************
  Symmetry u4_symm()
  {
    Symmetry u4_symm;
    Femto::Ivector S2_0;
    S2_0.push_back(0);
    S2_0.push_back(1);
    S2_0.push_back(2);
    S2_0.push_back(3);
    u4_symm.first.push_back(S2_0);
    u4_symm.second.push_back(1);

    return u4_symm;
  }

  // *********************************************************
  // Returns symmetry for t2 amplitude
  // *********************************************************
  Symmetry t2_symm()
  {
    Symmetry t2_symm;
    Femto::Ivector S2_0; // {0,1,2,3}
    S2_0.push_back(0);
    S2_0.push_back(1);
    S2_0.push_back(2);
    S2_0.push_back(3);

    Femto::Ivector S2_1; // {1.0,3,2}
    S2_1.push_back(1);
    S2_1.push_back(0);
    S2_1.push_back(3);
    S2_1.push_back(2);

    t2_symm.first.push_back(S2_0);
    t2_symm.first.push_back(S2_1);
    t2_symm.second.push_back(1);
    t2_symm.second.push_back(1);

    return t2_symm;
  }

  // *********************************************************
  // Returns unit symmetry for the general n-rank tensor
  // *********************************************************
  Symmetry un_symm(int n)
  {
    Symmetry un_symm;
    vector<int> temp;
    for(int num = 0;num < n;++num)  temp.push_back(num);
    un_symm.first.push_back(temp);
    un_symm.second.push_back(1);
    return un_symm;
  }

  // *********************************************************
  // Returns a stream of star_state
  // *********************************************************
  SQcont<char_state> GenerateInteractions(string input)
  {
    if(*(input.begin()) != '[') { cout << " >>>> Femto::GenerateInteractions produces syntax error [0] <<<< " << endl; abort(); }
    if(*(input.end()-1) != ']') { cout << " >>>> Femto::GenerateInteractions produces syntax error [1] <<<< " << endl; abort(); }
    input.erase(  input.begin());
    input.erase(--input.end()  );

    typedef boost::char_separator<char> char_separator;
    typedef boost::tokenizer<char_separator> tokenizer;

    char_separator sep(",", "", boost::keep_empty_tokens);
    tokenizer states(input, sep);

    SQcont<char_state> output;
    for(auto s = states.begin();s != states.end();++s){

      vector<string> message;
      boost::split(message, *s, boost::is_any_of(" "));

      for(auto ss = message.begin();ss != message.end();++ss){
	// core
	if     (*ss == "c" || *ss == "C") { if(!output.count(Femto::core)) output <= Femto::core; }
	// active
	else if(*ss == "a" || *ss == "A" || 
		*ss == "o" || *ss == "O") { if(!output.count(Femto::act )) output <= Femto::act;  }
	// virtual
	else if(*ss == "v" || *ss == "V") { if(!output.count(Femto::virt)) output <= Femto::virt; }
	else if(*ss != ""){
	  cout << " >>>> Femto::GenerateInteractions produces syntax error [2] <<<< " << endl;
	  cout << " >>>> Can't recognize \"" << *ss << "\""<< endl;
	  abort();
	} // End else
      } // End ss
    } // End s
    return output;

  }

  // *********************************************************
  // Returns current date and time
  // *********************************************************
  string Femto_date()
  {
    time_t timer;
    time(&timer);
    return string(ctime(&timer));
  }

  // *********************************************************
  // Returns a logo
  // *********************************************************
  string Femto_logo(const string s)
  {
    // All these logos were generated from :
    // http://patorjk.com/software/taag/#p=display&f=Graffiti&t=Type%20Something%20
    string retval("");
    retval += s + "       :::::::::: ::::::::::   :::   ::: ::::::::::: :::::::: \n" ;
    retval += s + "      :+:        :+:         :+:+: :+:+:    :+:    :+:    :+: \n" ;
    retval += s + "     +:+        +:+        +:+ +:+:+ +:+   +:+    +:+    +:+  \n" ;
    retval += s + "    :#::+::#   +#++:++#   +#+  +:+  +#+   +#+    +#+    +:+   \n" ;   
    retval += s + "   +#+        +#+        +#+       +#+   +#+    +#+    +#+    \n" ;        
    retval += s + "  #+#        #+#        #+#       #+#   #+#    #+#    #+#     \n" ;         
    retval += s + " ###        ########## ###       ###   ###     ########       \n" ;
    retval += s + "                                                              \n" ; 
    retval += s + ">> An Integrated Toolset For The Automated Tensor Generation <<\n";
    retval += s + "                                                              \n" ; 
    retval += s + "  >> Ver. " +  Femto_ver()  +                                "\n" ;
    retval += s + "  >> Date " +  Femto_date() +                                "\n" ;
    retval += s + "                                                              \n" ; 
    return retval;
  }

//*OLD*   // *********************************************************
//*OLD*   // Returns a logo
//*OLD*   // *********************************************************
//*OLD*   string Femto_logo(const string s)
//*OLD*   {
//*OLD*     // All these logos were generated from :
//*OLD*     // http://patorjk.com/software/taag/#p=display&f=Graffiti&t=Type%20Something%20
//*OLD*     srand((unsigned)time(NULL));
//*OLD*     int num = rand() % 13;
//*OLD*     string retval("");
//*OLD*     if(num == 0){
//*OLD*       retval += s + " ___________                __               \n" ;
//*OLD*       retval += s + " \\_   _____/____    _____ _/  |_  ____      \n" ;
//*OLD*       retval += s + "  |    __)_/ __ \\  /     \\\\   __\\/  _ \\ \n" ;
//*OLD*       retval += s + "  |     \\ \\  ___/ |  Y Y  \\|  | (  <_> )  \n" ;
//*OLD*       retval += s + "  \\___  /  \\___  >|__|_|  /|__|  \\____/   \n" ;
//*OLD*       retval += s + "      \\/       \\/       \\/                \n" ;
//*OLD*     } // Graffiti
//*OLD*     else if(num == 1){
//*OLD*       retval += s + "     ______                  __           \n" ;
//*OLD*       retval += s + "    / ____/___   ____ ___   / /_ ____     \n" ;
//*OLD*       retval += s + "   / /_   / _ \\ / __ `__ \\ / __// __ \\ \n" ;
//*OLD*       retval += s + "  / __/  /  __// / / / / // /_ / /_/ /    \n" ;
//*OLD*       retval += s + " /_/     \\___//_/ /_/ /_/ \\__/ \\____/  \n" ;
//*OLD*     } // Slant
//*OLD*     else if(num == 2){
//*OLD*       retval += s + " __/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\____________________________________________________________________                                   \n" ;   
//*OLD*       retval += s + "  _\\/\\\\\\///////////_____________________________________________________________________                                             \n" ;    
//*OLD*       retval += s + "   _\\/\\\\\\_______________________________________________________/\\\\\\_____________________                                         \n" ;        
//*OLD*       retval += s + "    _\\/\\\\\\\\\\\\\\\\\\\\\\__________/\\\\\\\\\\\\\\\\______/\\\\\\\\\\__/\\\\\\\\\\_____/\\\\\\\\\\\\\\\\\\\\\\______/\\\\\\\\\\____ \n" ;          
//*OLD*       retval += s + "     _\\/\\\\\\///////_________/\\\\\\/////\\\\\\___/\\\\\\///\\\\\\\\\\///\\\\\\__\\////\\\\\\////_____/\\\\\\///\\\\\\__               \n" ;            
//*OLD*       retval += s + "      _\\/\\\\\\_______________/\\\\\\\\\\\\\\\\\\\\\\___\\/\\\\\\_\\//\\\\\\__\\/\\\\\\_____\\/\\\\\\________/\\\\\\__\\//\\\\\\_       \n" ;         
//*OLD*       retval += s + "       _\\/\\\\\\______________\\//\\\\///////____\\/\\\\\\__\\/\\\\\\__\\/\\\\\\_____\\/\\\\\\_/\\\\___\\//\\\\\\__/\\\\\\__            \n" ;            
//*OLD*       retval += s + "        _\\/\\\\\\_______________\\//\\\\\\\\\\\\\\\\\\\\__\\/\\\\\\__\\/\\\\\\__\\/\\\\\\_____\\//\\\\\\\\\\_____\\///\\\\\\\\\\/___    \n" ;              
//*OLD*       retval += s + "         _\\///_________________\\//////////___\\///___\\///___\\///_______\\/////________\\/////_____                                   \n" ;         
//*OLD*     } // Slant Relief
//*OLD*     else if(num == 3){
//*OLD*       retval += s + "   o__ __o__/_                            o                         \n" ;
//*OLD*       retval += s + "  <|    v                                <|>                        \n" ;
//*OLD*       retval += s + "  < >                                    < >                        \n" ;
//*OLD*       retval += s + "   |         o__  __o   \\o__ __o__ __o    |        o__ __o         \n" ;
//*OLD*       retval += s + "   o__/_    /v      |>   |     |     |>   o__/_   /v     v\\        \n" ;
//*OLD*       retval += s + "   |       />      //   / \\   / \\   / \\   |      />       <\\    \n" ;
//*OLD*       retval += s + "  <o>      \\o    o/     \\o/   \\o/   \\o/   |      \\         /   \n" ;
//*OLD*       retval += s + "   |        v\\  /v __o   |     |     |    o       o       o        \n" ;
//*OLD*       retval += s + "  / \\        <\\/> __/>  / \\   / \\   / \\   <\\__    <\\__ __/>  \n" ;
//*OLD*     } // Acrobatic
//*OLD*     else if(num == 4){
//*OLD*       retval += s + "  .----------------.  .----------------.  .----------------.  .----------------.  .----------------.       \n" ;
//*OLD*       retval += s + " | .--------------. || .--------------. || .--------------. || .--------------. || .--------------. |      \n" ;
//*OLD*       retval += s + " | |  _________   | || |  _________   | || | ____    ____ | || |  _________   | || |     ____     | |      \n" ;
//*OLD*       retval += s + " | | |_   ___  |  | || | |_   ___  |  | || ||_   \\  /   _|| || | |  _   _  |  | || |   .'    `.   | |     \n" ;
//*OLD*       retval += s + " | |   | |_  \\_|  | || |   | |_  \\_|  | || |  |   \\/   |  | || | |_/ | | \\_|  | || |  /  .--.  \\  | | \n" ;
//*OLD*       retval += s + " | |   |  _|      | || |   |  _|  _   | || |  | |\\  /| |  | || |     | |      | || |  | |    | |  | |     \n" ;
//*OLD*       retval += s + " | |  _| |_       | || |  _| |___/ |  | || | _| |_\\/_| |_ | || |    _| |_     | || |  \\  `--'  /  | |    \n" ;
//*OLD*       retval += s + " | | |_____|      | || | |_________|  | || ||_____||_____|| || |   |_____|    | || |   `.____.'   | |      \n" ;
//*OLD*       retval += s + " | |              | || |              | || |              | || |              | || |              | |      \n" ;
//*OLD*       retval += s + " | '--------------' || '--------------' || '--------------' || '--------------' || '--------------' |      \n" ;
//*OLD*       retval += s + "  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'       \n" ;                                                            
//*OLD*     }
//*OLD*     else if(num == 5){
//*OLD*       retval += s + "                                                             \n" ;
//*OLD*       retval += s + "  _______________                                  ______    \n" ;
//*OLD*       retval += s + " |          |                 .'. .`. `````|`````.~      ~.  \n" ;
//*OLD*       retval += s + " |______    |______         .'   `   `.    |    |          | \n" ; 
//*OLD*       retval += s + " |          |             .'           `.  |    |          | \n" ;
//*OLD*       retval += s + " |          |___________.'               `.|     `.______.'  \n" ;
//*OLD*       retval += s + "                                                             \n" ;          
//*OLD*     }
//*OLD*     else if(num == 6){
//*OLD*       retval += s + "       :::::::::: ::::::::::   :::   ::: ::::::::::: :::::::: \n" ;
//*OLD*       retval += s + "      :+:        :+:         :+:+: :+:+:    :+:    :+:    :+: \n" ;
//*OLD*       retval += s + "     +:+        +:+        +:+ +:+:+ +:+   +:+    +:+    +:+  \n" ;
//*OLD*       retval += s + "    :#::+::#   +#++:++#   +#+  +:+  +#+   +#+    +#+    +:+   \n" ;   
//*OLD*       retval += s + "   +#+        +#+        +#+       +#+   +#+    +#+    +#+    \n" ;        
//*OLD*       retval += s + "  #+#        #+#        #+#       #+#   #+#    #+#    #+#     \n" ;         
//*OLD*       retval += s + " ###        ########## ###       ###   ###     ########       \n" ;
//*OLD*       retval += s + "                                                              \n" ; 
//*OLD*       retval += s + ">> An Integrated Toolset For The Automated Tensor Generation <<\n";
//*OLD*     }
//*OLD*     else if(num == 7){
//*OLD*       retval += s + " 8 8888888888   8 8888888888            ,8.       ,8.    8888888 8888888888 ,o888888o.     \n" ; 
//*OLD*       retval += s + " 8 8888         8 8888                 ,888.     ,888.         8 8888    . 8888     `88.   \n" ; 
//*OLD*       retval += s + " 8 8888         8 8888                .`8888.   .`8888.        8 8888   ,8 8888       `8b  \n" ;  
//*OLD*       retval += s + " 8 8888         8 8888               ,8.`8888. ,8.`8888.       8 8888   88 8888        `8b \n" ;   
//*OLD*       retval += s + " 8 888888888888 8 888888888888      ,8'8.`8888,8^8.`8888.      8 8888   88 8888         88 \n" ;   
//*OLD*       retval += s + " 8 8888         8 8888             ,8' `8.`8888' `8.`8888.     8 8888   88 8888         88 \n" ;   
//*OLD*       retval += s + " 8 8888         8 8888            ,8'   `8.`88'   `8.`8888.    8 8888   88 8888        ,8P \n" ;    
//*OLD*       retval += s + " 8 8888         8 8888           ,8'     `8.`'     `8.`8888.   8 8888   `8 8888       ,8P  \n" ;     
//*OLD*       retval += s + " 8 8888         8 8888          ,8'       `8        `8.`8888.  8 8888    ` 8888     ,88'   \n" ;          
//*OLD*       retval += s + " 8 8888         8 888888888888 ,8'         `         `8.`8888. 8 8888       `8888888P'     \n" ;          
//*OLD*     }
//*OLD*     else if(num == 8){
//*OLD*       retval += s + " 8888888888                     888                  \n" ;       
//*OLD*       retval += s + " 888                            888                  \n" ;           
//*OLD*       retval += s + " 888                            888                  \n" ;             
//*OLD*       retval += s + " 8888888  .d88b.  88888b.d88b.  888888  .d88b.       \n" ;                   
//*OLD*       retval += s + " 888     d8P  Y8b 888 \"888 \"88b 888    d88\"\"88b  \n" ;               
//*OLD*       retval += s + " 888     88888888 888  888  888 888    888  888      \n" ;               
//*OLD*       retval += s + " 888     Y8b.     888  888  888 Y88b.  Y88..88P      \n" ;       
//*OLD*       retval += s + " 888      \"Y8888  888  888  888  \"Y888  \"Y88P\"   \n" ;         
//*OLD*     }
//*OLD*     else if(num == 9){
//*OLD*       retval += s + " `MMMMMMM                                         \n" ;    
//*OLD*       retval += s + "  MM    \\                         /              \n" ;  
//*OLD*       retval += s + "  MM       ____  ___  __    __   /M      _____    \n" ;       
//*OLD*       retval += s + "  MM   ,  6MMMMb `MM 6MMb  6MMb /MMMMM  6MMMMMb   \n" ;    
//*OLD*       retval += s + "  MMMMMM 6M'  `Mb MM69 `MM69 `Mb MM    6M'   `Mb  \n" ;        
//*OLD*       retval += s + "  MM   ` MM    MM MM'   MM'   MM MM    MM     MM  \n" ;          
//*OLD*       retval += s + "  MM     MMMMMMMM MM    MM    MM MM    MM     MM  \n" ;           
//*OLD*       retval += s + "  MM     MM       MM    MM    MM MM    MM     MM  \n" ;          
//*OLD*       retval += s + "  MM     YM    d9 MM    MM    MM YM.  ,YM.   ,M9  \n" ;   
//*OLD*       retval += s + " _MM_     YMMMM9 _MM_  _MM_  _MM_ YMMM9 YMMMMM9   \n" ;      
//*OLD*     }
//*OLD*     else if(num == 10){
//*OLD*       retval += s + " `7MM\"\"\"YMM                         mm               \n" ;   
//*OLD*       retval += s + "   MM    `7                         MM                  \n" ;           
//*OLD*       retval += s + "   MM   d  .gP\"Ya `7MMpMMMb.pMMMb.mmMMmm ,pW\"Wq.      \n" ;            
//*OLD*       retval += s + "   MM\"\"MM ,M\'   Yb  MM    MM    MM  MM  6W\'   `Wb   \n" ;            
//*OLD*       retval += s + "   MM   Y 8M\"\"\"\"\"\"  MM    MM    MM  MM  8M     M8 \n" ;             
//*OLD*       retval += s + "   MM     YM.    ,  MM    MM    MM  MM  YA.   ,A9       \n" ;                
//*OLD*       retval += s + " .JMML.    `Mbmmd'.JMML  JMML  JMML.`Mbmo`Ybmd9'        \n" ;                 
//*OLD*     }
//*OLD*     else if(num == 11){
//*OLD*       retval += s + "      #                #########      #     #   # \n" ;
//*OLD*       retval += s + " ########## ##########         #   #######  #   # \n" ;
//*OLD*       retval += s + "     #    #         #          #    # #     #   # \n" ;
//*OLD*       retval += s + "     #    #        #   ########     # #     #   # \n" ;
//*OLD*       retval += s + "    #     #     # #           #  ##########    #  \n" ;   
//*OLD*       retval += s + "   #   # #       #            #       #       #   \n" ;     
//*OLD*       retval += s + "  #     #         #    ########       #     ##    \n" ;     
//*OLD*     }
//*OLD*     else if(num == 12){
//*OLD*       retval += s + "     _/_/_/_/                            _/             \n" ;                
//*OLD*       retval += s + "    _/        _/_/    _/_/_/  _/_/    _/_/_/_/    _/_/  \n" ;                  
//*OLD*       retval += s + "   _/_/_/  _/_/_/_/  _/    _/    _/    _/      _/    _/ \n" ;                 
//*OLD*       retval += s + "  _/      _/        _/    _/    _/    _/      _/    _/  \n" ;                    
//*OLD*       retval += s + " _/        _/_/_/  _/    _/    _/      _/_/    _/_/     \n" ;                      
//*OLD*     }
//*OLD*     return retval;
//*OLD*   }


//*OLD*   // *********************************************************
//*OLD*   // Returns a logo
//*OLD*   // *********************************************************
//*OLD*   void Femto_logo()
//*OLD*   {
//*OLD*     // All these logos were generated from :
//*OLD*     // http://patorjk.com/software/taag/#p=display&f=Graffiti&t=Type%20Something%20
//*OLD*     srand((unsigned)time(NULL));
//*OLD*     int num = rand() % 13;
//*OLD*     if(num == 0){
//*OLD*       cout << "! ___________                __          " << endl;
//*OLD*       cout << "! \\_   _____/____    _____ _/  |_  ____  " << endl;
//*OLD*       cout << "!  |    __)_/ __ \\  /     \\\\   __\\/  _ \\ " << endl;
//*OLD*       cout << "!  |     \\ \\  ___/ |  Y Y  \\|  | (  <_> )" << endl;
//*OLD*       cout << "!  \\___  /  \\___  >|__|_|  /|__|  \\____/ " << endl;
//*OLD*       cout << "!      \\/       \\/       \\/              " << endl;
//*OLD*     } // Graffiti
//*OLD*     else if(num == 1){
//*OLD*       cout << "!     ______                  __        " << endl;
//*OLD*       cout << "!    / ____/___   ____ ___   / /_ ____  " << endl;
//*OLD*       cout << "!   / /_   / _ \\ / __ `__ \\ / __// __ \\ " << endl;
//*OLD*       cout << "!  / __/  /  __// / / / / // /_ / /_/ / " << endl;
//*OLD*       cout << "! /_/     \\___//_/ /_/ /_/ \\__/ \\____/  " << endl;
//*OLD*     } // Slant
//*OLD*     else if(num == 2){
//*OLD*       cout << "! __/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\____________________________________________________________________                                   " << endl;   
//*OLD*       cout << "!  _\\/\\\\\\///////////_____________________________________________________________________                                             " << endl;    
//*OLD*       cout << "!   _\\/\\\\\\_______________________________________________________/\\\\\\_____________________                                         " << endl;        
//*OLD*       cout << "!    _\\/\\\\\\\\\\\\\\\\\\\\\\__________/\\\\\\\\\\\\\\\\______/\\\\\\\\\\__/\\\\\\\\\\_____/\\\\\\\\\\\\\\\\\\\\\\______/\\\\\\\\\\____ " << endl;          
//*OLD*       cout << "!     _\\/\\\\\\///////_________/\\\\\\/////\\\\\\___/\\\\\\///\\\\\\\\\\///\\\\\\__\\////\\\\\\////_____/\\\\\\///\\\\\\__               " << endl;            
//*OLD*       cout << "!      _\\/\\\\\\_______________/\\\\\\\\\\\\\\\\\\\\\\___\\/\\\\\\_\\//\\\\\\__\\/\\\\\\_____\\/\\\\\\________/\\\\\\__\\//\\\\\\_       " << endl;         
//*OLD*       cout << "!       _\\/\\\\\\______________\\//\\\\///////____\\/\\\\\\__\\/\\\\\\__\\/\\\\\\_____\\/\\\\\\_/\\\\___\\//\\\\\\__/\\\\\\__            " << endl;            
//*OLD*       cout << "!        _\\/\\\\\\_______________\\//\\\\\\\\\\\\\\\\\\\\__\\/\\\\\\__\\/\\\\\\__\\/\\\\\\_____\\//\\\\\\\\\\_____\\///\\\\\\\\\\/___    " << endl;              
//*OLD*       cout << "!         _\\///_________________\\//////////___\\///___\\///___\\///_______\\/////________\\/////_____                                   " << endl;         
//*OLD*     } // Slant Relief
//*OLD*     else if(num == 3){
//*OLD*       cout << "!   o__ __o__/_                            o                  " << endl;
//*OLD*       cout << "!  <|    v                                <|>                 " << endl;
//*OLD*       cout << "!  < >                                    < >                 " << endl;
//*OLD*       cout << "!   |         o__  __o   \\o__ __o__ __o    |        o__ __o   " << endl;
//*OLD*       cout << "!   o__/_    /v      |>   |     |     |>   o__/_   /v     v\\  " << endl;
//*OLD*       cout << "!   |       />      //   / \\   / \\   / \\   |      />       <\\ " << endl;
//*OLD*       cout << "!  <o>      \\o    o/     \\o/   \\o/   \\o/   |      \\         / " << endl;
//*OLD*       cout << "!   |        v\\  /v __o   |     |     |    o       o       o  " << endl;
//*OLD*       cout << "!  / \\        <\\/> __/>  / \\   / \\   / \\   <\\__    <\\__ __/>  " << endl;
//*OLD*     } // Acrobatic
//*OLD*     else if(num == 4){
//*OLD*       cout << "!  .----------------.  .----------------.  .----------------.  .----------------.  .----------------.  " << endl;
//*OLD*       cout << "! | .--------------. || .--------------. || .--------------. || .--------------. || .--------------. | " << endl;
//*OLD*       cout << "! | |  _________   | || |  _________   | || | ____    ____ | || |  _________   | || |     ____     | | " << endl;
//*OLD*       cout << "! | | |_   ___  |  | || | |_   ___  |  | || ||_   \\  /   _|| || | |  _   _  |  | || |   .'    `.   | | " << endl;
//*OLD*       cout << "! | |   | |_  \\_|  | || |   | |_  \\_|  | || |  |   \\/   |  | || | |_/ | | \\_|  | || |  /  .--.  \\  | | " << endl;
//*OLD*       cout << "! | |   |  _|      | || |   |  _|  _   | || |  | |\\  /| |  | || |     | |      | || |  | |    | |  | | " << endl;
//*OLD*       cout << "! | |  _| |_       | || |  _| |___/ |  | || | _| |_\\/_| |_ | || |    _| |_     | || |  \\  `--'  /  | | " << endl;
//*OLD*       cout << "! | | |_____|      | || | |_________|  | || ||_____||_____|| || |   |_____|    | || |   `.____.'   | | " << endl;
//*OLD*       cout << "! | |              | || |              | || |              | || |              | || |              | | " << endl;
//*OLD*       cout << "! | '--------------' || '--------------' || '--------------' || '--------------' || '--------------' | " << endl;
//*OLD*       cout << "!  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  " << endl;                                                            
//*OLD*     }
//*OLD*     else if(num == 5){
//*OLD*       cout << "!                                                             " << endl;
//*OLD*       cout << "!  _______________                                  ______    " << endl;
//*OLD*       cout << "! |          |                 .'. .`. `````|`````.~      ~.  " << endl;
//*OLD*       cout << "! |______    |______         .'   `   `.    |    |          | " << endl; 
//*OLD*       cout << "! |          |             .'           `.  |    |          | " << endl;
//*OLD*       cout << "! |          |___________.'               `.|     `.______.'  " << endl;
//*OLD*       cout << "!                                                             " << endl;          
//*OLD*     }
//*OLD*     else if(num == 6){
//*OLD*       cout << "!       :::::::::: ::::::::::   :::   ::: ::::::::::: :::::::: " << endl;
//*OLD*       cout << "!      :+:        :+:         :+:+: :+:+:    :+:    :+:    :+: " << endl;
//*OLD*       cout << "!     +:+        +:+        +:+ +:+:+ +:+   +:+    +:+    +:+  " << endl;
//*OLD*       cout << "!    :#::+::#   +#++:++#   +#+  +:+  +#+   +#+    +#+    +:+   " << endl;   
//*OLD*       cout << "!   +#+        +#+        +#+       +#+   +#+    +#+    +#+    " << endl;        
//*OLD*       cout << "!  #+#        #+#        #+#       #+#   #+#    #+#    #+#     " << endl;         
//*OLD*       cout << "! ###        ########## ###       ###   ###     ########       " << endl; 
//*OLD*     }
//*OLD*     else if(num == 7){
//*OLD*       cout << "! 8 8888888888   8 8888888888            ,8.       ,8.    8888888 8888888888 ,o888888o.     " << endl; 
//*OLD*       cout << "! 8 8888         8 8888                 ,888.     ,888.         8 8888    . 8888     `88.   " << endl; 
//*OLD*       cout << "! 8 8888         8 8888                .`8888.   .`8888.        8 8888   ,8 8888       `8b  " << endl;  
//*OLD*       cout << "! 8 8888         8 8888               ,8.`8888. ,8.`8888.       8 8888   88 8888        `8b " << endl;   
//*OLD*       cout << "! 8 888888888888 8 888888888888      ,8'8.`8888,8^8.`8888.      8 8888   88 8888         88 " << endl;   
//*OLD*       cout << "! 8 8888         8 8888             ,8' `8.`8888' `8.`8888.     8 8888   88 8888         88 " << endl;   
//*OLD*       cout << "! 8 8888         8 8888            ,8'   `8.`88'   `8.`8888.    8 8888   88 8888        ,8P " << endl;    
//*OLD*       cout << "! 8 8888         8 8888           ,8'     `8.`'     `8.`8888.   8 8888   `8 8888       ,8P  " << endl;     
//*OLD*       cout << "! 8 8888         8 8888          ,8'       `8        `8.`8888.  8 8888    ` 8888     ,88'   " << endl;          
//*OLD*       cout << "! 8 8888         8 888888888888 ,8'         `         `8.`8888. 8 8888       `8888888P'     " << endl;          
//*OLD*     }
//*OLD*     else if(num == 8){
//*OLD*       cout << "! 8888888888                     888              " << endl;       
//*OLD*       cout << "! 888                            888              " << endl;           
//*OLD*       cout << "! 888                            888              " << endl;             
//*OLD*       cout << "! 8888888  .d88b.  88888b.d88b.  888888  .d88b.   " << endl;                   
//*OLD*       cout << "! 888     d8P  Y8b 888 \"888 \"88b 888    d88\"\"88b  " << endl;               
//*OLD*       cout << "! 888     88888888 888  888  888 888    888  888  " << endl;               
//*OLD*       cout << "! 888     Y8b.     888  888  888 Y88b.  Y88..88P  " << endl;       
//*OLD*       cout << "! 888      \"Y8888  888  888  888  \"Y888  \"Y88P\"   " << endl;         
//*OLD*     }
//*OLD*     else if(num == 9){
//*OLD*       cout << "! `MMMMMMM                                         " << endl;    
//*OLD*       cout << "!  MM    \\                         /              " << endl;  
//*OLD*       cout << "!  MM       ____  ___  __    __   /M      _____    " << endl;       
//*OLD*       cout << "!  MM   ,  6MMMMb `MM 6MMb  6MMb /MMMMM  6MMMMMb   " << endl;    
//*OLD*       cout << "!  MMMMMM 6M'  `Mb MM69 `MM69 `Mb MM    6M'   `Mb  " << endl;        
//*OLD*       cout << "!  MM   ` MM    MM MM'   MM'   MM MM    MM     MM  " << endl;          
//*OLD*       cout << "!  MM     MMMMMMMM MM    MM    MM MM    MM     MM  " << endl;           
//*OLD*       cout << "!  MM     MM       MM    MM    MM MM    MM     MM  " << endl;          
//*OLD*       cout << "!  MM     YM    d9 MM    MM    MM YM.  ,YM.   ,M9  " << endl;   
//*OLD*       cout << "! _MM_     YMMMM9 _MM_  _MM_  _MM_ YMMM9 YMMMMM9   " << endl;      
//*OLD*     }
//*OLD*     else if(num == 10){
//*OLD*       cout << "! `7MM\"\"\"YMM                         mm            " << endl;   
//*OLD*       cout << "!   MM    `7                         MM            " << endl;           
//*OLD*       cout << "!   MM   d  .gP\"Ya `7MMpMMMb.pMMMb.mmMMmm ,pW\"Wq.  " << endl;            
//*OLD*       cout << "!   MM\"\"MM ,M\'   Yb  MM    MM    MM  MM  6W\'   `Wb " << endl;            
//*OLD*       cout << "!   MM   Y 8M\"\"\"\"\"\"  MM    MM    MM  MM  8M     M8 " << endl;             
//*OLD*       cout << "!   MM     YM.    ,  MM    MM    MM  MM  YA.   ,A9 " << endl;                
//*OLD*       cout << "! .JMML.    `Mbmmd'.JMML  JMML  JMML.`Mbmo`Ybmd9'  " << endl;                 
//*OLD*     }
//*OLD*     else if(num == 11){
//*OLD*       cout << "!      #                #########      #     #   # " << endl;
//*OLD*       cout << "! ########## ##########         #   #######  #   # " << endl;
//*OLD*       cout << "!     #    #         #          #    # #     #   # " << endl;
//*OLD*       cout << "!     #    #        #   ########     # #     #   # " << endl;
//*OLD*       cout << "!    #     #     # #           #  ##########    #  " << endl;   
//*OLD*       cout << "!   #   # #       #            #       #       #   " << endl;     
//*OLD*       cout << "!  #     #         #    ########       #     ##    " << endl;     
//*OLD*     }
//*OLD*     else if(num == 12){
//*OLD*       cout << "!     _/_/_/_/                            _/             " << endl;                
//*OLD*       cout << "!    _/        _/_/    _/_/_/  _/_/    _/_/_/_/    _/_/  " << endl;                  
//*OLD*       cout << "!   _/_/_/  _/_/_/_/  _/    _/    _/    _/      _/    _/ " << endl;                 
//*OLD*       cout << "!  _/      _/        _/    _/    _/    _/      _/    _/  " << endl;                    
//*OLD*       cout << "! _/        _/_/_/  _/    _/    _/      _/_/    _/_/     " << endl;                      
//*OLD*     }
//*OLD*   }


} //Femto::

