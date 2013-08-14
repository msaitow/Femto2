/* 

	Cadabra: a field-theory motivated computer algebra system.
	Copyright (C) 2001-2009  Kasper Peeters <kasper.peeters@aei.mpg.de>

   This program is free software: you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation, either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
*/

#include <algorithm>
#include <string>
#include <iostream>
#include "tree_util.hh"
#include "tree.hh"

using namespace std;

int main(int, char **)
{
  tree<string> tr;
  tree<string>::iterator top, one, two, loc, banana;
  
  top=tr.begin();
  one=tr.insert(top, "one");
  two=tr.append_child(one, "two");
  tr.append_child(two, "apple");
  banana=tr.append_child(two, "banana");
  tree<string>::iterator cherry(tr.append_child(banana,"cherry"));
  tr.append_child(two, "peach");
  tr.append_child(one,"three");
  tr.append_child(cherry, "mine");

  tree<string>::iterator here(tr.insert(top, "One"));
  here = tr.append_child(here, "The child");
  tr.append_child(here, "China");
  tree<string>::iterator mom(tr.append_child(here, "Japan"));
  tr.append_child(mom, "MOM");
 
  for(tree<string>::iterator s = tr.begin();s != tr.end();++s)
    cout << ",, ,, " << *s << endl;

  loc=find(tr.begin(), tr.end(), "one");
  if(loc!=tr.end()) {
    tree<string>::sibling_iterator sib=tr.begin(loc);
    for(;sib != tr.end(loc);++sib) 
      { cout << *sib << ", " << tr.depth(sib) << endl; }

    while(sib!=tr.end(loc)) {
      cout << (*sib) << endl;
      ++sib;
    }
    cout << endl;

    tree<string>::iterator sib2(tr.begin(loc));
    for(;sib2 != tr.end(loc);++sib2){
      for(int i=0; i<tr.depth(sib2); ++i) cout << " ";
      cout << (*sib2) << endl;
    } 

//*NOT*     tree<string>::iterator sib2=tr.begin(loc);
//*NOT*     tree<string>::iterator end2=tr.end(loc);
//*NOT*     while(sib2!=end2) {
//*NOT*       for(int i=0; i<tr.depth(sib2); ++i) cout << " ";
//*NOT*       cout << (*sib2) << endl;
//*NOT*       ++sib2;
//*NOT*     }
    cout << endl;

  }

  /////////////////TEST////////////////////////
  cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
  //loc=find(tr.begin(), tr.end(), "one");
  loc=tr.begin();
  while(tr.is_valid(loc)){
    cout << *loc << " <<-- Head!" << endl;
    tree<string>::iterator sib2(tr.begin(loc));
    for(;sib2 != tr.end(loc);++sib2){
      for(int i=0; i<tr.depth(sib2); ++i) cout << " ";
      cout << (*sib2) << endl;
    } 
    loc = tr.next_at_same_depth(loc); // Move to the next head!!
  } // End while
  cout << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
  /////////////////////////////////////////////////

  cout << " ::::::::::::::::::::::::::::::::::::::::::: " << endl;
  kptree::print_tree_bracketed(tr);
  cout << endl;
  cout << " ::::::::::::::::::::::::::::::::::::::::::: " << endl << endl;

  cout << ">>>>>>>>>>>>>> Merge <<<<<<<<<<<<<<<<<" << endl;
  // If one == One, do like this to merge contarction tree!
  loc=find(tr.begin(), tr.end(), "one");
  tree<string>::iterator Loc = find(tr.begin(), tr.end(), "One");
  tr.merge(tr.begin(loc), tr.end(loc), tr.begin(Loc), tr.end(Loc));

  loc=find(tr.begin(), tr.end(), "one");
  if(loc!=tr.end()) {
    tree<string>::sibling_iterator sib=tr.begin(loc);
    for(;sib != tr.end(loc);++sib) 
      { cout << *sib << ", " << tr.depth(sib) << endl; }

    while(sib!=tr.end(loc)) {
      cout << (*sib) << endl;
      ++sib;
    }
    cout << endl;

    tree<string>::iterator sib2(tr.begin(loc));
    for(;sib2 != tr.end(loc);++sib2){
      for(int i=0; i<tr.depth(sib2); ++i) cout << " ";
      cout << (*sib2) << endl;
    } 

//*NOT*     tree<string>::iterator sib2=tr.begin(loc);
//*NOT*     tree<string>::iterator end2=tr.end(loc);
//*NOT*     while(sib2!=end2) {
//*NOT*       for(int i=0; i<tr.depth(sib2); ++i) cout << " ";
//*NOT*       cout << (*sib2) << endl;
//*NOT*       ++sib2;
//*NOT*     }
    cout << endl;
  }

  cout << ">>>>>>>>>>>>>>>>> After Merge <<<<<<<<<<<<<<<<" << endl;
  cout << endl;
  cout << " ::::::::::::::::::::::::::::::::::::::::::: " << endl;
  kptree::print_tree_bracketed(tr);
  cout << endl;
  cout << " ::::::::::::::::::::::::::::::::::::::::::: " << endl << endl;

  cout << ">>>>>>>>>>>>>>>>> Erase <<<<<<<<<<<<<<<<<< " << endl;
  tree<string>::iterator LocE = find(tr.begin(), tr.end(), "One");
  tr.erase(LocE);
  cout << ">>>>>>>>>>>>>>>>> After Erase <<<<<<<<<<<<<<<<" << endl;
  cout << endl;
  cout << " ::::::::::::::::::::::::::::::::::::::::::: " << endl;
  kptree::print_tree_bracketed(tr);
  cout << endl;
  cout << " ::::::::::::::::::::::::::::::::::::::::::: " << endl << endl;

}
