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
#include "tree.hh"
#include "tree_util.hh"

using namespace std;

class Node{
public:
  Node() {}
  Node(const double i) { val_ = i; }
  double value() { return val_; }
  friend std::ostream &operator<<(std::ostream &os, const Node &i)
  { os << i.val_; return os; }
 
private:
  double val_;
};


int main(int, char **)
{

  tree<Node> gameTree;
  tree<Node>::iterator root = gameTree.insert(gameTree.begin(),Node(14));
  tree<Node>::iterator first = gameTree.append_child(root,Node(32.0));
  tree<Node>::iterator second = gameTree.append_child(root,Node(64.0));
  tree<Node>::iterator root2 = gameTree.insert(gameTree.begin(),Node(37));
  tree<Node>::iterator first2 = gameTree.append_child(root2,Node(43));
  gameTree.append_child(second,Node(21.0));
  gameTree.append_child(second,Node(24.0));

  tree<Node>::leaf_iterator begin = gameTree.begin_leaf();
  tree<Node>::leaf_iterator end = gameTree.end_leaf();

//*WRONG*   tree<Node>::iterator begin = gameTree.begin_leaf();
//*WRONG*   tree<Node>::iterator end = gameTree.end_leaf();

  int x = 0;
  while (begin != end) 
    { 
      cout << begin->value() << endl; 
      begin++;
    }

  cout << endl;
  cout << " ::::::::::::::::::::::::::::::::::::::::::: " << endl;
  kptree::print_tree_bracketed(gameTree);
  cout << endl;
  cout << " ::::::::::::::::::::::::::::::::::::::::::: " << endl << endl;


  // Iterate over all members that belong to the same depth like this!!!!!
  tree<Node>::fixed_depth_iterator itr_bgn = gameTree.begin_fixed(gameTree.begin(), 1);
  while(gameTree.is_valid(itr_bgn)){
    cout << *itr_bgn << endl;
    itr_bgn = gameTree.next_at_same_depth(itr_bgn);
  }

}
