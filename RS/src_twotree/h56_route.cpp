// dev.cpp
// Single-source PACE heuristic solver generated from h56_targeted_medium.cpp.
// Local solver sources and Whidden headers are inlined here for submission readability.
// External dependency: HiGHS headers/library at build time.
#ifndef H56_ROUTE_SECONDS_VALUE
#define H56_ROUTE_SECONDS_VALUE 300.0
#endif
#ifndef H46_HIGHS_SECONDS_VALUE
#define H46_HIGHS_SECONDS_VALUE 900.0
#endif
#ifndef H56_TOTAL_SECONDS_VALUE
#define H56_TOTAL_SECONDS_VALUE 1500.0
#endif
#ifndef H56_OUTPUT_SAFETY_SECONDS_VALUE
#define H56_OUTPUT_SAFETY_SECONDS_VALUE 3.0
#endif
#ifndef H56_HIGHS_SETUP_SECONDS_VALUE
#define H56_HIGHS_SETUP_SECONDS_VALUE 6.0
#endif
#ifndef H56_MIN_HIGHS_SECONDS_VALUE
#define H56_MIN_HIGHS_SECONDS_VALUE 8.0
#endif
#ifndef H56_ILP_PASSES_VALUE
#define H56_ILP_PASSES_VALUE 2
#endif
#ifndef H56_ENABLE_SWAP_REFINEMENT_VALUE
#define H56_ENABLE_SWAP_REFINEMENT_VALUE 0
#endif
#ifndef H56_REENRICH_SECONDS_VALUE
#define H56_REENRICH_SECONDS_VALUE 6.0
#endif
#ifndef H56_MAX_POOL_VALUE
#define H56_MAX_POOL_VALUE 12000
#endif
#ifndef H56_LOCAL_POOL_LIMIT_VALUE
#define H56_LOCAL_POOL_LIMIT_VALUE 5000
#endif
#ifndef H56_MEDIUM_INSTANCE_MIN_LEAVES_VALUE
#define H56_MEDIUM_INSTANCE_MIN_LEAVES_VALUE 900
#endif
#ifndef H56_MEDIUM_MAX_POOL_VALUE
#define H56_MEDIUM_MAX_POOL_VALUE 11000
#endif
#ifndef H56_MEDIUM_CORE_POOL_VALUE
#define H56_MEDIUM_CORE_POOL_VALUE 9000
#endif
#ifndef H56_MEDIUM_ENABLE_EXTERNAL_POOL_VALUE
#define H56_MEDIUM_ENABLE_EXTERNAL_POOL_VALUE 0
#endif
#ifndef H56_ENABLE_SMALL_REGION_SUBSETS_VALUE
#define H56_ENABLE_SMALL_REGION_SUBSETS_VALUE 0
#endif

// ===== BEGIN heuristic/clean_refactor/h56_targeted_medium.cpp =====
#define H46_NO_MAIN

// ===== BEGIN heuristic/clean_refactor/h46.cpp =====
#define main h45_lane_main

// ===== BEGIN heuristic/clean_refactor/h45.cpp =====
#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <chrono>
#include <climits>
#include <csignal>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <fcntl.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

// h25: singleton greedy remains the global incumbent. Inside Whidden, our
// greedy replaces the initial whole-tree and per-cluster approximation used
// to choose clustering and the upward-search starting target.
// Whidden uses C rand() to choose among equally good forests. libc's
// rand() sequence differs between macOS and Linux, so use one portable
// Park-Miller stream. This is also the sequence used by the macOS build that
// produced the stronger development results.
static uint32_t g_portable_rand_state = 1;

static void portable_srand(unsigned seed) {
    constexpr uint32_t modulus = 2147483647u;
    g_portable_rand_state = seed % modulus;
    if (g_portable_rand_state == 0) g_portable_rand_state = 1;
}

static int portable_rand() {
    constexpr uint64_t multiplier = 16807u;
    constexpr uint64_t modulus = 2147483647u;
    g_portable_rand_state =
        (uint32_t)((multiplier * g_portable_rand_state) % modulus);
    return (int)g_portable_rand_state;
}

#define srand portable_srand
#define rand portable_rand

#ifndef H25_GREEDY_SELECTOR_WEIGHT
#define H25_GREEDY_SELECTOR_WEIGHT 2
#endif

class Forest;
class Node;
class ClusterForest;
static int H25_GLOBAL_GREEDY_DISTANCE = -1;
static int h25_local_greedy_distance(
    Forest *first, Forest *second,
    Forest **out_first = nullptr, Forest **out_second = nullptr);
static int h25_greedy_subtree_score(
    Node *candidate, Forest *first, Forest *second);
#ifdef WC4_CLUSTER_ILP_SPLIT
static bool wc4_try_split_ilp(
    Forest *base_first, Forest *base_second, void *candidate_list,
    std::map<int, std::string> *reverse_label_map,
    int cluster_id, int split_id, Forest **out_first, Forest **out_second);
static bool wc4_try_cluster_ilp(
    Forest *base_first, Forest *base_second,
    std::map<int, std::string> *reverse_label_map,
    int cluster_id, ClusterForest *global_first, ClusterForest *global_second,
    int cluster_loc, bool final_cluster,
    Forest **out_first, Forest **out_second);
#endif


// ===== BEGIN heuristic/clean_refactor/h36_whidden/rspr.h =====
/*******************************************************************************
rspr.h

Calculate approximate and exact Subtree Prune and Regraft (rSPR)
distances and the associated maximum agreement forests (MAFs) between pairs
of rooted binary trees.
Supports arbitrary labels. See the
README for more information.

Copyright 2009-2014 Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
April 29, 2014
Version 1.2.2

This file is part of rspr.

rspr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

rspr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rspr.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************/

#define RSPR
//#define DEBUG 1
//#define DEBUG_CONTRACTED 1
//#define DEBUG_APPROX 1
//#define DEBUG_CLUSTERS 1
//#define DEBUG_SYNC 1
//#define DEBUG_UNDO 1
//#define DEBUG_DEPTHS 1
//#define DEBUG_CASE_COUNTER 1
//#define MULT_PICK_LARGEST_GROUP 1

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <climits>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <algorithm>

// ===== BEGIN heuristic/clean_refactor/h36_whidden/Forest.h =====
/*******************************************************************************
Forest.h

Data structure for a forest of binary trees

Copyright 2009-2014 Chris Whidden
cwhidden@dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
March 3, 2014
Version 1.2.1

This file is part of rspr.

rspr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

rspr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rspr.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************/
#ifndef INCLUDE_FOREST

#define INCLUDE_FOREST
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
#include <list>
#include <deque>

// ===== BEGIN heuristic/clean_refactor/h36_whidden/Node.h =====
/*******************************************************************************
Node.h

Data structure for a node of a binary tree
Contains methods to recursively work on a node's subtree

Copyright 2009-2014 Chris Whidden
cwhidden@dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
March 3, 2014
Version 1.2.1

This file is part of rspr.

rspr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

rspr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rspr.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#ifndef INCLUDE_NODE

#define INCLUDE_NODE

#define COPY_CONTRACTED
#define DEBUG_SUPPORT 0

#include <cstdio>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <set>
// skipped already inlined: heuristic/clean_refactor/h36_whidden/Forest.h

using namespace std;

bool IGNORE_MULTI = false;
double REQUIRED_SUPPORT = 0.0;
double MIN_LENGTH = -1.0;

struct StringCompare {
	bool operator() (const string &a, const string &b) const {
		return strcmp(a.c_str(), b.c_str()) < 0;
	}
};


// representation of a component with no leaves
#define DEAD_COMPONENT "*"
//#define DEBUG_PROTECTED "@"
/*void find_sibling_pairs_hlpr(Node *node, list<Node *> &sibling_pairs);
void find_leaves_hlpr(Node *node, vector<Node *> &leaves);
void str_subtree_hlpr(string *s);
list<Node *> find_sibling_pairs(Node *node);
vector<Node *> find_leaves(Node *node);
*/
int stomini(string s);


class Forest;

class Node {
	private:
	//Node *lc;			// left child
	//Node *rc;			// right child
	list<Node *> children; // children
	Node *p;			// parent
	list<Node *>:: iterator p_link;		// location in parents list
	Node *twin;			// counterpart in another tree
	string name;		// label
	int depth;			//distance from root
	int pre_num;	// preorder number
	int edge_pre_start;
	int edge_pre_end;

	int component_number;
	list <Node *> active_descendants;
	list <Node *> root_lcas;
	list<list <Node *>::iterator> removable_descendants;
	list <Node *>::iterator sibling_pair_loc;
	int sibling_pair_status;
	int num_clustered_children;
	Forest *forest;
        list <Node *> contracted_children;
	Node *contracted_lc;
	Node *contracted_rc;
	bool contracted;
	bool edge_protected;
	int max_merge_depth;
	bool allow_sibling;
	int lost_children;
	double support;
	double support_normalization;
        // TODO (Ben): put initialization in correct constructor
        int non_leaf_children = 0;

	public:
	Node() {
		init(NULL, NULL, NULL, "", 0);
	}
	Node(string n) {
		init(NULL, NULL, NULL, n, 0);
	}
	Node(string n, int d) {
		init(NULL, NULL, NULL, n, d);
	}
	Node(Node *lc, Node *rc, Node *p, string n, int d) {
		init(lc, rc, p, n, d);
	}
	void init(Node *lc, Node *rc, Node *p, string n, int d) {
//		this->lc = lc;
//		this->rc = rc;
		this->p = p;
		this->name = string(n);
		this->twin = NULL;
		this->depth = d;
		this->pre_num = -1;
		this->edge_pre_start = -1;
		this->edge_pre_end = -1;
		this->component_number = -2;
		this->active_descendants = list <Node *>();
		this->root_lcas = list <Node *>();
		this->removable_descendants = list< list<Node *>::iterator>();
		this->sibling_pair_loc = list<Node *>::iterator(); 
		this->sibling_pair_status = 0;
		this->num_clustered_children = 0;
		this->forest = NULL;
		this->contracted_lc = NULL;
		this->contracted_rc = NULL;
		this->contracted = false;
		this->edge_protected = false;
		this->allow_sibling = true;
		this->lost_children = 0;
		this->max_merge_depth = -1;
		this->support = -1;
		this->support_normalization = -1;
		this->children = list<Node *>();
		if (lc != NULL)
			add_child(lc);
		if (rc != NULL)
			add_child(rc);
	}
	// copy constructor
	Node(const Node &n) {
		p = NULL;
		name = n.name.c_str();
		twin = n.twin;
		depth = n.depth;
//		depth = 0;
		pre_num = n.pre_num;
		edge_pre_start = n.edge_pre_start;
		edge_pre_end = n.edge_pre_end;
		component_number = n.component_number;
		this->active_descendants = list <Node *>();
		this->removable_descendants = list< list<Node *>::iterator>();
		this->root_lcas = list <Node *>();
		//sibling_pair_loc = n.sibling_pair_loc;
		//sibling_pair_status = n.sibling_pair_status;
		this->sibling_pair_loc = list<Node *>::iterator(); 
		this->sibling_pair_status = 0;
		this->num_clustered_children = 0;
		this->forest = NULL;
		list<Node *>::const_iterator c;
		this->children = list<Node *>();
		for(c = n.children.begin(); c != n.children.end(); c++) {
			add_child(new Node(**c));
		}
#ifdef COPY_CONTRACTED
		if (n.contracted_lc == NULL)
			contracted_lc = NULL;
		else
			contracted_lc = new Node(*(n.contracted_lc), this);
		if (n.contracted_rc == NULL)
			contracted_rc = NULL;
		else
			contracted_rc = new Node(*(n.contracted_rc), this);
		this->contracted = n.contracted;
		for (auto i = n.contracted_children.begin(); i != n.contracted_children.end(); i++) {
		  Node* contracted_node = new Node(**i, this);
		  add_contracted_child(contracted_node);
		  contracted_node->set_parent(this);
		} 
#else
		this->contracted_lc = n.contracted_lc;
		this->contracted_rc = n.contracted_rc;
		this->contracted = n.contracted;
		for (auto i = n.contracted_children.begin(); i != n.contracted_children.end(); i++) {
		  add_contracted_child(**i);
		  contracted_node->set_parent(this);
		} 

#endif
		this->edge_protected = n.edge_protected;
		this->allow_sibling = n.allow_sibling;
		this->lost_children = n.lost_children;
		this->max_merge_depth = n.max_merge_depth;
		this->support = n.support;
		this->support_normalization = n.support_normalization;
		this->non_leaf_children = n.non_leaf_children;
	}

        /* Temporary until rSPR_branch_and_bound_mult_hlpr uses the UndoMachine */
        // copy constructor
        Node(const Node &n, map<Node*, Node*> *node_map) {
		p = NULL;
		name = n.name.c_str();
//		twin = n.twin;
		depth = n.depth;
//		depth = 0;
		pre_num = n.pre_num;
		edge_pre_start = n.edge_pre_start;
		edge_pre_end = n.edge_pre_end;
		component_number = n.component_number;
		this->active_descendants = list <Node *>();
		this->removable_descendants = list< list<Node *>::iterator>();
		this->root_lcas = list <Node *>();
		//sibling_pair_loc = n.sibling_pair_loc;
		//sibling_pair_status = n.sibling_pair_status;
		this->sibling_pair_loc = list<Node *>::iterator(); 
		this->sibling_pair_status = 0;
		this->num_clustered_children = 0;
		this->forest = NULL;
		list<Node *>::const_iterator c;
		this->children = list<Node *>();
		for(c = n.children.begin(); c != n.children.end(); c++) {
		        Node* new_node = new Node(**c, node_map);
			add_child(new_node);
			node_map->emplace(pair<Node*,Node*>(*c, new_node));
		}
#ifdef COPY_CONTRACTED
		if (n.contracted_lc == NULL)
			contracted_lc = NULL;
		else {
		  contracted_lc = new Node(*(n.contracted_lc), this);// node_map);
		  // contracted_lc->set_parent(this);
		  }
		if (n.contracted_rc == NULL)
			contracted_rc = NULL;
		else {
		  contracted_rc = new Node(*(n.contracted_rc), this);//node_map);
			//contracted_rc->set_parent(this);
		}
		this->contracted = n.contracted;
		for (auto i = n.contracted_children.begin(); i != n.contracted_children.end(); i++) {
		  Node* contracted_node = new Node(**i, node_map);
		  contracted_children.push_back(contracted_node);
		  contracted_node->set_parent(this);
		  //TODO: why does this change distance?
		  //node_map->emplace(pair<Node*,Node*>(*i, contracted_node));
		} 
#else
		this->contracted_lc = n.contracted_lc;
		this->contracted_rc = n.contracted_rc;
		this->contracted = n.contracted;
#endif
		this->edge_protected = n.edge_protected;
		this->allow_sibling = n.allow_sibling;
		this->lost_children = n.lost_children;
		this->max_merge_depth = n.max_merge_depth;
		this->support = n.support;
		this->support_normalization = n.support_normalization;
		this->non_leaf_children = n.non_leaf_children;
	}

  
	Node(const Node &n, Node *parent) {
		p = parent;
		name = n.name.c_str();
		twin = n.twin;
		if (p != NULL)
			depth = p->depth+1;
		else
		depth = n.depth;
		pre_num = n.pre_num;
		edge_pre_start = n.edge_pre_start;
		edge_pre_end = n.edge_pre_end;
		component_number = n.component_number;
		this->active_descendants = list <Node *>();
		this->removable_descendants = list< list<Node *>::iterator>();
		this->root_lcas = list <Node *>();
		//sibling_pair_loc = n.sibling_pair_loc;
		//sibling_pair_status = n.sibling_pair_status;
		this->sibling_pair_loc = list<Node *>::iterator(); 
		this->sibling_pair_status = 0;
		this->num_clustered_children = 0;
		this->forest = NULL;
		this->children = list<Node *>();
		list<Node *>::const_iterator c;
		for(c = n.children.begin(); c != n.children.end(); c++) {
			add_child(new Node(**c));
		}
#ifdef COPY_CONTRACTED
		if (n.contracted_lc == NULL)
			contracted_lc = NULL;
		else
			contracted_lc = new Node(*(n.contracted_lc), this);
		if (n.contracted_rc == NULL)
			contracted_rc = NULL;
		else
			contracted_rc = new Node(*(n.contracted_rc), this);
		this->contracted = n.contracted;
		for (auto i = n.contracted_children.begin(); i != n.contracted_children.end(); i++) {
		  Node* contracted_node = new Node(**i, this);
		  add_contracted_child(contracted_node);
		  contracted_node->set_parent(this);
		} 

#else
		this->contracted_lc = n.contracted_lc;
		this->contracted_rc = n.contracted_rc;
		this->contracted = n.contracted;
		for (auto i = n.contracted_children.begin(); i != n.contracted_children.end(); i++) {
		  add_contracted_child(**i);
		  contracted_node->set_parent(this);
		} 

#endif
		this->edge_protected = n.edge_protected;
		this->allow_sibling = n.allow_sibling;
		this->lost_children = n.lost_children;
		this->max_merge_depth = n.max_merge_depth;
		this->support = n.support;
		this->support_normalization = n.support_normalization;
	}
	// TODO: clear_parent function
	~Node() {
		list<Node *>::iterator c = children.begin();
		while(c!= children.end()) {
			Node *n = *c;
			c++;
			n->cut_parent();
		}
		cut_parent();
		active_descendants.clear();
		root_lcas.clear();
		removable_descendants.clear();
#ifdef COPY_CONTRACTED
		if (contracted_lc != NULL) {
			contracted_lc->delete_tree();
		}
		contracted_lc = NULL;
		if (contracted_rc != NULL) {
			contracted_rc->delete_tree();
		}
		contracted_rc = NULL;
#endif
	}
	// TODO: is this still useful?
	/*
	void fake_delete() {
		if (lc != NULL) {
			lc->p = NULL;
			//lc = NULL;
		}
		if (rc != NULL) {
			rc->p = NULL;
			//rc = NULL;
		}
		if (p != NULL) {
			p->delete_child(this);
			//p = NULL;
		}
		active_descendants.clear();
		root_lcas.clear();
		removable_descendants.clear();
	}
	*/

	// delete a subtree
	void delete_tree() {
		list<Node *>::iterator c = children.begin();
		while(c!= children.end()) {
			Node *n = *c;
			c++;
			n->delete_tree();
		}

		cut_parent();
		/*
		children.clear();
		*/

#ifdef COPY_CONTRACTED
		if (contracted_lc != NULL) {
			contracted_lc->delete_tree();
			contracted_lc = NULL;
		}
		if (contracted_rc != NULL) {
			contracted_rc->delete_tree();
			contracted_rc = NULL;
		}
		for (auto i = contracted_children.begin(); i != contracted_children.end(); i++) {
		  (*i)->delete_tree();
		}
		contracted_children.clear();
#endif
		delete this;
	}

	// cut edge between parent and child
	// should really be cut_child, no deleting occurs
	// TODO: is this useful? The only reason for this would be to
	// be sure it works when the child is not correctly set
	void delete_child(Node *n) {
		n->cut_parent();
	}

	// TODO: make sure this doesn't break things with >2 children
	// add a child
	void add_child(Node *n) {
		if (n->p != NULL)
			n->cut_parent();
		n->p_link = children.insert(children.end(),n);
		n->p = this;
		n->depth = depth+1;
		n->contracted = false;
	}

	// TODO: make sure this doesn't break things with >2 children
	// add a child
	void add_child_keep_depth(Node *n) {
		if (n->p != NULL)
			n->cut_parent();
		n->p_link = children.insert(children.end(),n);
		n->p = this;
		n->contracted = false;
	}


	/* TODO: need new method of putting a child in a specific spot
	   in the children list

		 maybe adjacent to an iterator?
	*/

	// insert a child before the given sibling
	 void insert_child(Node *sibling, Node *n) {
		if (n->p != NULL)
			n->cut_parent();
		n->p_link = children.insert(sibling->p_link, n);
		n->depth = depth+1;
		n->p = this;
		n->contracted = false;
	 }

	// insert a child before the given sibling
	 void insert_child_keep_depth(Node *sibling, Node *n) {
		if (n->p != NULL)
			n->cut_parent();
	 	n->p_link = children.insert(sibling->p_link, n);
		n->p = this;
		n->contracted = false;
	 }


	/* dangerous, not relevant to multifurcating
	Node *set_lchild(Node *n) {
		this->lc = n;
		if (n != NULL) {
			n->p = this;
			n->depth = depth+1;
		}
		return lc;
	}

	Node *set_lchild_keep_depth(Node *n) {
		this->lc = n;
		if (n != NULL) {
			n->p = this;
		}
		return lc;
	}
	Node *set_rchild(Node *n) {
		this->rc = n;
		if (n != NULL) {
			n->p = this;
			n->depth = depth+1;
		}
		return rc;
	}
	Node *set_rchild_keep_depth(Node *n) {
		this->rc = n;
		if (n != NULL) {
			n->p = this;
		}
		return rc;
	}
*/
	// potentially dangerous
	Node *set_parent(Node *n) {
		p = n;
		return p;
	}
  
        void add_children(list<Node *> *nodes) {
	    list<Node *>::iterator i;
	    for (i = nodes->begin(); i != nodes->end(); i++) {
	      //children.push_back((*i));
	      add_child(*i);
	    }
	}
	Node *set_twin(Node *n) {
		twin = n;
		return twin;
	}
	void set_name(string n) {
		name = string(n);
	}
	int set_depth(int d) {
		depth = d;
		return depth;
	}
	void fix_depths() {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->depth = depth+1;
			(*c)->fix_depths();
		}
	}
	int set_preorder_number(int p) {
		pre_num = p;
		return pre_num;
	}

	int set_edge_pre_start(int p) {
		edge_pre_start= p;
		return edge_pre_start;
	}

	int set_edge_pre_end(int p) {
		edge_pre_end= p;
		return edge_pre_end;
	}
        void set_non_leaf_children(int p) {
	        non_leaf_children = p;
	}

        void decrement_non_leaf_children() {
	        non_leaf_children--;
	}
        void increment_non_leaf_children() {
	        non_leaf_children++;
	}


	void copy_edge_pre_interval(Node *n) {
		if (n->edge_pre_start > -1) {
			edge_pre_start = n->edge_pre_start;
		}
		if (n->edge_pre_end > -1) {
			edge_pre_end = n->edge_pre_end;
		}
	}
	int set_component_number(int c) {
		component_number = c;
		return c;
	}
	list<Node *>& get_children() {
		return children;
	}
	Node *get_contracted_lc() {
		return contracted_lc;
	}
	Node *get_contracted_rc() {
		return contracted_rc;
	}

	void set_contracted_lc(Node *n) {
		contracted_lc = n;
	}
	void set_contracted_rc(Node *n) {
		contracted_rc = n;
	}


	bool is_protected() {
		return edge_protected;
	}
	bool is_contracted() {
		return contracted;
	}

	void protect_edge() {
		edge_protected = true;
	}

	void unprotect_edge() {
		edge_protected = false;
	}

	void unprotect_subtree() {
		unprotect_edge();
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->unprotect_edge();
		}
	}

	void protect_supported_edges() {
		if (support > 0)
			edge_protected = true;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->protect_supported_edges();
		}
	}

	bool can_be_sibling() {
		return allow_sibling;
	}

	void disallow_siblings() {
		allow_sibling = false;
	}

	void disallow_siblings_subtree() {
		allow_sibling = false;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->disallow_siblings_subtree();
		}
	}

	void allow_siblings() {
		allow_sibling = true;
	}

	void allow_siblings_subtree() {
		allow_sibling = true;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->allow_siblings_subtree();
		}
	}

        bool is_sibling_of(Node* other) {
	  if (parent() != NULL && other->parent() != NULL) {	    
	    return parent() == other->parent();
	  }
	  return false;
	}


	int num_lost_children() {
		return lost_children;
	}

	int get_max_merge_depth() {
		return max_merge_depth;
	}
	void set_max_merge_depth(int d) {
		max_merge_depth = d;
	}


	int count_lost_children_subtree() {
		int lost_children_count = lost_children;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			lost_children_count += (*c)->count_lost_children_subtree();
		}
		return lost_children_count;
	}
	int count_lost_subtree() {
		int lost_children_count = (lost_children > 0) ? 1 : 0;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			lost_children_count += (*c)->count_lost_subtree();
		}
		return lost_children_count;
	}


	void lost_child() {
		lost_children++;
	}

	void no_lost_children() {
		lost_children = 0;
	}

	void no_lost_children_subtree() {
		lost_children = 0;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->no_lost_children_subtree();
		}
	}

	double get_support() {
		return support;
	}
	void set_support(double s) {
		support = s;
	}
	void a_inc_support() {
#pragma omp atomic
		support += 1;
	}
	void a_dec_support() {
#pragma omp atomic
		support -= 1;
	}
	double get_support_normalization() {
		return support_normalization;
	}
	void set_support_normalization(double s) {
		support_normalization = s;
	}
	void a_inc_support_normalization() {
#pragma omp atomic
		support_normalization += 1;
	}
	void a_dec_support_normalization() {
#pragma omp atomic
		support_normalization -= 1;
	}

	void normalize_support() {
		if (support_normalization != 0)
			support /= support_normalization;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->normalize_support();
		}
	}

	int get_component_number() {
		return component_number;
	}
	void increase_clustered_children() {
		num_clustered_children++;
	}
	void decrease_clustered_children() {
		num_clustered_children--;
	}
	void set_num_clustered_children(int c) {
		num_clustered_children = c;
	}
	int get_num_clustered_children() {
		return num_clustered_children;
	}
	list <Node *> *get_active_descendants() {
		return &active_descendants;
	}
	list <Node *> *get_root_lcas(){
		return &root_lcas;
	}
	int get_sibling_pair_status(){
		return sibling_pair_status;
	}
	void set_sibling_pair_status(int s){
		sibling_pair_status = s;
	}
	void set_forest(Forest *f) {
		forest = f;
	}
	void set_forest_rec(Forest *f) {
		forest = f;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->set_forest_rec(f);
		}
	}
	Forest *get_forest() {
		return forest;
	}
	list<list <Node *>::iterator> *get_removable_descendants() {
		return &removable_descendants;
	}
	void initialize_component_number(int value) {
		component_number = value;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->initialize_component_number(value);
		}
	}
	void initialize_active_descendants(list <Node *> value) {
		active_descendants = value;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->initialize_active_descendants(value);
		}
	}
	void initialize_root_lcas(list <Node *> value) {
		root_lcas = value;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->initialize_root_lcas(value);
		}
	}
	void initialize_removable_descendants(list<list <Node *>::iterator> value) {
		removable_descendants = value;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->initialize_removable_descendants(value);
		}
	}


	/* contract:
	 * if this node has a parent and one child then contract it out
	 * if this node has a parent and no child then contract it out
	 *   parent will have one child so contract it as well.
	 * if this node has no parent and one child then take the
	 *	child's children and contract it out
	 * return the first degree two parent found or NULL if there
	 * was no contraction
	 */
	// TODO: check handling twins for interior nodes
	Node *contract(bool remove) {
		Node *parent = p;
		Node *child;
		Node *ret = NULL;
		// contract out this node and give child to parent
		if (parent != NULL) {
		//cout << p->str_subtree() << endl;
			if (children.size() == 1) {
					child = children.front();
					if (this == parent->children.back()) {
						parent->add_child_keep_depth(child);
						child->set_depth(depth);
					}
					else {
						list<Node *>::iterator sib = p_link;
						sib++;
						Node *sibling = *sib;
						parent->insert_child_keep_depth(sibling, child);
						child->set_depth(depth);
					}
					child->copy_edge_pre_interval(this);
					if (edge_protected && !child->is_protected())
						child->protect_edge();
					cut_parent();
					if (remove)
						delete this;
					ret = parent;
			}
			else if (children.empty()) {
				cut_parent();
				ret = parent->contract(remove);
				if (remove)
					delete this;
				//this->fake_delete();
			}
			else
				ret = this;

		}
		// if no parent then take children of single child and remove it
		else {

			// dead component or singleton, will be cleaned up by the forest
			if (children.empty()) {
				if (str() == "")
					name = DEAD_COMPONENT;
			}
			if (children.size() == 1) {
					child = children.front();
					child->cut_parent();

				/* cluster hack - if we delete a cluster node then
				 * we may try to use it later. This only happens once
				 * per cluster so we can spend linear time to update
				 * the forest
				 */
				if (child->num_clustered_children > 0) {
					// TODO: redundant?
					delete_child(child);
					if (remove)
						delete this;
					//this->fake_delete();
					ret = child;
				}
				else {
					// if child is a leaf then get rid of this so we don't lose refs
					// problem: if the child is not c, then we want to copy
					// otherwise we don't
					// copy other parameters and join the twin
					//to this if the child is a label
					Node *new_lc = child->lchild();
					Node *new_rc = child->rchild();
					if (child->is_leaf()) {
						if (child->get_twin() != NULL) {
							set_twin(child->get_twin());
							child->get_twin()->set_twin(this);
						}
						name = child->get_name().c_str();
//						name = child->str();
					}
					// TODO: redundant?
					child->cut_parent();
					list<Node *>::iterator c = child->children.begin();
					while(c!= child->children.end()) {
						Node *new_child = *c;
						c++;
						new_child->cut_parent();
						add_child(new_child);
					}
					if (child->contracted_lc != NULL)
						contracted_lc = child->contracted_lc;
					if (child->contracted_rc != NULL)
						contracted_rc = child->contracted_rc;
					pre_num = child->get_preorder_number();

					c = child->contracted_children.begin();
					while(c!= child->contracted_children.end()) {
						Node *new_child = *c;
						c++;
						new_child->set_parent(this);
						add_contracted_child(new_child);
						contracted = true;
					}
					child->contracted_children.clear();
					if (remove) {
						child->contracted_lc = NULL;
						child->contracted_rc = NULL;
						delete child;
					}
					ret = this;
				}
			}
		}

		return ret;
	}

	Node *contract() {
		return contract(false);
	}

	// TODO: binary only

	/* contract_sibling_pair:
	 * if this node has two child leaves then contract them out
	 * return true if contracted, otherwise false
	 */
	bool contract_sibling_pair() {
		if (lchild() != NULL && lchild()->is_leaf()
				&& rchild() != NULL && rchild()->is_leaf()) {
			#ifdef DEBUG
				string new_name = "<" + lchild()->str() + "," + rchild()->str() + ">";
			#else
				string new_name = "(" + lchild()->str() + "," + rchild()->str() + ")";
			#endif
			set_name(new_name);
			lchild()->cut_parent();
			rchild()->cut_parent();
			return true;
		}
		return false;
	}

	// TODO: binary only
	bool contract_sibling_pair_undoable() {
		if (lchild() != NULL && lchild()->is_leaf()
				&& rchild() != NULL && rchild()->is_leaf()) {
			/*
			#ifdef DEBUG
				string new_name = "<" + lc->str() + "," + rc->str() + ">";
			#else
				string new_name = "(" + lc->str() + "," + rc->str() + ")";
			#endif
			set_name(new_name);
			*/
			Node *lc = lchild();
			Node *rc = rchild();
			contracted_lc = lc;
			contracted_rc = rc;
			rc->cut_parent();
			lc->cut_parent();
			contracted_lc->contracted = true;
			contracted_rc->contracted = true;
			edge_protected = false;
			return true;
		}
		return false;
	}

	/* contract_sibling_pair_undoable
	 * works with multifurcating trees
	 * returns NULL for no contract, otherwise returns the contracted
	 * parent of the nodes
	 */
	Node *contract_sibling_pair_undoable(Node *child1, Node *child2) {

		if (child1->parent() != this ||
				child2->parent() != this)
			return NULL;
		if (children.size() == 2) {
			contract_sibling_pair_undoable();
			return this;
		}
		else {
			Node *new_child = new Node();
			// buggy
//			new_child->set_preorder_number(pre_num);
			if (child1->get_preorder_number() < child2->get_preorder_number()) {
				new_child->set_preorder_number(child1->get_preorder_number());
			}
			else {
				new_child->set_preorder_number(child2->get_preorder_number());
			}
			add_child(new_child);
			new_child->add_child(child1);
			new_child->add_child(child2);
			edge_protected = false;
			//new_child->contract_sibling_pair_undoable();
			return new_child;
		}
	}

	// TODO: binary only
	void undo_contract_sibling_pair() {
		// hacky, might hide problems
		if (contracted_lc != NULL)
			add_child(contracted_lc);
		if (contracted_rc != NULL)
			add_child(contracted_rc);
		contracted_lc = NULL;
		contracted_rc = NULL;
	}

        int recalculate_non_leaf_children() {
	        int total = 0;
		list<Node *>::iterator i;
		for (i = children.begin(); i != children.end(); i++) {
		  if (!(*i)->is_leaf()) {
			total++;
		  }
		}
		non_leaf_children = total;
		return total;
	}

        void add_contracted_child(Node *node) {
                contracted_children.push_back(node);
		node->contracted = true;
		node->set_parent(this);
	}

        list<Node *> *get_combined_children(list<Node *> *nodes) {
	        list<Node *> *combined = new list<Node *>();
		list<Node *>::iterator i;
		for (i = nodes->begin(); i != nodes->end(); i++) {
		  list<Node *> children = (*i)->get_children();
		  combined->insert(combined->end(), children.begin(), children.end());
		}
		return combined;
	}

  
        Node *contract_sibling_group(list<Node *> *nodes) {
	  //check to make sure all nodes are siblings
	  list<Node *>::iterator i;
	  for (i = nodes->begin(); i != nodes->end(); i++) {
	    if ((*i)->parent() != this) {
	      return NULL;
	    }
	  }
	  //contract all children into this
	  if (nodes->size() == children.size()) {
	    //TODO: Both have references to children, double free?
	    list<Node*> *new_children = get_combined_children(nodes);
	    children.clear();
	    for (i = new_children->begin(); i != new_children->end(); i++) {
	      children.push_back(*i);
	    }
	    delete new_children;
	    int min_pre_num = 0x7FFFFFFF; //max value constant?
	    for (i = nodes->begin(); i != nodes->end(); i++) {
	      //children.remove(*i); already removed in setting children
	      add_contracted_child(*i);
	      if ((*i)->get_preorder_number() < min_pre_num) {
		  min_pre_num = (*i)->get_preorder_number();
	      }
	    }
	    recalculate_non_leaf_children();
	    return this;
	  }
	  //Otherwise we create a new node and contract them into that
	  else {
	    Node *new_child = new Node();
	    add_child(new_child);
	    int min_pre_num = 0x7FFFFFFF; //max value constant?
	    for (i = nodes->begin(); i != nodes->end(); i++) {
	      new_child->add_contracted_child((*i));
	      children.remove((*i));
	      if ((*i)->get_preorder_number() < min_pre_num) {
		  min_pre_num = (*i)->get_preorder_number();
		}
	    }
	    new_child->set_preorder_number(min_pre_num);

	    
	    list<Node*> *new_children = get_combined_children(nodes);
	     new_child->add_children(new_children);
	    for (i = new_children->begin(); i != new_children->end(); i++) {
	      (*i)->set_parent(new_child);

	    }
	    delete new_children;
	    new_child->set_parent(this);
	    recalculate_non_leaf_children();
	    new_child->recalculate_non_leaf_children();
	    return new_child;
	  }
	}

    //Takes a list of nodes in the other tree to contract in this tree
    Node *contract_twin_group(list<Node *> *TO_nodes) {
	  list<Node *> TT_nodes;
	    list<Node *>::iterator i;
		for( i = TO_nodes->begin(); i != TO_nodes->end(); i++) {
		    TT_nodes.push_back((*i)->get_twin());			
		}
		return contract_sibling_group(&TT_nodes);		
    }
  
	void fix_contracted_order() {
		if (twin != NULL && twin->contracted_lc->twin != contracted_lc) {
			Node *swap = contracted_lc;
			contracted_lc = contracted_rc;
			contracted_rc = swap;
		}
	}


	// cut the edge between this node and its parent
	void cut_parent() {
		if (p != NULL) {
			// TODO hacky: fix this to use a multi list for contractions
			if (!contracted) {
				p->children.erase(p_link);
			}
			else {
				if (p->contracted_lc == this)
					p->contracted_lc = NULL;
				if (p->contracted_rc == this)
					p->contracted_rc = NULL;
			}
			p = NULL;
			p_link = children.end();
		}
	}

	// caution: destructive
	void contract_node() {
		if (p == NULL || is_leaf())
			return;
		list<Node *>::iterator c = children.begin();
		while(c!= children.end()) {
			Node *n = *c;
			c++;
			p->add_child(n);
		}

#ifdef COPY_CONTRACTED
		if (contracted_lc != NULL) {
			p->add_child(contracted_lc);
		}
		if (contracted_rc != NULL) {
			p->add_child(contracted_rc);
		}
#endif
		delete this;
	}

	Node *parent() {
		return p;
	}
	inline Node *lchild() {
		if (children.empty())
			return NULL;
		else
			return children.front();
	}

	Node *rchild() {
		if (children.empty())
				return NULL;
		list<Node *>::iterator c = ++(children.begin());
		if (c == children.end())
			return NULL;
		else
			return *c;
	}
	Node *get_twin() {
		return twin;
	}
	int get_depth() {
		return depth;
	}
	int get_preorder_number() {
		return pre_num;
	}
	int get_edge_pre_start() {
		return edge_pre_start;
	}
	int get_edge_pre_end() {
		return edge_pre_end;
	}
	string str() {
		string s = "";
		str_hlpr(&s);
		return s;
	}
	string get_name() {
		return name;
	}
  
	void str_hlpr(string *s) {
		if (!name.empty())
			*s += name.c_str();
		if (contracted_lc != NULL || contracted_rc != NULL || !contracted_children.empty()) {
			#ifdef DEBUG_CONTRACTED
				*s += "<";
			#else
				*s += "(";
			#endif
			if (contracted_lc != NULL) {
				contracted_lc->str_c_subtree_hlpr(s);
				*s += ",";
			}
			for (auto c = contracted_children.begin(); c != contracted_children.end(); c++) {
			  (*c)->str_c_subtree_hlpr(s);
			  *s += ",";
			}

			if (contracted_rc != NULL) {
				contracted_rc->str_c_subtree_hlpr(s);
			}
			else {
			  s->pop_back();
			}
			#ifdef DEBUG_CONTRACTED
				*s += ">";
			#else
				*s += ")";
			#endif
		}
	}

	string str_subtree() {
		string s = "";
		str_subtree_hlpr(&s);
		return s;
	}

	void str_subtree_hlpr(string *s) {
		str_hlpr(s);
		if (!is_leaf()) {
			*s += "(";
			list<Node *>::iterator c;
			for(c = children.begin(); c != children.end(); c++) {
				if (c != children.begin())
					*s += ",";
				(*c)->str_subtree_hlpr(s);
				if ((*c)->parent() != this)
					cout << "#";
			}
			*s += ")";
		}
#ifdef DEBUG_DEPTHS
		*s+= ":";
	stringstream ss;
	string a;
	ss << depth;
	a = ss.str();
		*s+= a;
#endif
#ifdef DEBUG_PROTECTED
		if (edge_protected)
			*s += DEBUG_PROTECTED;
#endif
	}

	string str_support_subtree(bool allow_negative) {
		string s = "";
		str_support_subtree_hlpr(&s, allow_negative);
		return s;
	}

	string str_support_subtree() {
		return str_support_subtree(false);
	}

	void str_support_subtree_hlpr(string *s, bool allow_negative) {
		str_hlpr(s);
		if (!is_leaf()) {
			*s += "(";
			list<Node *>::iterator c;
			for(c = children.begin(); c != children.end(); c++) {
				if (c != children.begin())
					*s += ",";
				(*c)->str_support_subtree_hlpr(s, allow_negative);
				if ((*c)->parent() != this)
					cout << "#";
			}
			*s += ")";
			if (get_support() > -1 || allow_negative) {
				stringstream ss;
				ss << setprecision (2) << get_support();
				*s+= ss.str();
				if (get_support_normalization() > -1 || allow_negative) {
					stringstream ss;
					ss << "#";
					ss << get_support_normalization();
					*s+= ss.str();
				}
			}
		}
	}

	string str_edge_pre_interval_subtree() {
		string s = "";
		str_edge_pre_interval_subtree_hlpr(&s);
		return s;
	}

	void str_edge_pre_interval_subtree_hlpr(string *s) {
		str_hlpr(s);
		if (!is_leaf()) {
			*s += "(";
			list<Node *>::iterator c;
			for(c = children.begin(); c != children.end(); c++) {
				if (c != children.begin())
					*s += ",";
				(*c)->str_edge_pre_interval_subtree_hlpr(s);
				if ((*c)->parent() != this)
					cout << "#";
			}
			*s += ")";
		}
		if (get_preorder_number() > -1) {
			stringstream ss;
			ss << ":";
			if (get_edge_pre_start() > -1) {
				ss << get_edge_pre_start();
				ss << ";";
			}
			ss << get_preorder_number();
			if (get_edge_pre_end() > -1) {
				ss << ";";
				ss << get_edge_pre_end();
			}
			*s+= ss.str();
		}
	}

	void str_c_subtree_hlpr(string *s) {
		str_hlpr(s);
		if (!is_leaf()) {
			#ifdef DEBUG_CONTRACTED
				*s += "<";
			#else
				*s += "(";
			#endif
			list<Node *>::iterator c;
			for(c = children.begin(); c != children.end(); c++) {
				if (c != children.begin())
					*s += ",";
				(*c)->str_c_subtree_hlpr(s);
				if ((*c)->parent() != this)
					cout << "#";
			}
			#ifdef DEBUG_CONTRACTED
				*s += ">";
			#else
				*s += ")";
			#endif
		}
	}

	string str_subtree_twin() {
		string s = "";
		str_subtree_twin_hlpr(&s);
		return s;
	}

	void str_subtree_twin_hlpr(string *s) {
		*s += name.c_str();//str_hlpr(s);
		if (twin != NULL) {
			*s += "{";
				twin->str_subtree_hlpr(s);
			*s += "}";
		}
		if (!is_leaf()) {
			*s += "(";
			list<Node *>::iterator c;
			for(c = children.begin(); c != children.end(); c++) {
				if (c != children.begin())
					*s += ",";
				(*c)->str_subtree_hlpr(s);
				if ((*c)->parent() != this)
					cout << "#";
			}
			*s += ")";
		}
	}


	void print() {
		cout << name;
	}
	void print_subtree() {
		cout << str_subtree();
		cout << endl;
	}

  	void print_subtree_hlpr() {
		cout << str_subtree();
	}
	void print_subtree_twin_hlpr() {
		cout << str_subtree_twin();
	}

	bool is_leaf() {
		return children.empty();
	}

	// TODO: binary only
	bool is_sibling_pair() {
		return (lchild() != NULL && lchild()->is_leaf()
				&& rchild() != NULL && rchild()->is_leaf());
	}

	bool is_temp_parent(){
		return (non_leaf_children == 0 && !is_leaf() 
				&& lchild() != NULL && pre_num == lchild()->get_preorder_number());
	}

	bool is_temp_target_parent(){
		return (non_leaf_children == 0 && !is_leaf() 
				&& pre_num == get_child_min_preorder());
	}

	int get_child_min_preorder(){
		int min_pre_num = INT_MAX;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			if(min_pre_num > (*c)->get_preorder_number()){
				min_pre_num = (*c)->get_preorder_number();
			}
		}
		return min_pre_num;
	}

	void get_all_children_preorder(list<int> *children_preorders) {
	    list<Node *>::iterator c;
		if(is_leaf()){
			children_preorders->push_back(pre_num);
			return;
		}
		for (c = children.begin(); c != children.end(); c++) {
		  	if (!(*c)->is_leaf()) {
				(*c)->get_all_children_preorder(children_preorders);
		  	}
		}		
    }

    bool is_sibling_group() {
	  return non_leaf_children == 0 && !is_leaf();
    }

	bool is_singleton() {
		return (parent() == NULL && is_leaf());
	}

	// TODO: binary only
	void find_sibling_pairs_hlpr(list<Node *> *sibling_pairs) {
		Node *lchild = this->lchild();
		Node *rchild = this->rchild();
		bool lchild_leaf = false;
		bool rchild_leaf = false;
		if (lchild != NULL) {
			if (lchild->is_leaf())
				lchild_leaf = true;
			else
				lchild->find_sibling_pairs_hlpr(sibling_pairs);
		}
		if (rchild != NULL) {
			if (rchild->is_leaf())
				rchild_leaf = true;
			else
				rchild->find_sibling_pairs_hlpr(sibling_pairs);
		}
		if (lchild_leaf && rchild_leaf) {
			sibling_pairs->push_back(lchild);
			sibling_pairs->push_back(rchild);
			//lchild->add_to_sibling_pairs(sibling_pairs, 1);
			//rchild->add_to_sibling_pairs(sibling_pairs, 2);
		}
	}
	
	// find the sibling pairs in this node's subtree
	void append_sibling_pairs(list<Node *> *sibling_pairs) {
		find_sibling_pairs_hlpr(sibling_pairs);
	}

	// find the sibling pairs in this node's subtree
	list<Node *> *find_sibling_pairs() {
		list<Node *> *sibling_pairs = new list<Node *>();
		find_sibling_pairs_hlpr(sibling_pairs);
		return sibling_pairs;
	}
  
	void find_leaves_hlpr(vector<Node *> &leaves) {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			if ((*c)->is_leaf())
				leaves.push_back(*c);
			else
				(*c)->find_leaves_hlpr(leaves);
		}

	}
      // finds the parents of the sibling groups in this node's subtree
    void append_sibling_groups(list<Node *> *sibling_groups) {
	  find_sibling_groups_hlpr(sibling_groups);	
    }

    // finds the parents of the sibling groups in this node's subtree
    list<Node *> *find_sibling_groups() {
	    list<Node *> *sibling_groups = new list<Node *>();
		find_sibling_groups_hlpr(sibling_groups);
		return sibling_groups;
    }
    void find_sibling_groups_hlpr(list<Node*> *sibling_groups) {
	    list<Node *>::iterator c;
		int total_non_leaves = 0;
		for (c = children.begin(); c != children.end(); c++) {
		  if (!(*c)->is_leaf()) {
			total_non_leaves++;
		    (*c)->find_sibling_groups_hlpr(sibling_groups);
		  }
		}
		non_leaf_children = total_non_leaves;
		if (is_sibling_group())
		  sibling_groups->push_back(this);		
    }

    // Assumes this is the parent of the sibling group in T1
    //outputs a list of list of nodes of siblings in T2 that have the same parents in T2
    void find_identical_sibling_groups(list<list<Node *>> *identical_sibling_groups) {
                identical_sibling_groups->clear();
		list<Node *>::iterator c;
		map<Node *, list<Node*>> parent_to_children = map<Node *, list<Node *>>();
		for (c = children.begin(); c != children.end(); c++) {
		  Node *T1_a = *c;
		  Node *T2_a = T1_a->get_twin();
		  if (T2_a->parent() == NULL)	{
			continue;
		  }
		  auto search = parent_to_children.find(T2_a->parent());
		  //if it does not equal the end, then we already have this parent in the map
		  if (search != parent_to_children.end()) {
			search->second.push_back(T2_a);
		  }
		  //otherwise we add it
		  else {
			parent_to_children.insert({T2_a->parent(), {T2_a}});
		  }
		}
		//Check which parents have more than one identical sibling
		map<Node *, list<Node*>>::iterator parents;
		for (parents = parent_to_children.begin(); parents != parent_to_children.end(); parents++) {
		  if (parents->second.size() > 1) {
			identical_sibling_groups->push_back(parents->second);
		  }
		}
    }
  
	// find the leaves in this node's subtree
	vector<Node *> find_leaves() {
		vector<Node *> leaves = vector<Node *>();
		if (is_leaf())
			leaves.push_back(this);
		else
			find_leaves_hlpr(leaves);
		return leaves;
	}
        void find_nodes_in_subtree_hlpr(vector<Node*> *nodes) {
	  for  (auto i = children.begin(); i != children.end(); i++) {
	    (*i)->find_nodes_in_subtree_hlpr(nodes);
	  }
	  nodes->push_back(this);
	}
        vector<Node*> find_nodes_in_subtree() {
                vector<Node*> nodes = vector<Node*>();
		for (auto i = children.begin(); i != children.end(); i++) {
		  (*i)->find_nodes_in_subtree_hlpr(&nodes);
		}
		return nodes;
	}

      
  //find the preorder numbers of the leaves in this node's subtree
        vector<int> find_leaf_preorders() {
	  vector<Node*> leaves = find_leaves();
	  vector<int> preorders = vector<int>();
	  for (int i = 0; i < leaves.size(); i++) {
	    preorders.push_back(leaves[i]->get_preorder_number());
	  }
	  
	  return preorders;
	}

	void find_interior_hlpr(vector<Node *> &interior) {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			if (!(*c)->is_leaf()) {
				interior.push_back(*c);
				(*c)->find_interior_hlpr(interior);
			}
		}

	}
  
	// find the interior nodes in this node's subtree
	// does not include this node
	vector<Node *> find_interior() {
		vector<Node *> interior = vector<Node *>();
		if (!is_leaf())
			find_interior_hlpr(interior);
		return interior;
	}
  /*
        int find_arbitrary_lca_hlpr(vector<Node *> &sibling_descendants_counter) {
		list<Node *>::iterator i;
		int total = 0;
		for (i = children->begin(); i != children->end(); i++) {
		  total += (*i)->find_arbitrary_lca_hlpr(sibling_descendants_counter);
		}
		sibling_descendants_counter[get_preorder_number()] += total;
		return sibling_descendants_counter[get_preorder_number()];
	}
  
        // Finds one of the lca's in T2 based on a sibling group parent of T1
        vector<Node *> find_arbitrary_lca(Node* T1_parent) {
	        vector<Node *> sibling_descendants_counter = vector<Node *>();
		list<Node *>::iterator i;
		list<Node *> siblings = T1_parent->get_children();
		for (i = siblings->begin(); i != siblings->end(); i++)
		  {
		    Node* twin = (*i)->get_twin();
		    if (twin != NULL) {
		      //1 it
		    }
		    else {
		      //0 it
		    }
		  }
		for () {
		     //i in components of t2 forest)
		  find_arbitrary_lca_hlpr(sibling_descendants_counter);
		}

	}
  */

        int get_deepest_siblings_hlpr(vector<int> &descendants, Node** to_be_placed) {
	  list<Node*>::iterator i;
	  for (i = children.begin(); i != children.end(); i++) {
	    int descendant_value = descendants[(*i)->get_preorder_number()];
	    if (descendant_value == -1) {
	      *to_be_placed = (*i);
	      return 2;
	    }
	    else if (descendant_value == 1) {
	      return 1 + (*i)->get_deepest_siblings_hlpr(descendants, to_be_placed);
	    }
	  }
	  return -300; //shouldn't get here, it should break when it tries to index -300.
	  //should really assert
	}

        void get_deepest_siblings(vector<int> &descendants, vector<vector<Node*>> &siblings_by_depth) {
	  get_deepest_siblings_strt(descendants, siblings_by_depth);
	}
        vector<vector<Node*>> get_deepest_siblings(vector<int> &descendants) {
	  vector<vector<Node *>> siblings_by_depth = vector<vector<Node *>>(10);
	  get_deepest_siblings_strt(descendants, siblings_by_depth);
	  return siblings_by_depth;
	}
  
        /*
	  Takes a vector of ints with the number of siblings that descend 
	  from that node indexed by their pre-order number
	  Assumes siblings are marked with -1
	  TODO (Ben): proper counting sort
	 */
         void get_deepest_siblings_strt(vector<int> &descendants, vector<vector<Node*>> &siblings_by_depth) {
	   if (get_preorder_number() == -1) { return; } //if -1 then this is rho
	   //This is the sibling, and is a component
	   if (descendants[get_preorder_number()] == -1) {
	     siblings_by_depth[0].push_back(this);
	     return;
	   }
	   list<Node *>::iterator i;
	  //Go down each path of 1s and take note of their path length to this node
	  //There should be the same amount of 1 paths as nodes.size()
	  for (i = children.begin(); i != children.end(); i++) {

	    if (descendants[(*i)->get_preorder_number()] == 1) {
	      Node* placed_node;
	      int depth = (*i)->get_deepest_siblings_hlpr(descendants, &placed_node);
	      while (siblings_by_depth.size() - 1 < depth)
		{
		  siblings_by_depth.resize(siblings_by_depth.size() * 2);
		}
	      if (placed_node != NULL){
		siblings_by_depth[depth].push_back(placed_node);
	      }
	    }	
	    //if it is -1 then this child is the node itself
	    else if (descendants[(*i)->get_preorder_number()] == -1) {
	      siblings_by_depth[1].push_back(*i);
	    }
	  }
	}

        
  
        Node *find_arbitrary_lca_hlpr(vector<int> &descendants) {
	        if (descendants[get_preorder_number()] < 2) { return NULL; }
		list<Node*>::iterator i;
	        for (i = children.begin(); i != children.end(); i++) {
		  if (descendants[(*i)->get_preorder_number()] > 1) {
		    return (*i)->find_arbitrary_lca_hlpr(descendants);
		  }
		}
		return this;
	}
        /*
	  Takes a vector of components in T2, 
	  and a count of how many sibling descendants the node 
	  has in the index of its preorder number
	  Returns lca of the nodes
	  Returns NULL if there is no path between any of the two siblings
	  
	  NOTE: this can be a function
        */
        Node *find_arbitrary_lca(vector<Node *> components, vector<int> &descendants) {	        
		for (int i = 0; i < components.size(); i++) {
		  Node* root = components[i];
		  if (root->get_preorder_number() == -1) {continue;}
		  if (descendants[root->get_preorder_number()] > 1) {
		    return root->find_arbitrary_lca_hlpr(descendants);
		  }
		}
		return NULL;
	}
  /*This is potentially slower than just doing a traversal of the whole tree, 
    since we could be visiting the same path many times. 
    It does mean if this is a short branch and the others are
    long that we dont visit those branches though
    TODO: traverse tree and sum on return. Label sibling group with 1s rest with 0
  */
        void find_pseudo_lca_descendant_count_hlpr(vector <int> &desc_counter) {
		if (p != NULL) {
		  desc_counter[p->get_preorder_number()]++;
		  p->find_pseudo_lca_descendant_count_hlpr(desc_counter);
		}
	}
        vector<int> find_pseudo_lca_descendant_count(int max_preorder) {
	        vector<int> sibling_descendants_counter = vector<int>(max_preorder, 0);
		list<Node *>::iterator i;
		for (i = children.begin(); i != children.end(); i++)
		  {
		    Node* twin = (*i)->get_twin();
		    if (twin != NULL) {
		      //set the twin itself to -1
		      //This will be used later when traversing down to notify we have hit the node
		      sibling_descendants_counter[twin->get_preorder_number()] = -1;
		      twin->find_pseudo_lca_descendant_count_hlpr(sibling_descendants_counter);
		    }
		  }
		return sibling_descendants_counter;
	}
  
	void find_descendants_hlpr(vector<Node *> &descendants) {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
				descendants.push_back(*c);
			if (!(*c)->is_leaf()) {
				(*c)->find_descendants_hlpr(descendants);
			}
		}

	}
	
	// find the descendants nodes in this node's subtree
	// does not include this node
	vector<Node *> find_descendants() {
		vector<Node *> descendants = vector<Node *>();
		if (!is_leaf())
			find_descendants_hlpr(descendants);
		return descendants;
	}

	bool contains_leaf(int number) {
		if (stomini(get_name()) == number)
			return true;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			bool ans = (*c)->contains_leaf(number);
			if (ans == true)
				return true;
		}
		vector<Node *> descendants = vector<Node *>();
		return false;
	}

	// make twins point to this tree in this node's subtree
	void resync() {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->resync();
		}
		if (twin != NULL)
			twin->set_twin(this);
	}

	// remove all twins
	void unsync() {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->unsync();
		}
		twin = NULL;
	}

	// remove all interior twins
	void unsync_interior() {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->unsync();
		}
		if(!is_leaf())
			twin = NULL;
	}

	void sync_af_twins() {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->sync_af_twins();
		}
		if (!is_leaf()) {
			twin = lchild()->get_twin()->parent();
			twin->set_twin(this);
		}
	}

	// find the root of this node's tree
	Node *find_root() {
		Node *root = this;
		while (root->parent() != NULL)
			root = root->parent();
		return root;
	}

	bool same_component(Node *n, int &lca_depth, int &path_length) {
		Node *a = this;
		Node *b = n;
		path_length = 0;
		while(a != b) {
			if ((b == NULL) || (a != NULL && a->get_depth() > b->get_depth()))
				a = a->parent();
			else
				b = b->parent();
			path_length++;
		}
		if (a == NULL) {
			path_length = -1;
			return false;
		}
		lca_depth = a->get_depth();
		return true;
	}

	bool same_component(Node *n) {
		int a,b;
		return same_component(n, a, b);
	}

	bool same_component(Node *n, int &lca_depth) {
		int a;
		return same_component(n, lca_depth, a);
	}

 
	void labels_to_numbers(map<string, int> *label_map, map<int, string> *reverse_label_map) {
		if (name != "") {
			map<string, int>::iterator i = label_map->find(name);
			stringstream ss;
			if (i != label_map->end()) {
				ss << i->second;
				name = ss.str();
			}
			else {
				int num = label_map->size();
				ss << num;
				label_map->insert(make_pair(name, num));
				reverse_label_map->insert(make_pair(num, name));
				name = ss.str();
			}
		}
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->labels_to_numbers(label_map, reverse_label_map);
		}
		if (contracted_lc != NULL)
			contracted_lc->labels_to_numbers(label_map, reverse_label_map);
		if (contracted_rc != NULL)
			contracted_rc->labels_to_numbers(label_map, reverse_label_map);
		for(c = contracted_children.begin(); c != contracted_children.end(); c++) {
			(*c)->labels_to_numbers(label_map, reverse_label_map);
		}
	}


  	void numbers_to_labels(map<int, string> *reverse_label_map) {
		if (name != "") {
			string converted_name = "";
			string current_num = "";
			string::iterator i = name.begin();
			size_t old_loc = 0;
			size_t loc = 0;
			while ((loc = name.find_first_of("0123456789", old_loc)) != string::npos) {
				converted_name.append(name.substr(old_loc, loc - old_loc)); 
				old_loc = loc;
				loc = name.find_first_not_of("0123456789", old_loc);
				string label = "";
				if (loc == string::npos)
					loc = name.size();
				label = name.substr(old_loc, loc - old_loc);
				map<int, string>::iterator j = reverse_label_map->find(atoi(label.c_str()));
				if (j != reverse_label_map->end()) {
					stringstream ss;
					ss << j->second;
					label = ss.str();
				}
				converted_name.append(label);
				old_loc = loc;
			}
			converted_name.append(name.substr(old_loc, name.size() - old_loc)); 
			name = converted_name;



		}
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->numbers_to_labels(reverse_label_map);
		}
		if (contracted_lc != NULL)
			contracted_lc->numbers_to_labels(reverse_label_map);
		if (contracted_rc != NULL)
			contracted_rc->numbers_to_labels(reverse_label_map);
		for(c = contracted_children.begin(); c != contracted_children.end(); c++) {
			(*c)->numbers_to_labels(reverse_label_map);
		}

	}

	void build_name_to_pre_map(map<string, int> *name_to_pre) {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->build_name_to_pre_map(name_to_pre);
		}
		if (is_leaf())
			name_to_pre->insert(make_pair(get_name(), get_preorder_number()));
	}

	void count_numbered_labels(vector<int> *label_counts) {
		if (name != "") {
			int label = stomini(name);
			if (label_counts->size() <= label)
				label_counts->resize(label+1,0);
			(*label_counts)[label]++;
		}

		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->count_numbered_labels(label_counts);
		}
		if (contracted_lc != NULL)
			contracted_lc->count_numbered_labels(label_counts);
		if (contracted_rc != NULL)
			contracted_rc->count_numbered_labels(label_counts);
	}

  int get_max_preorder_number(int next) {
    if (next < get_preorder_number()) {
      next = get_preorder_number();
    }
    for (auto c = children.begin(); c != children.end(); c++) {
      next = (*c)->get_max_preorder_number(next);
    }
    return next;
  }
  
	void print_preorder_number() {
		cout << name << " - " << pre_num << " - " << edge_pre_start << " - " << edge_pre_end << endl;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->print_preorder_number();
		}
	}

	void preorder_number() {
		preorder_number(0);
	}
  

	int preorder_number(int next) {
		set_preorder_number(next);
		next++;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			next = (*c)->preorder_number(next);
		}
		return next;
	}

	void edge_preorder_interval() {
		edge_pre_start = pre_num;
		if (is_leaf()) {
			edge_pre_end = pre_num;
		}
		else {
			list<Node *>::iterator c;
			for(c = children.begin(); c != children.end(); c++) {
				(*c)->edge_preorder_interval();
				if (edge_pre_end == -1 || (*c)->edge_pre_end > edge_pre_end)
					edge_pre_end = (*c)->edge_pre_end;
			}
		}
	}

	int size() {
		int s  = 1;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			s += (*c)->size();
		}
		return s;
	}

	int size_using_prenum() {
		if (is_leaf())
			return get_preorder_number();
		else
			return children.back()->size_using_prenum();
	}

	int max_depth() {
		int d  = 0;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			int c_d = (*c)->max_depth();
			if (c_d > d)
				d = c_d;
		}
		return d+1;
	}

	int max_degree() {
		int d = children.size();
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			int c_d = (*c)->max_degree();
			if (c_d > d)
				d = c_d;
		}
		return d;
	}

	// TODO: binary only
	// these will potentially be removed

void add_to_front_sibling_pairs(list<Node *> *sibling_pairs, int status) {
	sibling_pairs->push_front(this);
	clear_sibling_pair(sibling_pairs);
	sibling_pair_status = status;
	sibling_pair_loc = sibling_pairs->begin();
}

void add_to_sibling_pairs(list<Node *> *sibling_pairs, int status) {
	sibling_pairs->push_back(this);
	clear_sibling_pair(sibling_pairs);
	sibling_pair_status = status;
	sibling_pair_loc = sibling_pairs->end();
	sibling_pair_loc--;
}

void remove_sibling_pair(list<Node *> *sibling_pairs) {
	if (sibling_pair_status > 0) {
		list<Node *>::iterator loc = sibling_pair_loc;
		list<Node *>::iterator sibling_loc = loc;
		if (sibling_pair_status == 1)
			sibling_loc++;
		else if (sibling_pair_status == 2)
			sibling_loc--;

		if (sibling_loc != sibling_pairs->end()) {
			Node *old_sibling = *sibling_loc;
			old_sibling->sibling_pair_status = 0;
			sibling_pairs->erase(sibling_loc);
		}
		sibling_pairs->erase(loc);
		sibling_pair_status = 0;
	}
}

void clear_sibling_pair(list<Node *> *sibling_pairs) {
	if (sibling_pair_status > 0) {
		list<Node *>::iterator loc = sibling_pair_loc;
		list<Node *>::iterator sibling_loc = loc;
		if (sibling_pair_status == 1)
			sibling_loc++;
		else if (sibling_pair_status == 2)
			sibling_loc--;

		if (sibling_loc != sibling_pairs->end()) {
			Node *old_sibling = *sibling_loc;
			old_sibling->sibling_pair_status = 0;
		}
		sibling_pair_status = 0;
	}
}

Node *get_sibling() {
	Node *ret = NULL;
	if (p != NULL && p->children.size() > 1) {
		list<Node *>::iterator s = p_link;
		if (p_link == p->children.begin())
			s++;
		else
			s--;
		ret = *s;
	}
	return ret;
}

Node *get_right_sibling() {
	Node *ret = NULL;
	if (p != NULL && p->children.size() > 1) {
		list<Node *>::iterator s = p_link;
		s++;
		if (s != children.end())
			ret = *s;
	}
	return ret;
}

// TODO: binary only
Node *get_sibling(list<Node *> *sibling_pairs) {
	if (sibling_pair_status > 0) {
		list<Node *>::iterator loc = sibling_pair_loc;
		list<Node *>::iterator sibling_loc = loc;
		if (sibling_pair_status == 1)
			sibling_loc++;
		else if (sibling_pair_status == 2)
			sibling_loc--;
		return *sibling_loc;
	}
	else
		return NULL;
}

void set_sibling(Node *sibling) {
	if (sibling->sibling_pair_status > 0) {
		sibling_pair_loc = sibling->sibling_pair_loc;
		if (sibling->sibling_pair_status == 1)
			sibling_pair_loc++;
		else if (sibling->sibling_pair_status == 2)
			sibling_pair_loc--;
	}
}

void clear_sibling_pair_status() {
	sibling_pair_status = 0;
}

// fix parents
void fix_parents() {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			if ((*c)->parent() != this) {
				(*c)->p = this;
				(*c)->p_link = c;
			}
			(*c)->fix_parents();
		}
}

// TODO: binary only
/* TODO: think of a better way of trying all rootings. Maybe create
	 a function that copies a tree and then reroots it with a given
	 subtree as a child of the root.
*/

void left_rotate() {
	if (lchild()->lchild() != NULL) {
		Node *new_lc = lchild()->lchild();
		Node *new_rc = lchild();
		Node *new_rc_lc = lchild()->rchild();
		Node *new_rc_rc = rchild();

		new_lc->cut_parent();
		new_rc->cut_parent();
		new_rc_lc->cut_parent();
		new_rc_rc->cut_parent();
		add_child(new_lc);
		add_child(new_rc);
		rchild()->add_child(new_rc_lc);
		rchild()->add_child(new_rc_rc);
	}
}

void right_rotate() {
	if (rchild()->lchild() != NULL) {
		Node *new_lc = rchild()->lchild();
		Node *new_rc = rchild();
		Node *new_rc_lc = rchild()->rchild();
		Node *new_rc_rc = lchild();

		new_lc->cut_parent();
		new_rc->cut_parent();
		new_rc_lc->cut_parent();
		new_rc_rc->cut_parent();
		add_child(new_lc);
		add_child(new_rc);
		rchild()->add_child(new_rc_lc);
		rchild()->add_child(new_rc_rc);
	}
}

void next_rooting() {
	if (lchild()->lchild() != NULL)
		left_rotate();
	else if (rchild()->lchild() != NULL)
			right_rotate();
	else
		return;
	if (lchild()->pre_num < rchild()->pre_num && (lchild()->pre_num != 1 || rchild()->pre_num != 2 || !lchild()->is_leaf()))
		next_rooting();
}

/* reroot
 * rerooots the tree between new_lc and its parent. Maintains this
 * Node as the root
 * assumes that new_lc is a descendant of this Node
 * assumes this node is the root and has two children
 * other nodes may be multifurcating
*/
void reroot(Node *new_lc) {
	Node *new_rc = new_lc->parent();
	if (new_rc == this || new_rc == NULL) {
		return;
	}
	Node *prev = new_rc;
	Node *next = new_rc->parent();
//	Node *old_lc = lchild();
//	Node *old_rc_rc = rchild();
	new_lc->cut_parent();
	new_rc->cut_parent();
	while(next != NULL) {
		Node *current = next;
		next = current->parent();
		prev->add_child(current);
		prev = current;
	}
	Node *root = prev;
	while(root->get_children().size() > 0)
		root->parent()->add_child(root->lchild());
//	root->parent()->add_child(root->lchild());
	root->cut_parent();
	root->add_child(new_lc);
	root->add_child(new_rc);
}

void reroot_clean(Node *new_lc) {
	this->reroot(new_lc);
	this->set_depth(0);
	this->fix_depths();
	this->preorder_number();
	this->edge_preorder_interval();
}

// make the root binay again
void fixroot() {
	if (get_children().size() > 2) {
		Node *new_lc = new Node();
		while(get_children().size() > 1) {
			new_lc->add_child(get_children().back());
		}
		add_child(new_lc);
	}
}

  //Currently only used on T2, if you want to use it on T1, update non_leaf_children
Node *expand_children_out(list<Node*> nodes) {
  Node *new_child = new Node();
  //new_child->set_parent(this);
  int min_preorder = 0x7FFFFFFF;
  list<Node *>::iterator i;
  for (i = nodes.begin(); i != nodes.end(); i++) {
    new_child->add_child(*i);
    (*i)->set_parent(new_child);//not needed?
    children.remove(*i);
    if ((*i)->get_preorder_number() < min_preorder) {
      min_preorder = (*i)->get_preorder_number();
    }
  }
  //super hacky, may overlap with others. Currently this is being used only when the other node is being contracted anyways, but be careful
  new_child->set_preorder_number(min_preorder);
  add_child(new_child);
  return new_child;
}

Node *expand_parent_edge(Node *n) {
	if (p != NULL) {
		Node *old_p = p;
		cut_parent();
		Node *new_p = new Node();
		new_p->add_child(n);
		old_p->add_child(new_p);
		return p;
	}
	else {
		Node *new_child = new Node(name);
		Node *old_lc = lchild();
		Node *old_rc = rchild();
		old_lc->cut_parent();
		old_rc->cut_parent();

		new_child->add_child(old_lc);
		new_child->add_child(old_rc);
		new_child->contracted_lc = contracted_lc;
		if (contracted_lc != NULL)
			contracted_lc->p = new_child;
		new_child->contracted_rc = contracted_rc;
		if (contracted_rc != NULL)
		contracted_rc->p = new_child;
		name = "";
		contracted_lc = NULL;
		contracted_rc = NULL;
		add_child(new_child);
		return this;
	}
}

Node *undo_expand_parent_edge() {
	if (p != NULL) {
		Node *old_p = p;
		cut_parent();
		Node *child = children.front();
		child->cut_parent();
		old_p->add_child(child);
		return this;
	}
	else {
		Node *child = children.front();
		name = child->name;
		Node *new_lc = child->lchild();
		Node *new_rc = child->rchild();
		new_lc->cut_parent();
		new_rc->cut_parent();
		contracted_lc = child->contracted_lc;
		if (contracted_lc != NULL)
			contracted_lc->p = this;
		contracted_rc = child->contracted_rc;
		if (contracted_rc != NULL)
			contracted_rc->p = this;
		child->contracted_lc = NULL;
		child->contracted_rc = NULL;
		child->name = "";
		add_child(new_lc);
		add_child(new_rc);
		return child;
	}
}

/*  apply an SPR operation to move this subtree to be a
 *	sibling of new_sibling
 *
 *  Note: problems will occur if new_sibling is a descendant of this
 *  Note: does nothing if this is the root
 *
 *	Returns a node that will reverse the spr (old_sibling unless it
 *		was moved to maintain a root, in which case the root is returned,
 *		NULL if no spr occured)
 *
 * Note: binary only
 * TODO (ben): multifurcating
 */
Node *spr(Node *new_sibling, int &which_child) {
	Node *reverse;
	int prev_child_loc = 0;
	if (p == NULL || new_sibling == NULL)
		return NULL;
	Node *old_sibling = get_sibling();
	if (old_sibling == new_sibling)
		return NULL;
	Node *grandparent = p->p;
	if (p->lchild() == this)
		prev_child_loc = 1;
	else
		prev_child_loc = 2;
	// Prune
	if (grandparent != NULL) {
		bool leftc = false;
		if (old_sibling->parent() == grandparent->lchild())
			leftc = true;

		old_sibling->cut_parent();
		//p->delete_child(old_sibling);
		p->cut_parent(); //LEAK?
//		grandparent->delete_child(p);
		if (leftc && !grandparent->is_leaf()) {
				Node *ns = grandparent->children.front();
				grandparent->insert_child(ns,old_sibling);
		}
		else
			grandparent->add_child(old_sibling);

		reverse = old_sibling;
	}
	else {
		if (old_sibling->is_leaf())
			return NULL;
		Node *root = p;
		bool leftc = false;
		if (root->lchild() == this)
			leftc = true;
		Node *lc = old_sibling->lchild();
		Node *rc = old_sibling->rchild();
		root->delete_child(this);
		root->delete_child(old_sibling) ;
		if (lc != NULL)
			root->add_child(lc);
		if (rc != NULL)
			root->add_child(rc);

		//old_sibling->delete_child(old_sibling->lchild());
		//old_sibling->delete_child(old_sibling->rchild());
		if (leftc) {
			if (old_sibling->is_leaf())
				old_sibling->add_child(this);
			else {
				Node *new_sibling = old_sibling->children.front();
				old_sibling->insert_child(new_sibling,this);
			}
		}
		else
			old_sibling->add_child(this);
		reverse = root;
	}


	// Regraft
	if (new_sibling->p != NULL) {
		grandparent = new_sibling->p;
//		grandparent->delete_child(new_sibling);
		bool leftc = false;
		if (new_sibling == grandparent->lchild())
			leftc = true;
		if (leftc && !grandparent->is_leaf()) {
				Node *new_sibling = grandparent->children.front();
				grandparent->insert_child(new_sibling,p);
		}
		else {
			grandparent->add_child(p);
		}

		if (which_child == 1)
			p->add_child(new_sibling);
		else {
				Node *ns = p->children.front();
				p->insert_child(ns,new_sibling);
		}

	}
	else {
		Node *root = new_sibling;
		new_sibling = p;
		// TODO: still broken
		Node *lc = root->lchild();
		Node *rc = root->rchild();
		p->add_child(lc);
		p->add_child(rc);
//		p->delete_child(this);
		// problem here
		if (which_child == 0)
			which_child = prev_child_loc;
		if (which_child == 1) {
			root->add_child(this);
			root->add_child(new_sibling);
		}
		else {
			root->add_child(new_sibling);
			root->add_child(this);
		}


	}

	which_child = prev_child_loc;
	return reverse;
}

Node *spr(Node *new_sibling) {
	int na = 0;
	return spr(new_sibling, na);
}

  /*
    prune this subtree and attach as a sibling of new_sibiling
   */

void spr_mult(Node *new_sibling) {
  //prune
  Node* parent = p;
  if (parent != NULL) {
    cut_parent();
    if (parent->get_children().size() == 1) {
      parent->contract(true);
    }
  }
  else {
    return;
  }
  //regraft
  if (new_sibling->parent() != NULL) {
    Node* new_parent = new Node();
    new_sibling->parent()->add_child(new_parent);
    //new_parent->set_parent(new_sibling->parent());
    new_parent->add_child(new_sibling);
    new_parent->add_child(this);
    }
  else {
    Node* new_parent = new_sibling;//new Node();
    Node* replace = new Node(*new_sibling);
    new_parent->get_children().clear();
    new_parent->add_child(replace);
    new_parent->add_child(this);
  }
}

  
void find_descendant_counts_hlpr(vector<int> *dc) {
	int num_descendants = 0;
	list<Node *>::iterator c;
	for(c = children.begin(); c != children.end(); c++) {
		(*c)->find_descendant_counts_hlpr(dc);
		num_descendants += (*dc)[(*c)->get_preorder_number()];
		num_descendants += 1;
	}
	if (dc->size() <= get_preorder_number())
		dc->resize(get_preorder_number() + 1, -1);
	(*dc)[get_preorder_number()] = num_descendants;
}

vector<int> *find_descendant_counts() {
	vector<int> *dc = new vector<int>();
	find_descendant_counts_hlpr(dc);
	return dc;
}

void find_leaf_counts_hlpr(vector<int> *lc) {
	int num_leaves = 0;
	list<Node *>::iterator c;
	for(c = children.begin(); c != children.end(); c++) {
		(*c)->find_leaf_counts_hlpr(lc);
		num_leaves += (*lc)[(*c)->get_preorder_number()];
	}
	if (is_leaf())
		num_leaves++;
	if (lc->size() <= get_preorder_number())
		lc->resize(get_preorder_number() + 1, -1);
	(*lc)[get_preorder_number()] = num_leaves;
}

vector<int> *find_leaf_counts() {
	vector<int> *lc = new vector<int>();
	find_leaf_counts_hlpr(lc);
	return lc;
}

Node *find_median() {
	vector<int> *dc = find_descendant_counts();
	return find_median_hlpr(dc, (*dc)[get_preorder_number()] / 2);
}

Node *find_subtree_of_size(double percentage) {
	vector<int> *dc = find_descendant_counts();
	return find_median_hlpr(dc, (int)((*dc)[get_preorder_number()] * percentage));
}

Node *find_median_hlpr(vector<int> *dc, int target_size) {
	Node *largest_child_subtree = NULL;
	int lcs_size = 0;
	list<Node *>::iterator c;
	for(c = children.begin(); c != children.end(); c++) {
		int cs_size = (*dc)[(*c)->get_preorder_number()] + 1; 
		if (cs_size > lcs_size) {
			largest_child_subtree = *c;
			lcs_size = cs_size;
		}
	}
	if (lcs_size > target_size)
		return largest_child_subtree->find_median_hlpr(dc, target_size);
	else
		return this;
}

int any_leaf_preorder_number() {
	if (is_leaf()) {
		if (contracted_lc != NULL)
			return contracted_lc->any_leaf_preorder_number();
		else return get_preorder_number();
	}
	else
		return (*children.begin())->any_leaf_preorder_number();
}

Node *find_by_prenum(int prenum) {
//	cout << "find_by_prenum: " << str_subtree() << endl;
	if (prenum == get_preorder_number())
		return this;
	Node *search_child = NULL;
	int best_prenum = -1;
	list<Node *>::iterator c;
	for(c = children.begin(); c != children.end(); c++) {
//		cout << "\tchild: " << (*c)->str_subtree() << endl;
//		cout << "\tpre: " << (*c)->get_preorder_number() << endl;
		int p = (*c)->get_preorder_number();
		if (p == prenum)
			return *c;
		else if (p < prenum && p > best_prenum) {
			best_prenum = p;
			search_child = *c;
		}
	}
	if (search_child != NULL)
		return search_child->find_by_prenum(prenum);
	else
		return NULL;
}

Node *find_by_prenum_full(int prenum) {
	if (prenum == get_preorder_number())
		return this;
	list<Node *>::iterator c;
	for(c = children.begin(); c != children.end(); c++) {
		Node *test = (*c)->find_by_prenum_full(prenum);
		if (test != NULL)
			return test;
	}
	return NULL;
}

// expand all contracted nodes of a subtree starting at n
void expand_contracted_nodes() {
  	list<Node *>::iterator c;
	if (is_leaf()) {
		if (contracted_lc != NULL) {
			add_child(contracted_lc);
			contracted_lc = NULL;
		}
		if (contracted_rc != NULL) {
			add_child(contracted_rc);
			contracted_rc = NULL;
		}
		for(c = contracted_children.begin(); c != contracted_children.end(); c++) {
		        add_child(*c);
		}
		contracted_children.clear();

	}
	for(c = get_children().begin(); c != get_children().end(); c++) {
		(*c)->expand_contracted_nodes();
	}
	/*
	//if -1 then this is a node we created for a contraction, and it should go up one
	if (get_edge_pre_end() == -1) {
	  if (p != NULL) {
	    while (!children.empty()) {
	      p->add_child(children.front());
	    }
	    delete this;
	    return;
	  }
	  }*/
}

int get_name_num() {
	return atoi(get_name().c_str());
}


};

// function prototypes
Node *build_tree(string s);
Node *build_tree(string s, set<string, StringCompare> *include_only);
Node *build_tree(string s, int start_depth);
Node *build_tree(string s, int start_depth, set<string, StringCompare> *include_only);

int build_tree_helper(int start, const string& s, Node *parent,
		bool &valid, set<string, StringCompare> *include_only);
//void preorder_number(Node *node);
//int preorder_number(Node *node, int next);
string strip_newick_name(string &T);
vector<Node *> contract_deepest_siblings(vector<vector<Node *>> &siblings_by_depth);


// build a tree from a newick string
Node *build_tree(string s) {
	return build_tree(s, 0, NULL);
}
Node *build_tree(string s, int start_depth) {
	return build_tree(s, start_depth, NULL);
}
Node *build_tree(string s, set<string, StringCompare> *include_only) {
	return build_tree(s, 0, include_only);
}
Node *build_tree(string s, int start_depth, set<string, StringCompare> *include_only) {
	if (s == "")
		return new Node();
	Node *dummy_head = new Node("p", start_depth-1);
	bool valid = true;
	build_tree_helper(0, s, dummy_head, valid, include_only);
	Node *head = dummy_head->lchild();
	if (valid && head != NULL) {
		delete dummy_head;
		return head;
	}
	else {
		if (head != NULL)
			head->delete_tree();
		return dummy_head;
	}

}

// build_tree recursive helper function
int build_tree_helper(int start, const string& s, Node *parent,
		bool &valid, set<string, StringCompare> *include_only) {
	int loc = s.find_first_of("(,)", start);
	if (loc == string::npos) {
		string name = s.substr(start, s.size() - start);
		int name_end = name.find(':');
		if (name_end != string::npos)
			name = name.substr(0, name_end);
		if (include_only == NULL ||
				include_only->find(name) != include_only->end()) {
			Node *node = new Node(name);
			parent->add_child(node);
		}
		loc = s.size()-1;
		return loc;
	}
	while(s[start] == ' ' || s[start] == '\t')
		start++;
	int end = loc;
	while(s[end] == ' ' || s[end] == '\t')
		end--;
	string name = s.substr(start, end - start);
	int name_end = name.find(':');
	if (name_end != string::npos)
		name = name.substr(0, name_end);
	Node *node = NULL;
	if (include_only == NULL ||
			include_only->find(name) != include_only->end()) {
		node = new Node(name);
		parent->add_child(node);
	}

	int count = 1;
	if (s[loc] == '(') {
			loc = build_tree_helper(loc + 1, s, node, valid, include_only);
			while(s[loc] == ',') {
				loc = build_tree_helper(loc + 1, s, node, valid, include_only);
				count++;
			}
//			int loc_check = s.find_first_of("(,)", loc);
//			if (loc_check != string::npos &&
//					s[loc_check] == ','
			if (s[loc] != ')'
					|| (IGNORE_MULTI && count > 2)) {
				valid = false;
					return s.size()-1;
			}
			// TODO: get the support values here (and branch lengths?)
			// contract_node() if support is less than a threshold
			loc++;
			if (s[loc-1] == ')') {
				int numc = node->get_children().size();
				bool contracted = false;
				int next = s.find_first_of(",)", loc);
				if (next != string::npos) {
					if (next > loc && (REQUIRED_SUPPORT > 0 || MIN_LENGTH >= 0)) {
						string info = s.substr(loc, next - loc);
						// support values (a number after a close bracket)
						if (info[0] != ':') {
							double support = atof(info.c_str());
							#if DEBUG_SUPPORT
							cout << "support=" << support << endl;
							#endif
							if (REQUIRED_SUPPORT > 0 && support < REQUIRED_SUPPORT && numc > 0) {
							  #if DEBUG_SUPPORT
 								cout << "contracting (support)" << endl;
								#endif
								node->contract_node();
								contracted = true;
							}
						}
						// branch lengths (a number after a colon)
						int next_colon = info.find_first_of(":");
						if (next_colon != string::npos) {
							string length_str = info.substr(next_colon+1);
							double length = atof(length_str.c_str());
#if DEBUG_SUPPORT
							cout << "length=" << length << endl;
							#endif
							if (MIN_LENGTH >= 0 && !contracted && length <= MIN_LENGTH && numc > 0) {
#if DEBUG_SUPPORT
								cout << "contracting (length)" << endl;
								#endif
								node->contract_node();
								contracted = true;
							}
						}
						// support values (a number in square brackets)
						int square_bracket = info.find_first_of("[");
						if (square_bracket != string::npos) {
							string support_str = info.substr(square_bracket+1);
							double support = atof(support_str.c_str()) / 100.0;
							#if DEBUG_SUPPORT
							cout << "support=" << support << endl;
							#endif
							if (REQUIRED_SUPPORT > 0 && !contracted && support < REQUIRED_SUPPORT && numc > 0) {
							  #if DEBUG_SUPPORT
 								cout << "contracting (support)" << endl;
								#endif
								node->contract_node();
								contracted = true;
							}
						}

					}
					loc=next;
				}
				if (!contracted) {
					if (numc == 1)
						node->contract_node();
					else if (numc == 0 && name == "") {
						node->cut_parent();
						delete node;
					}
				}
			}
	}
	return loc;
}

// swap two nodes
void swap(Node **a, Node **b) {
	Node *temp = *a;
	*a = *b;
	*b = temp;
}


/*
void preorder_number(Node *node) {
	preorder_number(node, 0);
}
int preorder_number(Node *node, int next) {
	node->set_preorder_number(next);
	next++;
	if(node->lchild() != NULL) {
		next = preorder_number(node->lchild(), next);
	}
	if(node->rchild() != NULL) {
		next = preorder_number(node->rchild(), next);
	}
	return next;
}
*/
// return the smallest number in s

int stomini(string s) {
//	cout << "stomini" << endl;
	string number_characters = "+-0123456789";
	int i = 0;
	int min = INT_MAX;
	string current = "";
	for(int i = 0; i < s.size(); i++) {
		if (number_characters.find(s[i]) != string::npos) {
			current += s[i];
		}
		else if (current.size() > 0) {
			int num = atoi(current.c_str());
			if (num < min)
				min = num;
			current = "";
		}
	}
	if (current.size() > 0) {
		int num = atoi(current.c_str());
		if (num < min)
			min = num;
		current = "";
	}
//	cout << "returning " << min << endl;
	return min;
}

// assumes that an unrooted tree is represented with a 3-way multifurcation
string root(string s) {
//	cout << "root(string s)" << endl;
//	cout << s << endl;
	string r = "";
	int i = 0;
	int depth = 0;
	int first_c = -1;
	int second_c = -1;
	int last_bracket = -1;
	for(int i = 0; i < s.size(); i++) {
		if (s[i] == '(')
			depth++;
		else if (s[i] == ')') {
			depth--;
			last_bracket = i;
		}
		else if (depth == 1 && s[i] == ',') {
			if (first_c == -1)
				first_c = i;
			else if  (second_c == -1)
				second_c = i;
		}
	}
	if (second_c == -1 || last_bracket == -1)
		return s;
	else {
		r.append(s.substr(0,first_c+1));
		r.append("(");
		r.append(s.substr(first_c+1,last_bracket-first_c));
		r.append(")");
		r.append(s.substr(last_bracket+1,string::npos));
	}
//	cout << r << endl;
	return r;
}

template <typename T> vector<T> &random_select(vector <T> &V, int n) {
	vector<T> *ret = new vector<T>;
	int end = V.size();
	for(int i = 0; i < n && end > 0; i++) {
		int x = rand() % end;
		ret->push_back(V[x]);
		V[x] = V[end-1];
		// Is it better to remove them or leave them?
		V.pop_back();
		end--;
	}
	return *ret;
}

template <typename T> void print_vector(vector <T> V) {
	int end = V.size();
	for(int i = 0; i < end; i++) {
		if (i > 0)
			cout << ",";
		cout << V[i];
	}
	cout << endl;
}

string strip_newick_name(string &line) {
	string name;
	size_t loc = line.find_first_of("(");
	if (loc != string::npos) {
		name = "";
		if (loc != 0) {
			name = line.substr(0,loc);
			line.erase(0,loc);
		}
	}
	return name;
}

void reroot_safe(Node **n, Node *new_lc) {
	(*n)->reroot(new_lc);
	(*n)->set_depth(0);
	(*n)->fix_depths();
	(*n)->preorder_number();
	(*n)->edge_preorder_interval();
	// hack TODO fix
	string str = (*n)->str_subtree();
	(*n)->delete_tree();
	(*n) = build_tree(str);
	(*n)->preorder_number();
	(*n)->edge_preorder_interval();
}

//TODO: remove references to this, replace with map version
vector<Node *> contract_deepest_siblings(vector<vector<Node *>> &siblings_by_depth) {	  
	  vector<Node*> deepest_siblings = vector<Node*>();
	  for (int j = siblings_by_depth.size() - 1; j >= 0; j--) {
	    if (siblings_by_depth[j].size() != 0) {
	      for (int k = 0; k < siblings_by_depth[j].size(); k++) {
		deepest_siblings.push_back(siblings_by_depth[j][k]);
	      }
	    }
	  }
	  return deepest_siblings; 
	}
vector<Node *> contract_deepest_siblings(vector<vector<Node *>> &siblings_by_depth, map<Node*, int> *s_map) {	  
	  vector<Node*> deepest_siblings = vector<Node*>();
	  for (int j = siblings_by_depth.size() - 1; j >= 0; j--) {
	    if (siblings_by_depth[j].size() != 0) {
	      for (int k = 0; k < siblings_by_depth[j].size(); k++) {
		deepest_siblings.push_back(siblings_by_depth[j][k]);
		s_map->insert({siblings_by_depth[j][k], j - 1});
	      }
	    }
	  }
	  
	  return deepest_siblings; 
	}





#endif

// ===== END heuristic/clean_refactor/h36_whidden/Node.h =====

// ===== BEGIN heuristic/clean_refactor/h36_whidden/LCA.h =====
/*******************************************************************************
LCA.h

Data structure for LCA computations on a binary tree
Implementation of the RMQ-based methods of Bender and Farach-Colton

Copyright 2010-2014 Chris Whidden
cwhidden@dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
March 3, 2014
Version 1.2.1

This file is part of rspr.

rspr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

rspr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rspr.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#ifndef INCLUDE_LCA

#define INCLUDE_LCA
#include <cstdio>
#include <string>
#include <iostream>
// skipped already inlined: heuristic/clean_refactor/h36_whidden/Node.h
#include <vector>
#include <cmath>
using namespace std;

int mylog2 (int val) {
	if (val <= 0)
		return -1;
	return 31 - __builtin_clz(val);
	cout << "start" << ",";
	cout << val << ",";
	cout <<  31 - __builtin_clz(val) << ",";;
    int ret = -1;
    while (val > 0) {
        val >>= 1;
        ret++;
    }
		cout << ret << endl;
		cout << endl;
    return ret;
}

class LCA {
	private:
	Node *tree;
	vector<int> E;		// preorder numbers of euler tour
	vector<int> L;		// levels of euler tour
	vector<int> H;		// first occurence of a preorder number in E
	vector<int> T;    // real preorder to internal preorder mapping
	vector<Node *> N;	// preorder to node mapping
	vector<vector<int> > RMQ;	// precomputed RMQ values

	public:
	LCA(Node *tree) {
		this->tree = tree;
		T = vector<int>();
		if (tree->get_preorder_number() == -1)
			tree->preorder_number();
		euler_tour(tree, 0);
		precompute_rmq();
	}

	LCA() {
		this->tree = NULL;
	}
	void euler_tour(Node *node, int depth) {
		// First visit
		int preorder_number = N.size();
		int euler_number = E.size();
		N.push_back(node);
		if (T.size() <= node->get_preorder_number())
			T.resize(node->get_preorder_number()+1,-1);
		T[node->get_preorder_number()] = preorder_number;

		//cout << preorder_number << "\t";
		//node->print_subtree();
		H.push_back(euler_number);
		L.push_back(depth);
		E.push_back(preorder_number);

		// TODO: check that this is correct for multifurcating trees
		
		list<Node *>::const_iterator c;
		for(c = node->get_children().begin(); c != node->get_children().end();
				c++) {
			euler_tour(*c, depth+1);
			// Middle/Last visit
			L.push_back(depth);
			E.push_back(preorder_number);
		}

	}

	void precompute_rmq() {
		// E gives queries of length 1
		RMQ.push_back(E);
		for(int length = 1; length < E.size(); length *= 2) {
			vector<int> V = vector<int>();
			for(int start = 0; start < E.size() - length; start++) {
				V.push_back(rmq(start, start + length));
			}
			RMQ.push_back(V);
		}
	}
	//vector<vector<int> > RMQ;	// precomputed RMQ values

	// find the index of the rmq between indices i and j of E
	int rmq(int i, int j) {
		if (i == j)
			return E[i];
		int length = j - i - 1;
		//cout << "  " << i << " " << j << endl;
		//cout << "  " << length << endl;
		//cout << "  " << mylog2(length) << endl;
		length = (mylog2(length));
		//cout << "  " << length << endl;
		//cout << "  " << j - (1 << (length)) << endl;
		int rmq1 = RMQ[length+1][i];
		int rmq2;
		if (length >= 0)
			rmq2 = RMQ[length+1][j - (1 << (length))];
		else
			rmq2 = RMQ[length+1][j];
		if (rmq1 < rmq2)
			return rmq1;
		return rmq2;
	}
	Node *get_lca(Node *a, Node *b) {
		int preorder_a = T[a->get_preorder_number()];
		int preorder_b = T[b->get_preorder_number()];
		int lca_index;
		if (preorder_a <= preorder_b)
			lca_index = rmq(H[preorder_a], H[preorder_b]);
		else
			lca_index = rmq(H[preorder_b], H[preorder_a]);
		return N[lca_index];
	}

	Node *get_tree() {
		return tree;
	}

	/* copy constructor
	LCA(const LCA &n) {
	}
	*/
	void debug() {
		for(int i = 0; i < E.size(); i++) {
			cout << " " << E[i];
		}
		cout << endl;
		cout << endl;
		for(int i = 0; i < L.size(); i++) {
			cout << " " << L[i];
		}
		cout << endl;
		cout << endl;
		for(int i = 0; i < H.size(); i++) {
			cout << " " << H[i];
		}
		cout << endl;
		cout << endl;
		cout << endl;
		for(int i = 0; i < RMQ.size(); i++) {
			for(int j = 0; j < RMQ[i].size(); j++) {
			cout << " " << RMQ[i][j];
			}
			cout << endl;
			cout << endl;
		}
	}
};

#endif

// ===== END heuristic/clean_refactor/h36_whidden/LCA.h =====
#include <map>
#include <limits>
//#include "ClusterInstance.h"

using namespace std;

bool MULTI_CLUSTER = false;

class ClusterInstance;
extern bool LEAF_REDUCTION2;

class Forest {
	public:
		vector<Node *> components;
		vector<Node *> deleted_nodes;
		bool rho;
		Forest *twin;
		ClusterInstance *cluster;
                int max_preorder;
	public:
	Forest() {
		init(vector<Node *>());
	}
	Forest(vector<Node *> components) {
		init(components);
	}
	Forest(Node *head) {
		components = vector<Node *>();
		components.push_back(new Node(*head));
		deleted_nodes = vector<Node *>();
		rho = false;
		twin = NULL;
		cluster = NULL;
	}
	Forest(const Forest &f) {
		components = vector<Node *>(f.components.size());
		for(int i = 0; i < f.components.size(); i++) {
			//if (f.components[i] != NULL)
			components[i] = new Node(*f.components[i]);
		}
		deleted_nodes = vector<Node *>();
		rho = f.rho;
		twin = NULL;
		cluster = f.cluster;
		//label_nodes_with_forest();
	}

	Forest(Forest *f) {
		components = vector<Node *>(f->components.size());
		for(int i = 0; i < f->components.size(); i++) {
			//if (f->components[i] != NULL)
			components[i] = new Node(*f->components[i]);
		}
		deleted_nodes = vector<Node *>();
		rho = f->rho;
		twin = NULL;
		cluster = f->cluster;
		//label_nodes_with_forest();
	}
        /* Temporary until rSPR_branch_and_bound_mult_hlpr uses the UndoMachine */
        Forest(Forest *f, map<Node*, Node*> *node_map) {
		components = vector<Node *>(f->components.size());
		for(int i = 0; i < f->components.size(); i++) {
			if (f->components[i] == NULL) {
				components[i] = NULL;
				continue;
			}
		        Node* new_node = new Node(*f->components[i], node_map);
		        components[i] = new_node;
			node_map->emplace(pair<Node*,Node*>(f->components[i], new_node));
		}
		deleted_nodes = vector<Node *>();
		rho = f->rho;
		twin = NULL;
		cluster = f->cluster;
		max_preorder = f->max_preorder;
		//label_nodes_with_forest();
	}

  
	Forest(Forest *f, bool b) {
		components = vector<Node *>(f->components.size());
		for(int i = 0; i < f->components.size(); i++) {
			//if (f->components[i] != NULL)
			components[i] = new Node(*f->components[i]);
		}
		deleted_nodes = vector<Node *>();
		rho = f->rho;
		twin = NULL;
		cluster = f->cluster;
		//label_nodes_with_forest();
	}


	void init(vector<Node *> components) {
		this->components = vector<Node *>(components);
		deleted_nodes = vector<Node *>();
		rho = false;
		for(int i = 0; i < components.size(); i++) {
				if (components[i]->str() == "p")
					rho = true;
		}
		twin = NULL;
		cluster = NULL;
	}
	~Forest() {
		for(int i = 0; i < components.size(); i++) {
			if (components[i] != NULL) {
				components[i]->delete_tree();
				components[i] = NULL;
			}
		}
		for(int i = 0; i < deleted_nodes.size(); i++) {
			//if (deleted_nodes[i] != NULL) {
				deleted_nodes[i]->delete_tree();
				deleted_nodes[i] = NULL;
			//}
		}
	} 

	// swap the contents of two forests
	void swap(Forest *f) {
		vector<Node *> components_temp = this->components;
		this->components = f->components;
		f->components = components_temp;
		
//		/*
		vector<Node *> deleted_nodes_temp = this->deleted_nodes;
		this->deleted_nodes = f->deleted_nodes;
		f->deleted_nodes = deleted_nodes_temp;
//		*/

		bool rho_temp = this->rho;
		this->rho = f->rho;
		f->rho = rho_temp;

		ClusterInstance *c_temp = this->cluster;
		this->cluster = f->cluster;
		f->cluster = c_temp;
	}

	// print the forest
	void print_components() {
		vector<Node *>::iterator it = components.begin();
		for(it = components.begin(); it != components.end(); it++) {
			Node *root = *it;
			if (root == NULL)
				cout << "!";
			else if (root->is_leaf() && root->str() == "")
				cout << "*";
			else
				root->print_subtree_hlpr();
			cout << " ";
		}
		cout << endl;
	}

	// print the components seperated by s
	void print_components(string s) {
		vector<Node *>::iterator it = components.begin();
		for(it = components.begin(); it != components.end(); it++) {
			Node *root = *it;
			if (root == NULL)
				cout << "!";
			else
				root->print_subtree_hlpr();
			cout << s;
		}
		cout << endl;
	}

	// print the forest showing twins
	void print_components_with_twins() {
		vector<Node *>::iterator it = components.begin();
		for(it = components.begin(); it != components.end(); it++) {
			Node *root = *it;
			if (root == NULL)
				cout << "!";
			else
				root->print_subtree_twin_hlpr();
			cout << " ";
		}
		cout << endl;
	}

	// print the forest showing edge intervals
	void print_components_with_edge_pre_interval() {
		vector<Node *>::iterator it = components.begin();
		for(it = components.begin(); it != components.end(); it++) {
			Node *root = *it;
			if (root == NULL)
				cout << "!";
			else
				cout << root->str_edge_pre_interval_subtree();
			cout << " ";
		}
		cout << endl;
	}

	// return the string for this forest
	string str() {
		string s = "";
		vector<Node *>::iterator it;
		for(it = components.begin(); it != components.end(); it++) {
			Node *root = *it;
			if (root == NULL)
				s += "!";
			else
				s += root->str_subtree();
			s += " ";
		}
		return s;
	}

	inline void add_component(Node *head) {
		components.push_back(head);
	}
	inline void add_component(int pos, Node *head) {
		components.insert(components.begin() + pos, head);
	}
	inline void add_deleted_node(Node *n) {
		deleted_nodes.push_back(n);
	}

	inline Node *get_component(int i) {
		return components[i];
	}
	void set_component(int i, Node *head) {
		components[i] = head;
	}
	inline int num_components() {
		return components.size();
	}
	void set_twin(Forest *f) {
		twin = f;
	}
	Forest *get_twin() {
		return twin;
	}
	ClusterInstance *get_cluster() {
		return cluster;
	}
	void set_cluster(ClusterInstance *c) {
		cluster = c;
	}
	void set_cluster(ClusterInstance &c) {
		cluster = &c;
	}
	void update_component(Node *old_c, Node *new_c) {
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			Node *component = *i;
			if (&(*component) == &(*old_c))
				*i = new_c;
		}
	}

	// return a list of the sibling pairs
	list<Node *> *find_sibling_pairs() {
		list<Node *> *sibling_pairs = new list<Node *>();
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			Node *component = *i;
			component->append_sibling_pairs(sibling_pairs);
		}
		return sibling_pairs;
	}

  	// return a list of the sibling groups
	list<Node *> *find_sibling_groups() {
		list<Node *> *sibling_groups = new list<Node *>();
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			Node *component = *i;
			component->append_sibling_groups(sibling_groups);
		}
		return sibling_groups;
	}

  
	// return a deque of the singleton leaves
	list<Node *> find_singletons() {
		list<Node *> singletons = list<Node *>();
		vector<Node *>::iterator i = components.begin();
		// TODO: is this a hack or correct? We don't want the first
		// component to be a singleton because of rho!
		i++;
		for(; i != components.end(); i++) {
			Node *component = *i;
			if (component->is_leaf()) {
				singletons.push_back(component);
			}
		}
		return singletons;
	}

	// make nodes pointed to in the forest point back
	void resync() {
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			(*i)->resync();
		}
	}
	// clear twin pointers
	void unsync() {
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			(*i)->unsync();
		}
	}
	// clear interior twin pointers
	void unsync_interior() {
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			(*i)->unsync_interior();
		}
	}
	// clear twin pointers
	Node *find_by_prenum(int prenum) {
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			Node *ans = (*i)->find_by_prenum(prenum);
			if (ans != NULL)
				return ans;
		}
		return NULL;
	}

void labels_to_numbers(map<string, int> *label_map, map<int, string> *reverse_label_map) {
	vector<Node *>::iterator i;
	for(i = components.begin(); i != components.end(); i++) {
		(*i)->labels_to_numbers(label_map, reverse_label_map);
	}
}

void numbers_to_labels(map<int, string> *reverse_label_map) {
	vector<Node *>::iterator i;
	for(i = components.begin(); i != components.end(); i++) {
		(*i)->numbers_to_labels(reverse_label_map);
	}
}

int size() {
	return components.size();
}

bool add_rho() {
	if (rho)
		return false;
	Node *T_p = new Node("p");
	T_p->set_preorder_number(-1);
	add_component(T_p);
	rho = true;
	return true;
}

void set_rho(bool b) {
	rho = b;
}

bool contains_rho() {
	return rho;
}

void erase_components(int start, int end) {
	components.erase(components.begin()+start, components.begin()+end);
}
void erase_components() {
	components.clear();
}

// tell the nodes what forest they are in
void label_nodes_with_forest() {
	for(int i = 0; i < size(); i++) {
		get_component(i)->set_forest_rec(this);
	}
}

void move_first_component_to_end() {
	Node *temp = components[0];
	components[0] = components[components.size()-1];
	components[components.size()-1] = temp;
}
void unprotect_edges() {
	for(int i = 0; i < size(); i++) {
		get_component(i)->unprotect_subtree();
	}
}

};

// Functions

vector<Node *> find_labels(vector<Node *> components);
bool sync_twins(Forest *T1, Forest *T2);
void sync_interior_twins(Forest *T1, Forest *T2);
void sync_interior_twins(Node *n, LCA *twin_LCA);
void sync_interior_twins(Node *n, vector<LCA> *F2_LCAs);
list<Node *> *find_cluster_points(Forest *F1, Forest *F2);
void find_cluster_points(Node *n, list<Node *> *cluster_points,
		vector<int> *leaf_counts_F1, vector<int> *leaf_counts_F2);
void delete_and_merge_LCAs(list<Node *> *active_descendants,
		vector<LCA> *F2_LCAs, list<Node *>:: iterator node1_location,
		list<Node *>:: iterator node2_location);
void delete_and_merge_LCAs(Node *n, list<Node *> *active_descendants,
		vector<LCA> *F2_LCAs);



// Make the leaves of two forests point to their twin in the other tree
// Note: removes unique leaves
bool sync_twins(Forest *T1, Forest *T2) {
	vector<Node *> T1_labels = vector<Node *>();
	vector<Node *> T2_labels = vector<Node *>();
	vector<Node *> T1_components = T1->components;
	vector<Node *> T2_components = T2->components;
	vector<Node *>::iterator i;
	Node *T1_rho = NULL;
	Node *T2_rho = NULL;
	for(i = T1_components.begin(); i != T1_components.end(); i++) {
		Node *component = *i;
		vector<Node *> unsorted_labels = component->find_leaves();
		vector<Node *>::iterator j;
		for(j = unsorted_labels.begin(); j != unsorted_labels.end(); j++) {
			Node *leaf = *j;
//			cout << "T1: " << leaf->str() << endl;
			if (leaf->str() == "p") {
				T1_rho = leaf;
			}
			else {
				// find smallest number contained in the label
				int number = stomini(leaf->str());
//				cout << "\t" << number << endl;
				if (number < INT_MAX) {
					if (number >= T1_labels.size())
						T1_labels.resize(number+1, 0);
					T1_labels[number] = leaf;
				}
			}
		}
	}
	for(i = T2_components.begin(); i != T2_components.end(); i++) {
		Node *component = *i;
		vector<Node *> unsorted_labels = component->find_leaves();
		vector<Node *>::iterator j;
		for(j = unsorted_labels.begin(); j != unsorted_labels.end(); j++) {
			Node *leaf = *j;
//			cout << "T2: " << leaf->str() << endl;
			if (leaf->str() == "p") {
				T2_rho = leaf;
			}
			else {
				// find smallest number contained in the label
				int number = stomini(leaf->str());
//				cout << "\t" << number << endl;
				if (number < INT_MAX) {
					if (number >= T2_labels.size())
						T2_labels.resize(number+1, 0);
					T2_labels[number] = leaf;
				}
			}
		}
	}
	T1_labels.resize(T1_labels.size()+1);
	T1_labels[T1_labels.size()-1]=T1_rho;
	T2_labels.resize(T2_labels.size()+1);
	T2_labels[T2_labels.size()-1]=T2_rho;

	int size = T1_labels.size();
	if (size > T2_labels.size())
		size = T2_labels.size();
//	cout << "Syncing Twins" << endl;
	for(int i = 0; i < size; i++) {
		Node *T1_a = T1_labels[i];
		Node *T2_a = T2_labels[i];
		if (T1_a == NULL && T2_a != NULL) {
			Node *node = T2_a->parent();
			if (node == NULL)
				return false;
			int numc = node->get_children().size();
			if (node->parent() == NULL && node->lchild()->is_leaf() &&
					(numc == 1 || (numc == 2 && node->rchild()->is_leaf()))) {
				return false;
				Node *sibling = node->lchild();
				if (sibling == T2_a)
						sibling = node->rchild();
				T2_labels[stomini(sibling->str())] = sibling;
			}
			delete T2_a;
			if (node->get_children().size() < 2) {
				if (node->get_children().size() == 1)
					node->lchild()->lost_child();
				node = node->contract(true);
			}
		}
		else if (T2_a == NULL && T1_a != NULL) {
			Node *node = T1_a->parent();
			if (node == NULL)
				return false;
			int numc = node->get_children().size();
			if (node->parent() == NULL && node->lchild()->is_leaf() &&
					(numc == 1 || (numc == 2 && node->rchild()->is_leaf()))) {
				return false;
				Node *sibling = node->lchild();
				if (sibling == T1_a)
						sibling = node->rchild();
				T1_labels[stomini(sibling->str())] = sibling;
			}
			delete T1_a;
			if (node->get_children().size() < 2) {
				if (node->get_children().size() == 1)
					node->lchild()->lost_child();
				node = node->contract(true);
			}
			
		}
		if (T1_a != NULL && T2_a != NULL) {
			T1_a->set_twin(T2_a);
			T2_a->set_twin(T1_a);
//			cout << T1_a->str() << endl;
		}
	}
	for(int i = size; i < T1_labels.size(); i++) {
		Node *T1_a = T1_labels[i];
		if (T1_a != NULL) {
			Node *node = T1_a->parent();
			if (node == NULL)
				return false;
			int numc = node->get_children().size();
			if (node->parent() == NULL && node->lchild()->is_leaf() &&
					(numc == 1 || (numc == 2 && node->rchild()->is_leaf()))) {
				return false;
				Node *sibling = node->lchild();
				if (sibling == T1_a)
						sibling = node->rchild();
				T1_labels[stomini(sibling->str())] = sibling;
			}
			delete T1_a;
			if (node->get_children().size() < 2) {
				if (node->get_children().size() == 1)
					node->lchild()->lost_child();
				node = node->contract(true);
			}
			
		}
	}
	for(int i = size; i < T2_labels.size(); i++) {
		Node *T2_a = T2_labels[i];
		if (T2_a != NULL) {
			Node *node = T2_a->parent();
			if (node == NULL)
				return false;
			int numc = node->get_children().size();
			if (node->parent() == NULL && node->lchild()->is_leaf() &&
					(numc == 1 || (numc == 2 && node->rchild()->is_leaf()))) {
				return false;
				Node *sibling = node->lchild();
				if (sibling == T2_a)
						sibling = node->rchild();
				T2_labels[stomini(sibling->str())] = sibling;
			}
			delete T2_a;
			if (node->get_children().size() < 2) {
				if (node->get_children().size() == 1)
					node->lchild()->lost_child();
				node = node->contract(true);
			}
		}
	}
//	if (T1_loss != NULL)
//		*T1_loss = T1->get_component(0)->count_lost_children_subtree();
//		- T1->get_component(0)->num_lost_children();;
//	cout << "T1_lost = " << T1_lost << endl;
//	if (T2_loss != NULL)
//		*T2_loss = T2->get_component(0)->count_lost_children_subtree();
//		- T2->get_component(0)->num_lost_children();;
//	cout << "T2_lost = " << T2_lost << endl;

//	cout << "Finished Syncing Twins" << endl;
	return true;
}

/* make interior nodes point to the lca of their descendants in the other tree
   assumes that sync_twins has already been called
   assumes that component 1 of T1 matches with 1 of T2
      NOTE: this isn't true during the algorithm so this will need to be changed
      if we want to interleave clustering. It should be just component 1 of T1
	  matching multiple components of T2 (The first several components?)
   */
void sync_interior_twins(Forest *T1, Forest *T2) {
	Node  *root1 = T1->get_component(0);
	Node  *root2 = T2->get_component(0);
	LCA T1_LCA = LCA(root1);
	LCA T2_LCA = LCA(root2);
	sync_interior_twins(root1, &T2_LCA);
	sync_interior_twins(root2, &T1_LCA);
}

/* make interior nodes point to the LCA of their descendants in the other
	 forest if there is one unambiguous LCA
	 * This is true for a node n of T1 if all leaves that are a descendant
	 	of T1 either map to a single component of F2 or are from another
	 	component of F2 such that the root of that component maps to
	 	a descendant of n (i.e. a finished component)
   * assumes that sync_twins has already been called
	 */
// TODO: initializing parameters seems to be slow
void sync_interior_twins_real(Forest *T1, Forest *F2) {
	Node  *T1_root = T1->get_component(0);
	LCA T1_LCA = LCA(T1_root);
	int T1_size = T1_root->size_using_prenum();
	// roots of F2
	vector<Node *> F2_roots = vector<Node *>();
	// LCA queries for F2
	vector<LCA> F2_LCAs = vector<LCA>();
	// lists of root nodes that map to a given T1 node
	T1_root->initialize_root_lcas(list<Node *>());
	// list of active descendants
	T1_root->initialize_active_descendants(list<Node *>());

	// should be fine.
	for(int i = 0; i < F2->num_components(); i++) {
//		cout << "starting i" << endl;
		F2_roots.push_back(F2->get_component(i));
#ifdef DEBUG_SYNC
		cout << "COMPONENT " << i << ": " << F2_roots[i]->str_subtree() << endl;
		cout << "foo" << endl;
		cout << F2_roots[i]->get_twin() << endl;
		if (F2_roots[i]->get_twin() != NULL) {
		cout << F2_roots[i]->get_twin()->str_subtree() << endl;
		cout << F2_roots[i]->get_twin()->parent() << endl;
		cout << "fooa" << endl;
		if (F2_roots[i]->get_twin()->parent() != NULL) {
		cout << F2_roots[i]->get_twin()->parent()->str_subtree() << endl;
		}
		cout << "foob" << endl;
		}
#endif
		// ignore finished components
//		if (F2_roots[i]->get_twin() != NULL && F2_roots[i]->get_twin()->parent() == NULL) {
//			cout << "fooc" << endl;
		//F2_LCAs.push_back(F2_roots[i]);
//			F2_LCAs.push_back(LCA(F2_roots[i]));
//			cout << "fooc" << endl;
//			continue;
//		}
//		cout << "foo" << endl;
		// ignore rho components
		if (F2_roots[i]->str() == "p") {
		//	F2_LCAs.push_back(LCA());
		F2_LCAs.push_back(F2_roots[i]);
		//F2_LCAs.push_back(NULL);//LCA(F2_roots[i]));
			continue;
		}
//		cout << "aa" << endl;
		F2_LCAs.push_back(LCA(F2_roots[i]));
		// number the component
		F2_roots[i]->initialize_component_number(i);
		// list of nodes that get deleted when a component is finished
		F2_roots[i]->initialize_removable_descendants(list<list<Node *>::iterator>());
		// sync the component with T1
		if (F2_roots[i]->str() != "p" &&
				!(F2_roots[i]->get_twin() != NULL && F2_roots[i]->get_twin()->parent() == NULL)) {
			sync_interior_twins(F2_roots[i], &T1_LCA);
		}
		// keep reverse pointer for the root's twin
//		cout << "a" << endl;
		/*
		cout << F2_roots[i] << endl;
		cout << F2_roots[i]->str_subtree() << endl;
		cout << F2_roots[i]->get_twin() << endl;
		cout << F2_roots[i]->get_twin()->str_subtree() << endl;
		cout << F2_roots[i]->get_twin()->get_parameter_ref(ACTIVE_DESCENDANTS) << endl;
		cout << F2_roots[i]->get_twin()->get_parameter_ref(ROOT_LCAS) << endl;
		cout << boost::any_cast<list<Node *> >(F2_roots[i]->get_twin()->get_parameter_ref(ROOT_LCAS))->size() << endl;
		*/
		if (i > 0 || T1->contains_rho())
			F2_roots[i]->get_twin()->get_root_lcas()->push_back(F2_roots[i]);
		else
			T1_root->get_root_lcas()->push_back(F2_roots[i]);
//		cout << "b" << endl;
	}
//	cout << "syncing" << endl;
	sync_interior_twins(T1_root, &F2_LCAs); 
}

/* make interior nodes point to the lca of their descendants in the other
 * tree
 * assumes that sync_twins has already been called
 */
void sync_interior_twins(Node *n, LCA *twin_LCA) {
	list<Node *>::iterator c = n->get_children().begin();
	if (c == n->get_children().end())
		return;
	sync_interior_twins(*c, twin_LCA);
	n->set_twin((*c)->get_twin());
	c++;
	while(c != n->get_children().end()) {
		sync_interior_twins(*c, twin_LCA);
		Node *twin = twin_LCA->get_lca(n->get_twin(), (*c)->get_twin());
		n->set_twin(twin);
		c++;
	}
}

void sync_interior_twins(Node *n, vector<LCA> *F2_LCAs) {
	Node *lc = n->lchild();
	Node *rc = n->rchild();
	list<Node *> *active_descendants = n->get_active_descendants();
	// visit children first
	list<Node *>::iterator c;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		sync_interior_twins(*c, F2_LCAs);
	}
	#ifdef DEBUG_SYNC
	cout << "SYNC_INTERIOR_TWINS()" << endl;
	cout << n->str_subtree() << endl;
	#endif
	if (n->get_children().size() == 0) {
//		cout << "leaf" << endl;
		active_descendants->push_back(n->get_twin());
		list<Node *>::iterator node_location = active_descendants->end();
		node_location--;
			n->get_twin()->get_removable_descendants()->push_back(node_location);
	}
	// no rc so propogate up
	if (n->get_children().size() == 1) {
		Node *lc = n->get_children().front();
//		cout << "no rc" << endl;
		n->set_twin(lc->get_twin());
		list<Node *> *lc_active_descendants = lc->get_active_descendants();
		active_descendants->splice(active_descendants->end(),*lc_active_descendants);
	}
	// TODO: generalize from here for 2 or more children
	// two children so put their info together
	else if (lc != NULL && rc != NULL) {
//		cout << "two children" << endl;
		list<Node *> *lc_active_descendants = lc->get_active_descendants();
		list<Node *> *rc_active_descendants = rc->get_active_descendants();

/*	#ifdef DEBUG_SYNC
	cout << "active_descendants lc" << endl;
	for(list<Node *>::iterator i =  lc_active_descendants-> begin(); i != lc_active_descendants->end(); i++) {
		cout << "\t" << (*i)->str_subtree() << endl;
	}
		cout << endl;
	cout << "active_descendants rc" << endl;
	for(list<Node *>::iterator i =  rc_active_descendants-> begin(); i != rc_active_descendants->end(); i++) {
		cout << "\t" << (*i)->str_subtree() << endl;
	}
		cout << endl;
	#endif
*/

		vector<list<Node *>::iterator> node_location =
				vector<list<Node *>::iterator>();
		list<Node *>::iterator node1_location;
		int nonempty_active_descendants_count = 0;
		for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
			if (!(*c)->get_active_descendants()->empty()) {
				nonempty_active_descendants_count++;
				if (nonempty_active_descendants_count > 1) {
					node1_location = active_descendants->end();
					node1_location--;
					node_location.push_back(node1_location);
				}
//		cout << active_descendants->size() << endl;
				active_descendants->splice(active_descendants->end(),
						*((*c)->get_active_descendants()));
//		cout << active_descendants->size() << endl;
			}
		}

		/* check the intersection points to see if we have two
			leaves from the same component
		*/
//		cout << "foo" << endl;
		#ifdef DEBUG_SYNC
		cout << active_descendants->size() << endl;
		#endif
		for(int i = 0; i < node_location.size(); i++) {
			list<Node *>::iterator node1_location = node_location[i];
			list<Node *>::iterator node2_location = node1_location;
			node2_location++;
			delete_and_merge_LCAs(active_descendants, F2_LCAs, node1_location,
					node2_location);
		}
		#ifdef DEBUG_SYNC
		cout << "done first merge" << endl;
		#endif

		/* check to see if n is twinned by a root of F2
			 if so, then remove each leaf twinned by that component
			 and check each of the new intersection points
		*/
		list<Node *> *root_lcas = n->get_root_lcas();
		while(!root_lcas->empty()) {

			Node *root_lca = root_lcas->front();
			root_lcas->pop_front();
			/* TODO: problem when n is a root
				 We don't care about this but it might mean there is a different
				 problem
				 */
			if (n->parent() != NULL) {
			#ifdef DEBUG_SYNC
			cout << "deleting from component " << endl;
			#endif
				delete_and_merge_LCAs(root_lca, active_descendants, F2_LCAs);
			}
		}
		#ifdef DEBUG_SYNC
		cout << "done checking component" << endl;
		#endif

		/* If we have a single element in n's active descendants
			 list then set twin pointers appropriately
		*/

	#ifdef DEBUG_SYNC
	cout << "active_descendants done" << endl;
	for(list<Node *>::iterator i =  active_descendants->begin(); i != active_descendants->end(); i++) {
		cout << "\t" << (*i)->str_subtree() << endl;
	}
	#endif
		if (active_descendants->size() == 1) {
			#ifdef DEBUG_SYNC
			cout << "found twin" << endl;
			#endif
			Node *twin = active_descendants->front();
			n->set_twin(twin);
		}
	}
	if (n->parent() == NULL)
		active_descendants->clear();
}

void sync_af_twins(Forest *F1, Forest *F2) {
	F1->unsync();
	F2->unsync();
	sync_twins(F1, F2);
	for(int i = 0; i < F1->num_components(); i++) {
		F1->get_component(i)->sync_af_twins();
	}
}

/* merge two nodes from a list into their LCA if they are from
	 the same component
	 */
void delete_and_merge_LCAs(list<Node *> *active_descendants,
		vector<LCA> *F2_LCAs, list<Node *>:: iterator node1_location,
		list<Node *>:: iterator node2_location) {
#ifdef DEBUG_SYNC
	cout << "DELETE_AND_MERGE FIRST" << endl;
	cout << "active_descendants before" << endl;
	for(list<Node *>::iterator i =  active_descendants-> begin(); i != active_descendants->end(); i++) {
		cout << "\t" << (*i)->str_subtree() << endl;
	}
	cout << endl;
#endif

//	while(active_descendants->size() > 1) {
	Node *n1 = *node1_location;
	Node *n2 = *node2_location;
#ifdef DEBUG_SYNC
	cout << "n1=" << n1->str_subtree() << endl;
	cout << "n2=" << n2->str_subtree() << endl;
#endif
	int component1 = n1->get_component_number();
	int component2 = n2->get_component_number();
#ifdef DEBUG_SYNC
	cout << "c1=" << component1 << endl;
	cout << "c2=" << component2 << endl;
#endif
	if (component1 == component2) {
#ifdef DEBUG_SYNC
		cout << "size=" << (*F2_LCAs).size() << endl;
		for(int i = 0; i < (*F2_LCAs).size(); i++) {
			if ((*F2_LCAs)[i].get_tree() == NULL)
				cout << "\t" << "NULL" << endl;
			else
				cout << "\t" << (*F2_LCAs)[i].get_tree()->str_subtree() << endl;
		}
		cout << (*F2_LCAs)[component1].get_tree()->str_subtree() << endl;
#endif
		Node *lca = (*F2_LCAs)[component1].get_lca(n1,n2);
#ifdef DEBUG_SYNC
		cout << lca->str_subtree() << endl;
#endif
//		cout << "xa" << endl;
		list<Node *>::iterator lca_location =
			active_descendants->insert(node1_location,lca);
//		cout << "xb" << endl;
//		active_descendants->erase(node1_location);
//		cout << "xc" << endl;

// TODO: could this be faster?
		bool remove = false;
		list<list<Node *>::iterator>::iterator i;
		for(i = n1->get_removable_descendants()->begin(); i != n1->get_removable_descendants()->end(); i++) {
			if (*i == node1_location) {
				//active_descendants->erase(*i);
				remove = true;
				break;
			}
		}
		if (remove) {
			active_descendants->erase(*i);
			n1->get_removable_descendants()->erase(i);
		}
//		n1->get_removable_descendants()->clear();
		// TODO: delete each when clearing?
//		cout << "xd" << endl;
//		active_descendants->erase(node2_location);
//		cout << "xe" << endl;
		remove = false;
		for(i = n2->get_removable_descendants()->begin(); i != n2->get_removable_descendants()->end(); i++) {
			if (*i == node2_location) {
				//active_descendants->erase(*i);
				remove = true;
				break;
			}
		}
		if (remove) {
			active_descendants->erase(*i);
			n2->get_removable_descendants()->erase(i);
		}
//		n2->get_removable_descendants()->clear();
//		cout << "xf" << endl;
		lca->get_removable_descendants()->push_back(lca_location);
//		cout << "xg" << endl;
	}
//		else {
//			break;
//		}
//	}
#ifdef DEBUG_SYNC
	cout << "active_descendants after" << endl;
	for(list<Node *>::iterator i =  active_descendants-> begin(); i != active_descendants->end(); i++) {
		cout << "\t" << (*i)->str_subtree() << endl;
	}
#endif
}

/* delete each leaf from the list that is twinned with the component
	 of n. For each such deleted node, merge its predecessor
	 and successor in the list into their LCA if they are from
	 the same component (other than n's component)
	 */
void delete_and_merge_LCAs(Node *n, list<Node *>
		*active_descendants, vector<LCA> *F2_LCAs) {
	int component = n->get_component_number();
	list<list<Node *>::iterator> *removable_descendants	=
			n->get_removable_descendants();

	if (n->lchild() != NULL)
		delete_and_merge_LCAs(n->lchild(), active_descendants, F2_LCAs);
	if (n->rchild() != NULL)
		delete_and_merge_LCAs(n->rchild(), active_descendants, F2_LCAs);
	#ifdef DEBUG_SYNC

	cout << n->str_subtree() << endl;
	cout << "removable_descendants" << endl;
	for(list<list<Node *>::iterator>::iterator i =  removable_descendants-> begin(); i != removable_descendants->end(); i++) {
		cout << "\t" << (**i)->str_subtree() << endl;
	}	
	#endif
//	cout << "foo" << endl;
	while (!removable_descendants->empty()) {
//		cout << "fooa" << endl;
		list<Node *>::iterator leaf_location = removable_descendants->front();
//		cout << "foob" << endl;
		removable_descendants->pop_front();
//		cout << "fooc" << endl;
//		cout << *leaf_location << endl;
//		cout << (*leaf_location)->str_subtree() << endl;
//		cout << "food" << endl;

// TODO: problem here
		if (leaf_location != active_descendants->begin() &&
				leaf_location != active_descendants->end() &&
				leaf_location != -- active_descendants->end()) {
//		if (active_descendants->front() != *leaf_location
//				&& active_descendants->back() != *leaf_location) {
			list<Node *>::iterator node1_location = leaf_location;
//		cout << "fooe" << endl;
			list<Node *>::iterator node2_location = leaf_location;
//		cout << "foof" << endl;
			node1_location--;
//		cout << "foog" << endl;
			node2_location++;
//		cout << "fooh" << endl;
			active_descendants->erase(leaf_location);
//		cout << "fooi" << endl;
			int node1_component = (*node1_location)->get_component_number();
//		cout << "fooj" << endl;
			if (component != node1_component)
				delete_and_merge_LCAs(active_descendants, F2_LCAs, node1_location,
						node2_location);
//		cout << "fook" << endl;
		}
		else {//if (active_descendants->size() > 1){
			active_descendants->erase(leaf_location);
		}

	}
//	cout << "foo end" << endl;
	// TODO: continue to lc and rc?
}

list<Node *> *find_cluster_points(Forest *F1, Forest *F2) {
	list<Node *> *cluster_points = new list<Node *>();
	vector<int> *leaf_counts_F1 = NULL;
	vector<int> *leaf_counts_F2 = NULL;
	if (MULTI_CLUSTER) {
		leaf_counts_F1 = F1->get_component(0)->find_leaf_counts(); 
		leaf_counts_F2 = F2->get_component(0)->find_leaf_counts(); 
	}
	find_cluster_points(F1->get_component(0), cluster_points, leaf_counts_F1,
			leaf_counts_F2);
	if (MULTI_CLUSTER) {
		delete leaf_counts_F1;
//		delete leaf_counts_F2;
	}
	//cout << "foo" << endl;
	return cluster_points;
}

// find the cluster points
void find_cluster_points(Node *n, list<Node *> *cluster_points,
		vector<int> *leaf_counts_F1, vector<int> *leaf_counts_F2) {
//	cout << "Start: " << n->str_subtree() << endl;
	list<Node *>::iterator c;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		find_cluster_points(*c, cluster_points, leaf_counts_F1,
				leaf_counts_F2);
	}
	/*
	cout << "here" << endl;
	cout << "n= " << n->str_subtree() << endl;
	cout << n->get_depth() << endl;
	if (n->get_twin() != NULL) {
		cout << "n_twin= " << n->get_twin()->str_subtree() << endl;
		cout << "n_twin_twin= " << n->get_twin()->get_twin()->str_subtree() << endl;
		cout << n->get_twin()->get_twin()->get_depth() << endl;
	}
	cout << n->parent() << endl;
	*/
	bool is_cluster = true;
	int num_clustered_children = 0;
	if (n->get_twin() == NULL ||
			n->parent() == NULL ||
			n->get_children().size() < 2 ||
			n->get_depth() > n->get_twin()->get_twin()->get_depth())
		is_cluster = false;
	else {
#ifdef RSPR
		if (!LEAF_REDUCTION2){
#endif
			for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
				if ((*c)->get_twin() != NULL &&
						(*c)->get_depth() <= (*c)->get_twin()->get_twin()->get_depth())
				num_clustered_children++;
			}
			if (num_clustered_children == n->get_children().size())
				is_cluster = false;
#ifdef RSPR
		}
#endif
	}
	if (is_cluster) {
//		cout << "added cluster_point" << endl;
		cluster_points->push_back(n);
	}
	// buggy, needs testing, doesn't seem worth it
	else if (MULTI_CLUSTER && n->get_twin() != NULL && n->parent() != NULL &&
			n->get_children().size() >= 2) {
		// TODO: use find_leaf_counts if this works
		Node *n_twin = n->get_twin();
		int num_leaves = (*leaf_counts_F1)[n->get_preorder_number()];
		vector<Node *> chosen = vector<Node *>();
		int chosen_leaves = 0;
		if (n_twin != NULL && n->get_edge_pre_start() > -1 && n->get_edge_pre_end() > -1 && n_twin->get_children().size() > 2) {
//			cout << "foo" << endl;
//			cout << n->str_subtree() << endl;
//			cout << num_leaves << endl;
//			cout << n->get_edge_pre_start() << endl;
//			cout << n->get_edge_pre_end() << endl;
			for(c = n_twin->get_children().begin(); c != n_twin->get_children().end(); c++) {
//				cout << "\t" << (*c)->str_subtree() << endl;
				int c_num_leaves = (*leaf_counts_F2)[(*c)->get_preorder_number()];
//				cout << "\t" << c_num_leaves << endl;
				int c_twin_pre = (*c)->get_twin()->get_preorder_number();
//				cout << "\t" << c_twin_pre << endl;
				if (c_twin_pre >= n->get_edge_pre_start() &&
						c_twin_pre <= n->get_edge_pre_end()) {
//					cout << "yes" << endl;
					chosen.push_back(*c);
					chosen_leaves += c_num_leaves;
				}
			}
			// PROBLEM: the new node should have its own preorder number
			// and its own size
			if (num_leaves == chosen_leaves) {
				Node *new_child = new Node();
				n_twin->add_child(new_child);
				new_child->set_preorder_number(n_twin->get_preorder_number());
				new_child->set_edge_pre_start(n_twin->get_edge_pre_start());
				new_child->set_edge_pre_end(n_twin->get_edge_pre_end());
				n->set_twin(new_child);
				new_child->set_twin(n);
				cluster_points->push_back(n);
				for(int i = 0; i < chosen.size(); i++) {
					new_child->add_child(chosen[i]);
				}
			}

		}
	}
//	cout << "End: " << n->str_subtree() << endl;
}

// swap two forests
void swap(Forest **a, Forest **b) {
	(*a)->swap(*b);
}

// expand all contracted nodes
void expand_contracted_nodes(Forest *F) {
	for(int i = 0; i < F->num_components(); i++) {
		F->get_component(i)->expand_contracted_nodes();
	}
}

Forest *build_finished_forest(string &name) {
	Forest *new_forest = new Forest();
	string::iterator i = name.begin();
		size_t old_loc = 0;
		size_t loc = 0;
		while ((loc = name.find(" ", old_loc)) != string::npos) {
			//cout << "old_loc=" << old_loc << endl;
			//cout << "loc=" << loc << endl;
			//cout << name.substr(old_loc,loc-old_loc) << endl;
			//new_forest->add_component(build_tree(name.substr(old_loc,loc-old_loc)));
			new_forest->add_component(new Node(name.substr(old_loc,loc-old_loc)));
			if (name.substr(old_loc,loc-old_loc) == "p")
				new_forest->set_rho(true);
			//new_forest->print_components();
			old_loc = loc+1;
		}
		new_forest->add_component(new Node(name.substr(old_loc,loc-old_loc)));
		return new_forest;
		//new_forest->add_component(build_tree(name.substr(old_loc,loc-old_loc)));
}
Forest *build_forest(string &name) {
	Forest *new_forest = new Forest();
	string::iterator i = name.begin();
		size_t old_loc = 0;
		size_t loc = 0;
		while ((loc = name.find(" ", old_loc)) != string::npos) {
			//cout << "old_loc=" << old_loc << endl;
			//cout << "loc=" << loc << endl;
			//cout << name.substr(old_loc,loc-old_loc) << endl;
			//new_forest->add_component(build_tree(name.substr(old_loc,loc-old_loc)));
			new_forest->add_component(build_tree(name.substr(old_loc,loc-old_loc)));
			if (name.substr(old_loc,loc-old_loc) == "p")
				new_forest->set_rho(true);
			//new_forest->print_components();
			old_loc = loc+1;
		}
		new_forest->add_component(build_tree(name.substr(old_loc,loc-old_loc)));
		return new_forest;
}
#endif

// ===== END heuristic/clean_refactor/h36_whidden/Forest.h =====

// ===== BEGIN heuristic/clean_refactor/h36_whidden/ClusterForest.h =====
/*******************************************************************************
ClusterForest.h


Copyright 2011-2014 Chris Whidden
cwhidden@dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
March 3, 2014
Version 1.2.1

This file is part of rspr.

rspr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

rspr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rspr.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************/
#ifndef INCLUDE_CLUSTERFOREST

#define INCLUDE_CLUSTERFOREST
// skipped already inlined: heuristic/clean_refactor/h36_whidden/Forest.h

using namespace std;

class ClusterForest: public Forest {
	public:
		vector<Node *> cluster_nodes;

	public:
	ClusterForest() : Forest(){
		init();
	}
	ClusterForest(vector<Node *> components) : Forest(components) {
		init();
	}
	ClusterForest(Node *head) : Forest(head){
		init();
	}
	ClusterForest(const ClusterForest &f) : Forest(f) {
		cluster_nodes = vector<Node *>(f.cluster_nodes);
	}
	ClusterForest(ClusterForest *f, map<Node*, Node*> *node_map)
			: Forest(f, node_map) {
		cluster_nodes = vector<Node *>();
		for (int i = 0; i < (int)f->cluster_nodes.size(); i++) {
			Node *cluster_node = f->cluster_nodes[i];
			if (cluster_node == NULL) {
				cluster_nodes.push_back(NULL);
				continue;
			}
			map<Node*, Node*>::iterator found = node_map->find(cluster_node);
			if (found == node_map->end())
				cluster_nodes.push_back(NULL);
			else
				cluster_nodes.push_back(found->second);
		}
	}

	void init() {
		cluster_nodes = vector<Node *>();
		cluster_nodes.push_back(NULL);
	}
	~ClusterForest() {
		/* don't delete cluster nodes because they are, by
		   definition, in the components
		*/
		cluster_nodes.clear();
	} 

	// add a new cluster at n
	void add_cluster(Node *n, string name) {
		Node *n_parent = n->parent();
		n->cut_parent();
		Node *n_child = new Node(name);
		n_parent->add_child(n_child);
		add_component(n);
		cluster_nodes.push_back(n_child);
	}

	// swap the contents of two forests
	void swap(ClusterForest *f) {
		vector<Node *> components_temp = this->components;
		this->components = f->components;
		f->components = components_temp;
		
		/*
		vector<Node *> deleted_nodes_temp = this->deleted_nodes;
		this->deleted_nodes = f->deleted_nodes;
		f->deleted_nodes = deleted_nodes_temp;
		*/

		vector<Node *> cluster_nodes_temp = this->cluster_nodes;
		this->cluster_nodes = f->cluster_nodes;
		f->cluster_nodes = cluster_nodes_temp;
	}

	void swap(Forest *f) {
		vector<Node *> components_temp = this->components;
		this->components = f->components;
		f->components = components_temp;

		if (contains_rho())
			f->rho = true;
	}

	inline Node *get_cluster_node(int i) {
		return cluster_nodes[i];
	}
	inline int num_clusters() {
		return cluster_nodes.size();
	}
	void join_cluster(int cluster_loc, Forest *solved_cluster) {
		Node *cluster_node = get_cluster_node(cluster_loc);
		Node *cluster_parent = cluster_node->parent();
		cluster_node->cut_parent();
//		delete cluster_node;
		cluster_node->delete_tree();
//		delete components[cluster_loc];
		components[cluster_loc]->delete_tree();
		components[cluster_loc] = NULL;
		int start = 0;
		if (solved_cluster->contains_rho()) {
			if (cluster_parent->parent() != NULL ||
					!((cluster_parent->lchild() && cluster_parent->lchild()->get_name() == "X")
					|| (cluster_parent->rchild() && cluster_parent->rchild()->get_name() == "X")))
				cluster_parent->contract(true);
		}
		else {
			cluster_parent->add_child(solved_cluster->get_component(0));
			//cluster_parent->add_child(new Node(*(solved_cluster->get_component(0))));
			start = 1;
			if (cluster_parent->get_children().size() == 1) {
				cluster_parent->contract(true);
			}
		}
		// should we add these to a finished_components or something?
		for(int i = start; i < solved_cluster->num_components(); i++) {
			if (solved_cluster->get_component(i)->str() != "p")
				add_component(solved_cluster->get_component(i));
				//add_component(new Node(*(solved_cluster->get_component(i))));
			else
				solved_cluster->get_component(i)->delete_tree();
		}
		solved_cluster->erase_components();

	}

	void join_cluster(Forest *solved_cluster) {
		int start = 0;
		if (solved_cluster->contains_rho()) {
			add_rho();
		}
		else {
			add_component(0,solved_cluster->get_component(0));
			//add_component(0,new Node(*(solved_cluster->get_component(0))));
			start = 1;
		}
		// should we add these to a finished_components or something?
		for(int i = start; i < solved_cluster->num_components(); i++) {
			if (solved_cluster->get_component(i)->str() != "p")
				//add_component(new Node(*(solved_cluster->get_component(i))));
				add_component(solved_cluster->get_component(i));
			else
				solved_cluster->get_component(i)->delete_tree();
		}
		solved_cluster->erase_components();
	}

};

// Functions

// swap two cluster forests
void swap(ClusterForest **a, ClusterForest **b) {
	(*a)->swap(*b);
}
#endif

// ===== END heuristic/clean_refactor/h36_whidden/ClusterForest.h =====
// skipped already inlined: heuristic/clean_refactor/h36_whidden/LCA.h

// ===== BEGIN heuristic/clean_refactor/h36_whidden/ClusterInstance.h =====
/*******************************************************************************
ClusterInstance.h


Copyright 2012-2014 Chris Whidden
cwhidden@dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
March 3, 2014
Version 1.2.1

This file is part of rspr.

rspr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

rspr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rspr.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************/
#ifndef INCLUDE_CLUSTERINSTANCE

#define INCLUDE_CLUSTERINSTANCE
// skipped already inlined: heuristic/clean_refactor/h36_whidden/Forest.h
// skipped already inlined: heuristic/clean_refactor/h36_whidden/Node.h

using namespace std;

class ClusterInstance {
	public:
		Forest *F1;
		Forest *F2;
		Node *F1_cluster_node;
		Node *F2_cluster_node;
		bool F2_has_component_zero;

	public:
	ClusterInstance(){
		init(NULL, NULL, NULL, NULL, false);
	}
	ClusterInstance(Forest *f1, Forest *f2, Node *f1_cluster_node,
			Node *f2_cluster_node,bool f2_has_component_zero) {
		init (f1, f2, f1_cluster_node, f2_cluster_node, f2_has_component_zero);
	}
	ClusterInstance(const ClusterInstance &c) {
		init(c.F1,c.F2,c.F1_cluster_node,c.F2_cluster_node,c.F2_has_component_zero);
	}
	~ClusterInstance() {
	}

	void init(Forest *f1, Forest *f2, Node *f1_cluster_node,
			Node *f2_cluster_node, bool f2_has_component_zero) {
		F1 = f1;
		F2 = f2;
		F1_cluster_node = f1_cluster_node;
		F2_cluster_node = f2_cluster_node;
		F2_has_component_zero = f2_has_component_zero;
	}

	bool is_original() {
		if (F1_cluster_node == NULL && F2_cluster_node == NULL)
			return true;
		else
			return false;
	}

	int join_cluster(Forest *original_F1, Forest *original_F2) {
		int adjustment = 0;
		int rho_1 = -1;
		int rho_2 = -1;
		int start = 0;
		Forest *F1_cluster_node_forest = NULL;
		Forest *F2_cluster_node_forest = NULL;
		#ifdef DEBUG_CLUSTERS
			cout << "join_cluster" << endl;
			cout << "original_F1: ";
			original_F1->print_components();
			cout << "original_F2: ";
			original_F2->print_components();
			cout << "cluster_F1: ";
			F1->print_components();
			cout << "cluster_F2: ";
			F2->print_components();
			if (F1_cluster_node != NULL) {
				cout << "F1_cluster_node: ";
				F1_cluster_node->print_subtree();
				cout << "F1_cluster_node_forest: ";
				F1_cluster_node->get_forest()->print_components();
				F1_cluster_node_forest = F1_cluster_node->get_forest();
			}
			if (F2_cluster_node != NULL) {
				cout << "F2_cluster_node: ";
				F2_cluster_node->print_subtree();
				cout << "F2_cluster_node_forest: ";
				F2_cluster_node->get_forest()->print_components();
				F2_cluster_node_forest = F2_cluster_node->get_forest();
			}
		#endif
		if (F1_cluster_node != NULL) {
			F1_cluster_node->decrease_clustered_children();
			F1_cluster_node_forest = F1_cluster_node->get_forest();
		}
		if (F2_cluster_node != NULL)
			F2_cluster_node->decrease_clustered_children();

		if (F1->contains_rho()) {
			bool contract = true;
			if (F1_cluster_node->is_leaf()) {
				if (F1_cluster_node->get_num_clustered_children() >= 1) {
					contract = false;
				}
				else {
					F1_cluster_node->add_child(F1->get_component(0));
					start = 1;
				}
			}
			else {
				F1_cluster_node->get_forest()->add_component(F1->get_component(0));
				start = 1;
			}
			if (F1_cluster_node->parent() == NULL) {
				if (F1_cluster_node->lchild() != NULL &&
					 F1_cluster_node->lchild()->get_num_clustered_children() > 0) {
					F1_cluster_node_forest->update_component(
							F1_cluster_node, F1_cluster_node->lchild());
				}
				if (F1_cluster_node->rchild() != NULL &&
					 F1_cluster_node->rchild()->get_num_clustered_children() > 0) {
					F1_cluster_node_forest->update_component(
							F1_cluster_node, F1_cluster_node->rchild());
				}
			}
			if (contract)
				F1_cluster_node->contract();
		}
		else {
			F1_cluster_node->add_child(F1->get_component(0));
			if (F1_cluster_node->get_num_clustered_children() <= 0)
				F1_cluster_node->contract();
			start = 1;
		}
		// should we add these to a finished_components or something?
		for(int i = start; i < F1->num_components(); i++) {
			if (F1->get_component(i)->str() != "p")
				original_F1->add_component(F1->get_component(i));
			else {
				rho_1 = i;
			}
		}

		bool skip = false;
		Node *F2_cluster = F1->get_component(0)->get_twin();
		// TODO: this is still not quite right. We should only add the cluster
		//to the front if that is where it came from
//		cout << "F1_rho=" << F1->contains_rho() << endl;
//		cout << "F2_rho=" << F2->contains_rho() << endl;
		if ((F1->contains_rho() || F2->contains_rho()) && !F2_has_component_zero && F2_cluster_node != NULL) {
//			if (F2_has_component_zero) {
//				cout << "PROBLEM!!" << endl;
				//if (!original_F1->contains_rho())
				//	original_F1->add_rho();
				//if (!original_F2->contains_rho())
					//original_F2->add_rho();
//			}
			bool contract = true;
			if (F2_cluster_node->is_leaf()) {
				if (F2_cluster_node->get_num_clustered_children() >= 1) {
					contract = false;
				}
				else {
					//F2_cluster_node->add_child(F2->get_component(0));
					F2_cluster_node->add_child(F2_cluster);
					skip = true;
				}
			}
			if (F2_cluster_node->parent() == NULL) {
				if (F2_cluster_node->lchild() != NULL &&
					 F2_cluster_node->lchild()->get_num_clustered_children() > 0) {
					F2_cluster_node->get_forest()->update_component(
							F2_cluster_node, F2_cluster_node->lchild());
				}
				if (F2_cluster_node->rchild() != NULL &&
					 F2_cluster_node->rchild()->get_num_clustered_children() > 0) {
					F2_cluster_node->get_forest()->update_component(
							F2_cluster_node, F2_cluster_node->rchild());
				}
			}
			if (contract) {
				F2_cluster_node->contract();
				if (F1_cluster_node_forest->get_twin()->get_component(0)->str() == F2_cluster->str() && F1_cluster_node_forest->get_component(0)->str() == F2_cluster->str() && !F1_cluster_node_forest->contains_rho()) {
					F1_cluster_node_forest->add_rho();
					F1_cluster_node_forest->get_twin()->add_rho();
					adjustment++;
				}
			}
		}
		else { // if (F2_cluster != NULL) {
			if (F2_cluster == NULL)
				F2_cluster = F2->get_component(0);
			skip = true;
//			cout << "skip=" << skip << endl;
			// TODO: this is not right yet
			if (F2_has_component_zero) {
//				cout << __LINE__ << endl;
				F1_cluster_node_forest->get_twin()->add_component(0, F2_cluster);
				if (F1->contains_rho() || F2->contains_rho()) {
					if (!F1_cluster_node_forest->contains_rho())
						F1_cluster_node_forest->add_rho();
					if (!F1_cluster_node_forest->get_twin()->contains_rho())
						F1_cluster_node_forest->get_twin()->add_rho();
				}
			}
			else if (F2_cluster_node == NULL) {
				skip = false;
			}
			else {
				// TODO: check this
				F2_cluster_node->add_child(F2->get_component(0));
				F2_cluster = F2->get_component(0);
				if ((F2_cluster_node->lchild() == NULL || F2_cluster_node->rchild() == NULL) && F2_cluster_node->get_num_clustered_children() <= 0) {
					F2_cluster_node->contract();
				}
			}
		}
		// should we add these to a finished_components or something?
		for(int i = 0; i < F2->num_components(); i++) {
			if (F2->get_component(i) == F2_cluster) {
				if (!skip) {
					if (F2->get_component(i)->str() != "p")
						F1_cluster_node_forest->get_twin()->add_component(F2->get_component(i));
					else
						rho_2 = i;
				}
			}
			else {
				if (F2->get_component(i)->str() != "p")
					original_F2->add_component(F2->get_component(i));
				else
					rho_2 = i;
			}
		}

		/*
		cout << "Joined" << endl;
		cout << F1_cluster_node << endl;
		cout << F1->get_component(0) << endl;
		if (!F1->contains_rho() && F1_cluster_node != NULL)
			cout << F1_cluster_node->str_subtree() << endl;
		else {
			cout << F1->get_component(0)->str_subtree() << endl;
		}
		if (!F1->contains_rho() && !F2->contains_rho() && F2_cluster_node != NULL)
			cout << F2_cluster_node->str_subtree() << endl;
		else
			cout << F2_cluster->str_subtree() << endl;
			*/

		if (rho_1 >= 0) {
			F1->get_component(rho_1)->delete_tree();
			F1->set_component(rho_1,NULL);
		}
		if (rho_2 >= 0) {
			F2->get_component(rho_2)->delete_tree();
			F2->set_component(rho_2,NULL);
		}
		F1->erase_components();
		F2->erase_components();

		#ifdef DEBUG_CLUSTERS
			cout << "join_cluster finished" << endl;
			cout << "original_F1: ";
			original_F1->print_components();
			cout << "original_F2: ";
			original_F2->print_components();
			if (F1_cluster_node_forest != NULL) {
				cout << "F1_cluster_node_forest: ";
				F1_cluster_node_forest->print_components();
			}
			if (F2_cluster_node_forest != NULL) {
				cout << "F2_cluster_node_forest: ";
				F2_cluster_node_forest->print_components();
			}
		#endif
		return adjustment;
	}

};

void cluster_reduction_find_components(Node *n,
				vector<bool> *F2_cluster_copy_components,
				vector<bool> *old_F2_keep_components,
				int cluster_component_number) {
	Node *lc = n->lchild();
	Node *rc = n->rchild();
	if (lc != NULL)
		cluster_reduction_find_components(lc, F2_cluster_copy_components,
				old_F2_keep_components, cluster_component_number);
	if (rc != NULL)
		cluster_reduction_find_components(rc, F2_cluster_copy_components,
				old_F2_keep_components, cluster_component_number);
	if (lc == NULL && rc == NULL) {
		int cnumber = n->get_twin()->get_component_number();
		if (cnumber != cluster_component_number) {
			(*F2_cluster_copy_components)[cnumber] = true;
			(*old_F2_keep_components)[cnumber] = false;
		}
	}
}

list<ClusterInstance> cluster_reduction(Forest *old_F1, Forest *old_F2,
		list<Node *> *cluster_points) {
	list<ClusterInstance> clusters = list<ClusterInstance>();
	vector<bool> old_F2_keep_components =
		vector<bool>(old_F2->num_components(), true);
	for(list<Node *>::iterator i = cluster_points->begin();
			i != cluster_points->end(); i++) {
		// Cluster F1
		Node *F1_root_node = *i;
		Node *F1_cluster_node = F1_root_node->parent();
//		Node *p = F1_cluster_node;
//		while (p->parent() != NULL)
//			p = p->parent();
//		cout << "root address=" << &(*p) << endl;

		F1_root_node->cut_parent();
		vector<Node *> cluster_F1_components = vector<Node *>();
		cluster_F1_components.push_back(F1_root_node);
		Forest *F1 = new Forest(cluster_F1_components);

		// Cluster F2
		Node *F2_root_node = F1_root_node->get_twin();
		Node *F2_cluster_node = F2_root_node->parent();
		vector<bool> F2_cluster_copy_components =
			vector<bool>(old_F2->num_components(), false);
		int cnumber = F2_root_node->get_component_number();
		bool F2_has_component_zero = false;
		bool skip_F2_cluster = false;
		if (F2_root_node->parent() != NULL) {
			F2_root_node->cut_parent();
		}
		else {
			if (old_F2_keep_components[cnumber] == false) {
				skip_F2_cluster = true;
			}
			else if (old_F2->get_component(cnumber) == F2_root_node){
				old_F2_keep_components[cnumber] = false;		
				if (cnumber == 0) {
					F2_has_component_zero = true;
				}
			}
			else {
//				if (F2_root_node->get_forest()->get_cluster()->F2_has_component_zero == false) {
					F2_root_node->get_forest()->get_cluster()
							->F2_has_component_zero = true;
//					if (F2_root_node->get_forest()->contains_rho() == false)
//						F2_root_node->get_forest()->add_rho();
					F2_cluster_node = F2_root_node->get_forest()->get_cluster()->F2_cluster_node;
//				}
				skip_F2_cluster = true;
			}
		}

		cluster_reduction_find_components(F1_root_node,
				&F2_cluster_copy_components, &old_F2_keep_components, cnumber);
		vector<Node *> cluster_F2_components = vector<Node *>();
		if (!skip_F2_cluster)
			cluster_F2_components.push_back(F2_root_node);
		for(int i = 0; i < F2_cluster_copy_components.size(); i++) {
			if (F2_cluster_copy_components[i] == true)
				cluster_F2_components.push_back(old_F2->get_component(i));
		}
		Forest *F2 = new Forest(cluster_F2_components);
		clusters.push_back(ClusterInstance(F1, F2, F1_cluster_node,
					F2_cluster_node, F2_has_component_zero));
		F1->set_cluster(clusters.back());
		F2->set_cluster(clusters.back());
		if (F1_cluster_node != NULL)
			F1_cluster_node->increase_clustered_children();
		if (F2_cluster_node != NULL)
			F2_cluster_node->increase_clustered_children();
		F1->label_nodes_with_forest();
		F2->label_nodes_with_forest();
		F1->set_twin(F2);
		F2->set_twin(F1);

		for(int i = 0; i < F2_cluster_copy_components.size(); i++) {
			if (F2_cluster_copy_components[i] == true) {
				cluster_F2_components.push_back(
						(old_F2->get_component(i)));
//				cout << "true " ;
			}
//			else
//				cout << "false " ;
		}
//		cout << endl;
	}
	// remove any clustered components from old_F2
	vector<Node *> old_F2_remaining_components = vector<Node *>();
//	cout << "size=" << old_F2->num_components() << endl;
	for(int i = 0; i < old_F2_keep_components.size(); i++) {
		if (old_F2_keep_components[i] == true) {
			old_F2_remaining_components.push_back(
					(old_F2->get_component(i)));
//			cout << "true " ;
		}
//		else
//			cout << "false " ;
	}
//	cout << endl;
//	cout << endl;
	Forest *replace_old_F2 = new Forest(old_F2_remaining_components);
	replace_old_F2->swap(old_F2);
	replace_old_F2->erase_components();
	delete replace_old_F2;


	clusters.push_back(ClusterInstance(old_F1, old_F2, NULL, NULL, true));
	old_F1->set_cluster(clusters.back());
	old_F2->set_cluster(clusters.back());
	old_F1->label_nodes_with_forest();
	old_F2->label_nodes_with_forest();
	old_F1->set_twin(old_F2);
	old_F2->set_twin(old_F1);

	return clusters;

}


#endif

// ===== END heuristic/clean_refactor/h36_whidden/ClusterInstance.h =====

// ===== BEGIN heuristic/clean_refactor/h36_whidden/SiblingPair.h =====
/*******************************************************************************
SiblingPair.h

Data structure for a sibling pair of a binary tree

Copyright 2012-2014 Chris Whidden
cwhidden@dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
March 3, 2014
Version 1.2.1

This file is part of rspr.

rspr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

rspr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rspr.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#ifndef INCLUDE_SIBLINGPAIR

#define INCLUDE_SIBLINGPAIR

#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <set>
// skipped already inlined: heuristic/clean_refactor/h36_whidden/Node.h

using namespace std;

class SiblingPair {
	public:
	Node *a;
	Node *c;
	int key;
	int key2;

	SiblingPair() {
		a = NULL;
		c = NULL;
		key = -1;
		key2 = -1;
	}
	SiblingPair(Node *A, Node *C) {
		a = A;
		c = C;
		key = a->get_preorder_number();
		if (c->get_preorder_number() < key)
			key = c->get_preorder_number();
	//	a->get_preorder_number();
	//	if (c->get_preorder_number() < key)
	//		key = c->get_preorder_number();
	}
	bool operator< (const SiblingPair &sp) const {
		return key < sp.key;
	}
};

	// TODO: binary only
	void find_sibling_pairs_set_hlpr(Node *n,
			set<SiblingPair> *sibling_pairs) {
		Node *lchild = n->lchild();
		Node *rchild = n->rchild();
		bool lchild_leaf = false;
		bool rchild_leaf = false;
		if (lchild != NULL) {
			if (lchild->is_leaf())
				lchild_leaf = true;
			else
				find_sibling_pairs_set_hlpr(lchild,sibling_pairs);
		}
		if (rchild != NULL) {
			if (rchild->is_leaf())
				rchild_leaf = true;
			else
				find_sibling_pairs_set_hlpr(rchild,sibling_pairs);
		}
		if (lchild_leaf && rchild_leaf) {
			sibling_pairs->insert(SiblingPair(lchild,rchild));
			//lchild->add_to_sibling_pairs(sibling_pairs, 1);
			//rchild->add_to_sibling_pairs(sibling_pairs, 2);
		}
	}
	
	// find the sibling pairs in this node's subtree
	void append_sibling_pairs_set(Node *n,set<SiblingPair> *sibling_pairs) {
		find_sibling_pairs_set_hlpr(n,sibling_pairs);
	}

	// find the sibling pairs in this node's subtree
	set<SiblingPair> *find_sibling_pairs_set(Node *n) {
		set<SiblingPair> *sibling_pairs = new set<SiblingPair>();
		find_sibling_pairs_set_hlpr(n,sibling_pairs);
		return sibling_pairs;
	}

	// return a set of the sibling pairs
	set<SiblingPair> *find_sibling_pairs_set(Forest *f) {
		set<SiblingPair> *sibling_pairs = new set<SiblingPair>();
		for(int i = 0; i < f->num_components(); i++) {
			Node *component = f->get_component(i);
			append_sibling_pairs_set(component,sibling_pairs);
		}
		return sibling_pairs;
	}
#endif

// ===== END heuristic/clean_refactor/h36_whidden/SiblingPair.h =====

// ===== BEGIN heuristic/clean_refactor/h36_whidden/UndoMachine.h =====
/*******************************************************************************
UndoMachine.h

Data structure for recording and undoing tree alterations

Copyright 2012-2014 Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
March 3, 2014
Version 1.2.1

This file is part of rspr.

rspr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

rspr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rspr.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/
#ifndef DEFINE_UNDOMACHINE

#define DEFINE_UNDOMACHINE

// skipped already inlined: heuristic/clean_refactor/h36_whidden/Node.h
// skipped already inlined: heuristic/clean_refactor/h36_whidden/Forest.h
#include <list>
#include <typeinfo>


class Node;

class Undoable {
	public:
	virtual ~Undoable() {}
	virtual void undo() = 0;
};

class UndoMachine {
	public:
	list<Undoable *> events;
	int size;

	UndoMachine() {
		events = list<Undoable *>();
		size = 0;
	}

	void add_event(Undoable *event) {
		events.push_back(event);
		size++;
	}

	void insert_event(list<Undoable *>::iterator i, Undoable *event) {
		list<Undoable *>::iterator j = i;
		j++;
		events.insert(j, event);
		size++;
	}

	 list<Undoable *>::iterator get_bookmark() {
		 list<Undoable *>::iterator bookmark = events.end();
		 if (size > 0)
		 	bookmark--;
		 else
			 bookmark = events.begin();
		 return bookmark;
	 }

	void undo() {
		if (!events.empty()) {
			Undoable *event = events.back();
#ifdef DEBUG_UNDO
			cout << typeid(*event).name() << endl;
#endif
			events.pop_back();
			event->undo();
			delete event;
			size--;
		}
	}

	void undo(int num) {
		int temp_num = num;
		while (temp_num > 0) {
			undo();
			temp_num--;
		}
	}

	void undo_all() {
		undo_to(0);
	}

	void undo_to(int to) {
		while (size > to)
			undo();
	}

	void clear_to(int to) {
		while (size > to) {
			Undoable *event = events.back();
			events.pop_back();
			delete event;
			size--;
		}

	}

	int num_events() {
		return size;
	}
};

void ContractEvent(UndoMachine *um, Node *n);

class AddRho : public Undoable {
	public:
	Forest *F;
	AddRho(Forest *f) {
		F = f;
	}
	void undo() {
		F->set_rho(false);
		F->get_component(F->num_components()-1)->delete_tree();
		F->erase_components(F->num_components()-1,F->num_components());
	}
};

class AddComponent : public Undoable {
	public:
	Forest *F;
	AddComponent(Forest *f) {
		F = f;
	}
	void undo() {
		F->erase_components(F->num_components()-1,F->num_components());
	}
};

class AddComponentToFront : public Undoable {
	public:
	Forest *F;
	AddComponentToFront(Forest *f) {
		F = f;
	}
	void undo() {
		F->erase_components(0,1);
	}
};

// TODO: use new insert_child function with a stored successor sibling
// does the end work? maybe a seperate variable for that?
class CutParent : public Undoable {
	public:
	Node *child;
	Node *parent;
	int branch;
	int depth;
	CutParent(Node *c) {
		child = c;
		parent = c->parent();
		branch = 0;
		depth = c->get_depth();
		if (parent != NULL) {
			if (parent->lchild() == child)
				branch = 1;
			else
				branch = 2;
		}
	}

	void undo() {
		if (branch == 1) {
			if (parent->is_leaf())
				parent->add_child(child);
			else
				parent->insert_child(parent->get_children().front(), child);
		}
		else if (branch == 2)
			parent->add_child(child);
		child->set_depth(depth);
	}
};

class ClearSiblingPair : public Undoable {
	public:
		Node *a;
		Node *c;
		int a_status;
		int c_status;
		ClearSiblingPair(Node *x, Node *y) {
			if (x->get_sibling_pair_status() == 1 ||
					y->get_sibling_pair_status() == 2) {
				a = x;
				c = y;
			}
			else {
				a = y;
				c = x;
			}
			a_status = a->get_sibling_pair_status();
			c_status = c->get_sibling_pair_status();
		}

		void undo() {
			a->set_sibling_pair_status(1);
			c->set_sibling_pair_status(2);
			if (a_status == 0) {
				c->set_sibling(a);
			}
			else if (c_status == 0) {
				a->set_sibling(c);
			}

		}
};

class PopClearedSiblingPair : public Undoable {
	public:
		Node *a;
		Node *c;
		int a_status;
		int c_status;

		list<Node *> *sibling_pairs;
		PopClearedSiblingPair(Node *x, Node *y, list<Node *> *s) {
			sibling_pairs = s;
			a = x;
			c = y;
		}

		void undo() {
			sibling_pairs->push_back(c);
			sibling_pairs->push_back(a);
		}
};

class PopSiblingPair : public Undoable {
	public:
		Node *a;
		Node *c;
		list<Node *> *sibling_pairs;
		PopSiblingPair(Node *x, Node *y, list<Node *> *s) {
			sibling_pairs = s;
				a = x;
				c = y;
		}

		void undo() {
			sibling_pairs->push_back(c);
			sibling_pairs->push_back(a);
		}
};

class ContractSiblingPair : public Undoable {
	public:
		Node *node;
		int c1_depth;
		int c2_depth;
		bool binary_node;
		bool node_protected;
		ContractSiblingPair(Node *n) {
			init(n);
			if (n->is_protected())
				node_protected = true;
			else
				node_protected = false;
		}
		ContractSiblingPair(Node *n, Node *child1, Node *child2,
				UndoMachine *um) {
			if (n->get_children().size() == 2)
				init(n);
			else {
				um->add_event(new CutParent(child1));
				um->add_event(new CutParent(child2));
				binary_node = false;
				node = n;
			}
			if (n->is_protected())
				node_protected = true;
			else
				node_protected = false;
		}

		void init(Node *n) {
			node = n;
			if (n->lchild() != NULL)
				c1_depth = n->lchild()->get_depth();
			else
				c1_depth = -1;
			if (n->rchild() != NULL)
				c2_depth = n->rchild()->get_depth();
			else
				c2_depth = -1;
			binary_node = true;
		}

		void undo() {
			if (binary_node) {
				node->undo_contract_sibling_pair();
				if (c1_depth > -1)
					node->lchild()->set_depth(c1_depth);
				if (c2_depth > -1)
					node->rchild()->set_depth(c2_depth);
			}
			if (node_protected)
				node->protect_edge();
		}
};

class AddToFrontSiblingPairs : public Undoable {
	public:
		list<Node *> *sibling_pairs;
		AddToFrontSiblingPairs(list<Node *> *s) {
			sibling_pairs = s;
		}
		void undo() {
			if (!sibling_pairs->empty()) {
				sibling_pairs->pop_front();
				sibling_pairs->pop_front();
			}
		}
};

class AddToSiblingPairs : public Undoable {
	public:
		list<Node *> *sibling_pairs;
		AddToSiblingPairs(list<Node *> *s) {
			sibling_pairs = s;
		}
		void undo() {
			if (!sibling_pairs->empty()) {
				sibling_pairs->pop_back();
				sibling_pairs->pop_back();
			}
		}
};

class AddToSetSiblingPairs : public Undoable {
	public:
		set<SiblingPair> *sibling_pairs;
		SiblingPair pair;
		AddToSetSiblingPairs(set<SiblingPair> *sp, SiblingPair p) {
			sibling_pairs = sp;
			//pair = SiblingPair(p->a,p->c);
			pair = p;
		}
		void undo() {
			if (!sibling_pairs->empty()) {
				sibling_pairs->erase(pair);
			}
		}
};

class RemoveSetSiblingPairs : public Undoable {
	public:
		set<SiblingPair> *sibling_pairs;
		SiblingPair pair;
		RemoveSetSiblingPairs(set<SiblingPair> *sp, SiblingPair p) {
			sibling_pairs = sp;
			pair = p;
		}
		void undo() {
			sibling_pairs->insert(pair);
		}
};

class AddInSiblingPairs : public Undoable {
	public:
		list<Node *> *sibling_pairs;
		int pos;
		AddInSiblingPairs(list<Node *> *s, int p) {
			sibling_pairs = s;
			pos = p;
		}
		void undo() {
			if (!sibling_pairs->empty()) {
				list<Node *>::iterator c = sibling_pairs->begin();
				for(int i = 0; i <= pos && c != sibling_pairs->end(); i++) {
					c++;
				}
				if (c != sibling_pairs->end()) {
					list<Node *>::iterator rem = c;
					c++;
					sibling_pairs->erase(rem);
					rem = c;
					c++;
					sibling_pairs->erase(rem);
				}
			}
		}
};

class SetTwin : public Undoable {
	public:
		Node *node;
		Node *twin;

		SetTwin(Node *n) {
			node = n;
			twin = n->get_twin();
		}

		void undo() {
			node->set_twin(twin);
		}
};

class ChangeName : public Undoable {
	public:
		Node *node;
		string name;

		ChangeName(Node *n) {
			node = n;
			name = string(n->get_name());
		}

		void undo() {
			node->set_name(name);
		}
};

class ChangeEdgePreInterval : public Undoable {
	public:
		Node *node;
		int start;
		int end;

		ChangeEdgePreInterval(Node *n) {
			start = n->get_edge_pre_start();
			end = n->get_edge_pre_end();
			node = n;
		}

		void undo() {
			node->set_edge_pre_start(start);
			node->set_edge_pre_end(end);
		}
};

class ChangePreNum : public Undoable {
	public:
		Node *node;
		int prenum;

		ChangePreNum(Node *n) {
			prenum = n->get_preorder_number();
			node = n;
		}

		void undo() {
			node->set_preorder_number(prenum);
		}
};

class ChangeRightChild : public Undoable {
	public:
		Node *node;
		Node *rchild;
		int rchild_depth;

		ChangeRightChild(Node *n) {
			node = n;
			rchild = n->rchild();
			if (rchild != NULL)
				rchild_depth = rchild->get_depth();
		}

		void undo() {
			if (rchild != NULL) {
				//node->add_child_keep_depth(rchild);
				node->add_child(rchild);
				rchild->set_depth(rchild_depth);
			}
			else
				if (node->rchild() != NULL)
					node->rchild()->cut_parent();
		}
};

class ChangeLeftChild : public Undoable {
	public:
		Node *node;
		Node *lchild;
		int lchild_depth;

		ChangeLeftChild(Node *n) {
			node = n;
			lchild = n->lchild();
			if (lchild != NULL)
				lchild_depth = lchild->get_depth();
		}

		void undo() {
			if (lchild != NULL) {
				//node->add_child_keep_depth(lchild);
				node->add_child(lchild);
				lchild->set_depth(lchild_depth);
			}
			else
				if (node->lchild() != NULL)
					node->lchild()->cut_parent();
		}
};

class AddChild : public Undoable {
	public:
		Node *child;
		int depth;

		AddChild(Node *c) {
			child = c;
			if (c != NULL)
				depth = c->get_depth();
		}

		void undo() {
			if (child != NULL) {
				child->cut_parent();
				child->set_depth(depth);
			}
		}
};

class AddContractedLC : public Undoable {
	public:
		Node *node;

		AddContractedLC(Node *n) {
			node = n;
		}

		void undo() {
			node->set_contracted_lc(NULL);
		}
};

class AddContractedRC : public Undoable {
	public:
		Node *node;

		AddContractedRC(Node *n) {
			node = n;
		}

		void undo() {
			node->set_contracted_rc(NULL);
		}
};

class CreateNode : public Undoable {
	public:
		Node *node;

		CreateNode(Node *n) {
			node = n;
		}

		void undo() {
			if (node != NULL)
				delete node;
		}
};

class ProtectEdge : public Undoable {
	public:
		Node *node;

		ProtectEdge(Node *n) {
			node = n;
		}

		void undo() {
			if (node != NULL)
				node->unprotect_edge();
		}
};

class UnprotectEdge : public Undoable {
	public:
		Node *node;

		UnprotectEdge(Node *n) {
			node = n;
		}

		void undo() {
			if (node != NULL)
				node->protect_edge();
		}
};

class ListPushBack : public Undoable {
	public:
		list<Node *> *l;

		ListPushBack(list<Node *> *x) {
			l = x;
		}

		void undo() {
			l->pop_back();
		}
};

class ListPopBack : public Undoable {
	public:
		list<Node *> *l;
		Node *node;

		ListPopBack(list<Node *> *x) {
			l = x;
			if (!l->empty())
				node = l->back();
			else
				node = NULL;
		}

		void undo() {
			if (node != NULL)
				l->push_back(node);
		}
};


void ContractEvent(UndoMachine *um, Node *n, list<Undoable *>::iterator
		bookmark) {
		Node *parent = n->parent();
		Node *child;
		Node *lc = n->lchild();
		Node *rc = n->rchild();
		Node *ret = NULL;
		// contract out this node and give child to parent
		if (parent != NULL) {
			if (lc && !rc) {
				child = lc;
				um->add_event(new ChangeEdgePreInterval(child));
				um->add_event(new CutParent(child));
				um->add_event(new CutParent(n));
				if (n->is_protected() && !child->is_protected())
					um->add_event(new ProtectEdge(child));
			}
			else if (rc && !lc) {
				child = rc;
				um->add_event(new ChangeEdgePreInterval(child));
				um->add_event(new CutParent(child));
				um->add_event(new CutParent(n));
				if (n->is_protected() && !child->is_protected())
					um->add_event(new ProtectEdge(child));
			}
			else if (lc == NULL && rc == NULL) {
				um->insert_event(bookmark, new CutParent(n));
				parent->delete_child(n);
				ContractEvent(um, parent);
				parent->add_child(n);
			}
		}
		// if no parent then take children of single child and remove it
		else {

			// dead component or singleton, will be cleaned up by the forest
			if (n->get_children().empty()) {
				um->add_event(new ChangeName(n));
			}
			else if (n->get_children().size() == 1) {
				child = n->get_children().front();
//				if (rc == NULL) {
//					um->add_event(new ChangeRightChild(n));
//					child = lc;
//				}
//				else {
//					um->add_event(new ChangeLeftChild(n));
//					child = rc;
//				}
				um->add_event(new CutParent(child));
				/* cluster hack - if we delete a cluster node then
				 * we may try to use it later. This only happens once
				 * per cluster so we can spend linear time to update
				 * the forest
				 */
				if (child->get_num_clustered_children() > 0) {
					//um->add_event(new CutParent(n));
				}
				if (child->get_num_clustered_children() <= 0) {
					// if child is a leaf then get rid of this so we don't lose refs
					// problem: if the child is not c, then we want to copy
					// otherwise we don't
					// copy other parameters and join the twin
					//to this if the child is a label

					Node *new_lc = child->lchild();
					Node *new_rc = child->rchild();
					if (child->is_leaf()) {
						if (child->get_twin() != NULL) {
							um->add_event(new SetTwin(n));
							um->add_event(new SetTwin(child->get_twin()));
						}
						um->add_event(new ChangeName(n));
					}
					um->add_event(new ChangePreNum(n));
					// redundant?
					//um->add_event(new CutParent(child));
					list<Node *>::iterator c;
					for(c = child->get_children().begin();
							c != child->get_children().end();
							c++) {
						um->add_event(new CutParent(*c));
					}
					if (child->get_contracted_lc() != NULL)
						um->add_event(new AddContractedLC(n));
					if (child->get_contracted_rc() != NULL)
						um->add_event(new AddContractedRC(n));
				}
			}
		}

	}

void ContractEvent(UndoMachine *um, Node *n) {
	list<Undoable *>::iterator bookmark = um->get_bookmark();
	ContractEvent(um, n, bookmark);
}

#endif

// ===== END heuristic/clean_refactor/h36_whidden/UndoMachine.h =====

using namespace std;

enum RELAXATION {STRICT, NEGATIVE_RELAXED, ALL_RELAXED};

const string whitespaces = " \t\f\v\n\r";

// note: not using undo
int rSPR_worse_3_mult_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests);
int rSPR_worse_3_mult_approx(Forest *T1, Forest *T2);
int rSPR_worse_3_mult_approx(Forest *T1, Forest *T2, bool sync);

int rSPR_branch_and_bound_mult_range(Forest *T1, Forest *T2, int start_k);
int rSPR_branch_and_bound_mult_range(Forest *T1, Forest *T2, int start_k, int end_k);
int rSPR_branch_and_bound_mult(Forest *T1, Forest *T2, int k);
int rSPR_branch_and_bound_mult_hlpr(Forest *T1, Forest *T2, int k, list<Node*> *sibling_groups, list<Node*> *singletons, Node *protected_stack, list<pair<Forest,Forest>> *AFs, int* num_ties);



int rSPR_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons,
		list<Node *> *sibling_pairs);
int rSPR_3_approx(Forest *T1, Forest *T2);
int rSPR_worse_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests);
int rSPR_worse_3_approx(Forest *T1, Forest *T2);
int rSPR_worse_3_approx(Forest *T1, Forest *T2, bool sync);
int rSPR_worse_3_approx(Node *subtree, Forest *T1, Forest *T2);
int rSPR_worse_3_approx(Node *subtree, Forest *T1, Forest *T2, bool sync);
int rSPR_worse_3_approx_binary_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests);
int rSPR_worse_3_approx_binary(Forest *T1, Forest *T2, bool sync);
int rSPR_worse_3_approx_binary(Forest *T1, Forest *T2);
int rSPR_branch_and_bound(Forest *T1, Forest *T2);
int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k);
int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k,
		map<string, int> *label_map,
		map<int, string> *reverse_label_map);

int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int end_k);
int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int start_k,
		int end_k);
int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
		set<SiblingPair> *sibling_pairs, list<Node *> *singletons, bool cut_b_only,
		list<pair<Forest,Forest> > *AFs, list<Node *> *protected_stack,
		int *num_ties);
int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
		set<SiblingPair> *sibling_pairs, list<Node *> *singletons, bool cut_b_only,
		list<pair<Forest,Forest> > *AFs, list<Node *> *protected_stack,
		int *num_ties, Node *prev_T1_a, Node *prev_T1_c);
int rSPR_total_approx_distance(Node *T1, vector<Node *> &gene_trees);
int rSPR_total_approx_distance(Node *T1, vector<Node *> &gene_trees,
		int threshold);
int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees);
int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees,
		vector<int> *original_scores);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, bool approx);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end, bool approx);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int max_spr);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int max_spr, int start, int end);
void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees);
void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end);
void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, bool approx);
void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end, bool approx);
int rf_total_distance(Node *T1, vector<Node *> &gene_trees);
int rf_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees);
void rf_pairwise_distance(Node *T1, vector<Node *> &gene_trees);
void rf_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end);
void rf_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees);
void rf_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end);
int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int threshold);
int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int threshold, vector<int> *original_scores);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, Forest **out_F1, Forest **out_F2);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map, int min_k, int max_k);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, int min_k, int max_k);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map, int min_k, int max_k, Forest **out_F1, Forest **out_F2);
void reduction_leaf_mult(Forest *T1, Forest *T2);
void reduction_leaf(Forest *T1, Forest *T2);
void reduction_leaf(Forest *T1, Forest *T2, UndoMachine *um);
bool chain_match(Node *T1_node, Node *T2_node, Node *T2_node_end);
Node *find_subtree_of_approx_distance(Node *n, Forest *F1, Forest *F2, int target_size);
Node *find_best_root(Node *T1, Node *T2);
double find_best_root_acc(Node *T1, Node *T2);
void find_best_root_hlpr(Node *T2, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc);
void find_best_root_hlpr(Node *n, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc,
		int *p_group_1_descendants, int *p_group_2_descendants, int *num_ties);
int rf_distance(Node *T1, Node *T2);
int count_differing_bipartitions(Node *n);
bool contains_bipartition(Node *n, int pre_start, int pre_end,
		int group_1_total, int group_2_total, int *p_group_1_descendants,
		int *p_group_2_descendants);
void modify_bipartition_support(Node *T1, Node *T2, enum RELAXATION relaxed);
void modify_bipartition_support(Node *n, Forest *F1, Forest *F2,
		Node *T1, Node *T2, vector<int> *F1_descendant_counts, enum RELAXATION);
void modify_bipartition_support(Forest *F1, Forest *F2, Node *n1);
bool is_nonbranching(Forest *T1, Forest *T2, Node *T1_a, Node *T1_c, Node *T2_a, Node *T2_c);
bool outgroup_root(Node *T, set<string, StringCompare> outgroup);
bool outgroup_root(Node *n, vector<int> &num_in, vector<int> &num_out);
bool outgroup_reroot(Node *n, vector<int> &num_in, vector<int> &num_out);
void count_in_out(Node *n, vector<int> &num_in, vector<int> &num_out,
		set<string, StringCompare> &outgroup);
void randomize_tree_with_spr(Node* T1, Node* T2, int count);
/*Joel's part*/
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map);
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2);
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2, bool verbose);
int rSPR_total_distance(Forest *T1, vector<Node *> &gene_trees);

bool BB = false;
bool APPROX_CHECK_COMPONENT = false;
bool APPROX_REVERSE_CUT_ONE_B = false;
bool APPROX_REVERSE_CUT_ONE_B_2 = false;
bool APPROX_CUT_ONE_B = false;
bool APPROX_CUT_TWO_B = false;
bool APPROX_CUT_TWO_B_ROOT = false;
bool APPROX_EDGE_PROTECTION = false;
bool CUT_ONE_B = false;
bool REVERSE_CUT_ONE_B = false;
bool REVERSE_CUT_ONE_B_2 = false;
bool REVERSE_CUT_ONE_B_3 = false;
bool CUT_TWO_B = false;
bool CUT_TWO_B_ROOT = false;
bool CUT_ALL_B = false;
bool CUT_AC_SEPARATE_COMPONENTS = false;
bool CUT_ONE_AB = false;
bool CLUSTER_REDUCTION = false;
bool PREFER_RHO = false;
bool MAIN_CALL = true;
bool MEMOIZE = false;
bool MULTIFURCATING = false;
bool MULT_4_BRANCH = false;
bool USE_CASE_7 = true;
bool ALL_MAFS = false;
int NUM_CLUSTERS = 0;
int MAX_CLUSTERS = -1;
bool UNROOTED_MIN_APPROX = false;
bool VERBOSE = false;
bool CLAMP = false;
int MAX_SPR = 1000;
int CLUSTER_MAX_SPR = MAX_SPR;
int MIN_SPR = 0;
bool FIND_RATE = false;
bool EDGE_PROTECTION = false;
bool EDGE_PROTECTION_TWO_B = false;
bool ABORT_AT_FIRST_SOLUTION = false;
bool PREORDER_SIBLING_PAIRS = false;
bool DEEPEST_ORDER = false;
bool DEEPEST_PROTECTED_ORDER = false;
bool NEAR_PREORDER_SIBLING_PAIRS = false;
bool LEAF_REDUCTION = false;
bool LEAF_REDUCTION2 = false;
bool SPLIT_APPROX = false;
bool IN_SPLIT_APPROX = false;
int SPLIT_APPROX_THRESHOLD = 25;
float INITIAL_TREE_FRACTION = 0.4;
bool COUNT_LOSSES = false;
bool CUT_LOST = false;
bool CHECK_MERGE_DEPTH = false;
bool check_all_pairs = true;
bool PREFER_NONBRANCHING = false;
int CLUSTER_TUNE = -1;
int SIMPLE_UNROOTED_LEAF = 0;
bool SHOW_CLUSTERS = false;

#ifdef RR4_CAPTURE_SPLIT_REGIONS
static vector<vector<int> > RR4_CAPTURED_SPLIT_REGIONS;
#endif

#ifdef WC1_CAPTURE_WHIDDEN_COMPONENTS
struct WC1CapturedWhiddenComponent {
	vector<int> leaves;
	int cluster_id;
	int split_id;
	int source_id;
};

static vector<vector<int> > WC1_CAPTURED_WHIDDEN_COMPONENTS;
static vector<WC1CapturedWhiddenComponent> WC1_CAPTURED_WHIDDEN_RECORDS;

static bool wc1_whidden_leaf_id(const string &name,
		map<int, string> *reverse_label_map, int &leaf_id) {
	if (name.empty() || name == "p")
		return false;
	string label = name;
	bool numeric_name = true;
	for (size_t i = 0; i < name.size(); i++) {
		if (!isdigit((unsigned char) name[i])) {
			numeric_name = false;
			break;
		}
	}
	if (numeric_name && reverse_label_map != NULL) {
		map<int, string>::iterator found =
			reverse_label_map->find(atoi(name.c_str()));
		if (found != reverse_label_map->end())
			label = found->second;
	}
	if (label.empty())
		return false;
	int value = 0;
	for (size_t i = 0; i < label.size(); i++) {
		if (!isdigit((unsigned char) label[i]))
			return false;
		value = value * 10 + (label[i] - '0');
	}
	if (value <= 0)
		return false;
	leaf_id = value;
	return true;
}

static int wc1_capture_source_id(const char *source) {
	if (source == NULL) return 0;
	if (strcmp(source, "cluster_greedy") == 0) return 1;
	if (strcmp(source, "cluster_exact") == 0) return 2;
	if (strcmp(source, "split_commit") == 0) return 3;
	return 0;
}

static void wc1_capture_whidden_forest(Forest *forest,
		map<int, string> *reverse_label_map, const char *source,
		int cluster_id, int split_id) {
	if (forest == NULL)
		return;
	int added = 0;
	for (int i = 0; i < forest->num_components(); i++) {
		Node *component = forest->get_component(i);
		if (component == NULL)
			continue;
		vector<Node *> leaves = component->find_leaves();
		vector<int> ids;
		ids.reserve(leaves.size());
		bool ok = true;
		for (size_t j = 0; j < leaves.size(); j++) {
			int leaf_id = 0;
			if (!wc1_whidden_leaf_id(leaves[j]->get_name(),
					reverse_label_map, leaf_id)) {
				ok = false;
				break;
			}
			ids.push_back(leaf_id);
		}
		if (!ok || ids.empty())
			continue;
		sort(ids.begin(), ids.end());
		ids.erase(unique(ids.begin(), ids.end()), ids.end());
		WC1_CAPTURED_WHIDDEN_COMPONENTS.push_back(ids);
		WC1CapturedWhiddenComponent record;
		record.leaves = ids;
		record.cluster_id = cluster_id;
		record.split_id = split_id;
		record.source_id = wc1_capture_source_id(source);
		WC1_CAPTURED_WHIDDEN_RECORDS.push_back(record);
		added++;
	}
	cerr << "wc1_capture source=" << source
		 << " cluster=" << cluster_id
		 << " split=" << split_id
		 << " forest_components=" << forest->num_components()
		 << " added=" << added
		 << " total=" << WC1_CAPTURED_WHIDDEN_COMPONENTS.size()
		 << endl;
}
#endif

#if defined(WC2_ALT_SPLIT_SELECT) || defined(WC3_CAPTURE_ALT_SPLIT_TRIALS) || defined(WC4_CLUSTER_ILP_SPLIT)
static void wc2_collect_split_candidates_hlpr(Node *node, Forest *F1,
		Forest *F2, int target_size, vector<pair<int, Node *> > &out) {
	if (node == NULL)
		return;
	vector<Node *> leaves = node->find_leaves();
	int leaf_count = (int) leaves.size();
	if (leaf_count > 1 && leaf_count <= target_size * 3) {
		int greedy = h25_greedy_subtree_score(node, F1, F2);
		if (greedy >= 0) {
			int target_distance = max(1, SPLIT_APPROX_THRESHOLD);
			int score = 100000
				- 2048 * abs(greedy - target_distance)
				- 8 * abs(leaf_count - target_size)
				+ greedy;
			out.push_back(make_pair(score, node));
		}
	}
	list<Node *> children = node->get_children();
	for (list<Node *>::iterator child = children.begin();
			child != children.end(); child++) {
		wc2_collect_split_candidates_hlpr(
				*child, F1, F2, target_size, out);
	}
}
#endif

#ifdef WC2_ALT_SPLIT_SELECT
static Node *wc2_select_split_candidate(Node *fallback, Forest *F1,
		Forest *F2, int target_size, int cluster_id, int split_id) {
	if (F1 == NULL || F1->num_components() == 0 ||
			F1->get_component(0) == NULL)
		return fallback;
	vector<pair<int, Node *> > candidates;
	wc2_collect_split_candidates_hlpr(
			F1->get_component(0), F1, F2, target_size, candidates);
	if (candidates.empty())
		return fallback;
	sort(candidates.begin(), candidates.end(),
			[](const pair<int, Node *> &a, const pair<int, Node *> &b) {
				if (a.first != b.first) return a.first > b.first;
				int as = (int) a.second->find_leaves().size();
				int bs = (int) b.second->find_leaves().size();
				return as < bs;
			});
	Node *chosen = candidates.front().second;
	cerr << "wc2_split_select cluster=" << cluster_id
		 << " split=" << split_id
		 << " candidates=" << candidates.size()
		 << " fallback_pre="
		 << (fallback == NULL ? -1 : fallback->get_preorder_number())
		 << " chosen_pre=" << chosen->get_preorder_number()
		 << " chosen_leaves=" << chosen->find_leaves().size()
		 << " score=" << candidates.front().first
		 << endl;
	return chosen;
}
#endif

#ifdef WC3_CAPTURE_ALT_SPLIT_TRIALS
#ifndef WC3_ALT_SPLIT_LIMIT_VALUE
#define WC3_ALT_SPLIT_LIMIT_VALUE 4
#endif
#ifndef WC3_ALT_SPLIT_K_WINDOW_VALUE
#define WC3_ALT_SPLIT_K_WINDOW_VALUE 8
#endif

static void wc3_capture_alt_split_trials(Node *normal_choice,
		Forest *base_f1, Forest *base_f2, int target_size,
		map<int, string> *reverse_label_map, int cluster_id, int split_id) {
	if (base_f1 == NULL || base_f2 == NULL || base_f1->num_components() == 0)
		return;
	vector<pair<int, Node *> > candidates;
	wc2_collect_split_candidates_hlpr(base_f1->get_component(0), base_f1,
			base_f2, target_size, candidates);
	if (candidates.empty())
		return;
	sort(candidates.begin(), candidates.end(),
			[](const pair<int, Node *> &a, const pair<int, Node *> &b) {
				if (a.first != b.first) return a.first > b.first;
				int as = (int) a.second->find_leaves().size();
				int bs = (int) b.second->find_leaves().size();
				return as < bs;
			});

	int normal_pre = normal_choice == NULL
		? -1 : normal_choice->get_preorder_number();
	int tried = 0;
	int captured = 0;
	for (size_t candidate_index = 0;
			candidate_index < candidates.size() &&
			tried < WC3_ALT_SPLIT_LIMIT_VALUE;
			candidate_index++) {
		Node *candidate = candidates[candidate_index].second;
		if (candidate == NULL)
			continue;
		int candidate_pre = candidate->get_preorder_number();
		if (candidate_pre == normal_pre)
			continue;
		tried++;

		Forest f1a = Forest(*base_f1);
		Forest f2a = Forest(*base_f2);
		Node *a_split_node = f1a.find_by_prenum(candidate_pre);
		if (a_split_node == NULL)
			continue;
		f1a.get_component(0)->disallow_siblings_subtree();
		a_split_node->allow_siblings_subtree();
		int start = rSPR_worse_3_approx(a_split_node, &f1a, &f2a);
		if (start == INT_MAX)
			start = 0;
		start /= 3;
		int end = min(base_f1->get_component(0)->size(),
				start + WC3_ALT_SPLIT_K_WINDOW_VALUE);
		bool found = false;
		for (int k = start; k <= end && !found; k++) {
			Forest f1s = Forest(*base_f1);
			Forest f2s = Forest(*base_f2);
			if (!sync_twins(&f1s, &f2s))
				break;
			Node *split_node = f1s.find_by_prenum(candidate_pre);
			if (split_node == NULL)
				break;
			f1s.get_component(0)->disallow_siblings_subtree();
			split_node->allow_siblings_subtree();
			set<SiblingPair > *sibling_pairs =
				find_sibling_pairs_set(split_node);
			list<Node *> singletons = f2s.find_singletons();
			list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();
			list<Node *> protected_stack = list<Node *>();
			int num_ties = 2;
			int split_k = rSPR_branch_and_bound_hlpr(&f1s, &f2s, k,
					sibling_pairs, &singletons, false, &AFs,
					&protected_stack, &num_ties);
			delete sibling_pairs;
			if (!AFs.empty()) {
				wc1_capture_whidden_forest(&AFs.front().first,
						reverse_label_map, "split_alt",
						cluster_id, split_id);
				captured++;
				found = true;
				cerr << "wc3_alt_split cluster=" << cluster_id
					 << " split=" << split_id
					 << " pre=" << candidate_pre
					 << " k=" << k
					 << " split_k=" << split_k
					 << " leaves=" << candidate->find_leaves().size()
					 << " captured=" << captured
					 << endl;
			}
		}
	}
	cerr << "wc3_alt_split_done cluster=" << cluster_id
		 << " split=" << split_id
		 << " tried=" << tried
		 << " captured=" << captured
		 << endl;
}
#endif

#ifdef WC4_CLUSTER_ILP_SPLIT
#ifndef WC4_SPLIT_TRIAL_LIMIT_VALUE
#define WC4_SPLIT_TRIAL_LIMIT_VALUE 6
#endif
#ifndef WC4_SPLIT_K_WINDOW_VALUE
#define WC4_SPLIT_K_WINDOW_VALUE 10
#endif
#ifndef WC4_SPLIT_FUTURE_GATE_VALUE
#define WC4_SPLIT_FUTURE_GATE_VALUE 0
#endif
#ifndef WC4_SPLIT_FUTURE_GATE_MAX_STEPS_VALUE
#define WC4_SPLIT_FUTURE_GATE_MAX_STEPS_VALUE 8
#endif

static bool wc4_collect_split_solution_by_pre(
		Forest *base_f1, Forest *base_f2, int split_pre,
		bool full_window, list<pair<Forest,Forest> > &out,
		int *out_k, int *out_split_k) {
	if (base_f1 == NULL || base_f2 == NULL || base_f1->num_components() == 0)
		return false;
	Forest f1a = Forest(*base_f1);
	Forest f2a = Forest(*base_f2);
	Node *a_split_node = f1a.find_by_prenum(split_pre);
	if (a_split_node == NULL)
		return false;
	f1a.get_component(0)->disallow_siblings_subtree();
	a_split_node->allow_siblings_subtree();
	int start = rSPR_worse_3_approx(a_split_node, &f1a, &f2a);
	if (start == INT_MAX)
		start = 0;
	start /= 3;
	int end = base_f1->get_component(0)->size();
	if (!full_window)
		end = min(end, start + WC4_SPLIT_K_WINDOW_VALUE);

	for (int k = start; k <= end; k++) {
		Forest f1s = Forest(*base_f1);
		Forest f2s = Forest(*base_f2);
		if (!sync_twins(&f1s, &f2s))
			return false;
		Node *split_node = f1s.find_by_prenum(split_pre);
		if (split_node == NULL)
			return false;
		f1s.get_component(0)->disallow_siblings_subtree();
		split_node->allow_siblings_subtree();
		set<SiblingPair > *sibling_pairs =
			find_sibling_pairs_set(split_node);
		list<Node *> singletons = f2s.find_singletons();
		list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();
		list<Node *> protected_stack = list<Node *>();
		int num_ties = 2;
		int split_k = rSPR_branch_and_bound_hlpr(&f1s, &f2s, k,
				sibling_pairs, &singletons, false, &AFs,
				&protected_stack, &num_ties);
		delete sibling_pairs;
		if (!AFs.empty()) {
			out.push_back(AFs.front());
			if (out_k != NULL) *out_k = k;
			if (out_split_k != NULL) *out_split_k = split_k;
			return true;
		}
	}
	return false;
}

static int wc4_boundary_greedy_distance(Forest *first, Forest *second) {
	if (first == NULL || second == NULL ||
			first->num_components() <= 0 ||
			second->num_components() <= 0 ||
			first->get_component(0) == NULL ||
			second->get_component(0) == NULL)
		return INT_MAX / 4;
	Forest first_boundary(first->get_component(0));
	Forest second_boundary(second->get_component(0));
	if (!sync_twins(&first_boundary, &second_boundary))
		return INT_MAX / 4;
	int distance = h25_local_greedy_distance(&first_boundary, &second_boundary);
	if (distance >= 0)
		return distance;
	return max((int) first_boundary.get_component(0)->find_leaves().size(),
			(int) second_boundary.get_component(0)->find_leaves().size());
}

static int wc4_split_commit_score(
		Forest *candidate_first, Forest *candidate_second,
		int old_components) {
	if (candidate_first == NULL || candidate_second == NULL)
		return INT_MAX / 4;
	int committed_delta =
		max(0, candidate_first->num_components() - old_components);
	int boundary_distance =
		wc4_boundary_greedy_distance(candidate_first, candidate_second);
	return committed_delta + boundary_distance;
}

static int wc4_project_split_finish_score(Forest *start_first,
		Forest *start_second) {
	uint32_t saved_rand_state = g_portable_rand_state;
	auto finish_score = [&](int score) {
		g_portable_rand_state = saved_rand_state;
		return score;
	};
	if (start_first == NULL || start_second == NULL)
		return finish_score(INT_MAX / 4);
	Forest f1 = Forest(*start_first);
	Forest f2 = Forest(*start_second);
	if (!sync_twins(&f1, &f2))
		return finish_score(INT_MAX / 4);
	int steps = 0;
	while (SPLIT_APPROX && steps < WC4_SPLIT_FUTURE_GATE_MAX_STEPS_VALUE) {
		Node *original_split_node = find_subtree_of_approx_distance(
				f1.get_component(0), &f1, &f2, SPLIT_APPROX_THRESHOLD * 2);
		if (original_split_node == NULL)
			break;
		if (original_split_node == f1.get_component(0) && steps > 0)
			break;
		Forest f1a = Forest(f1);
		Forest f2a = Forest(f2);
		Node *a_split_node =
			f1a.find_by_prenum(original_split_node->get_preorder_number());
		if (a_split_node == NULL)
			break;
		f1a.get_component(0)->disallow_siblings_subtree();
		a_split_node->allow_siblings_subtree();
		int start = rSPR_worse_3_approx(a_split_node, &f1a, &f2a);
		if (start == INT_MAX)
			start = 0;
		start /= 3;
		int end = f1.get_component(0)->size();
		bool committed = false;
		for (int k = start; k <= end; k++) {
			Forest f1s = Forest(f1);
			Forest f2s = Forest(f2);
			if (!sync_twins(&f1s, &f2s))
				break;
			Node *split_node =
				f1s.find_by_prenum(original_split_node->get_preorder_number());
			if (split_node == NULL)
				break;
			f1s.get_component(0)->disallow_siblings_subtree();
			split_node->allow_siblings_subtree();
			set<SiblingPair > *sibling_pairs =
				find_sibling_pairs_set(split_node);
			list<Node *> singletons = f2s.find_singletons();
			list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();
			list<Node *> protected_stack = list<Node *>();
			int num_ties = 2;
			rSPR_branch_and_bound_hlpr(&f1s, &f2s, k,
					sibling_pairs, &singletons, false, &AFs,
					&protected_stack, &num_ties);
			delete sibling_pairs;
			if (!AFs.empty()) {
				AFs.front().first.swap(&f1);
				AFs.front().second.swap(&f2);
				f2.unprotect_edges();
				f1.get_component(0)->allow_siblings_subtree();
				committed = true;
				break;
			}
		}
		if (!committed)
			break;
		steps++;
	}
	return finish_score(
			f1.num_components() + wc4_boundary_greedy_distance(&f1, &f2));
}

static bool wc4_try_split_ilp_commit(Node *normal_choice,
		Forest *f1, Forest *f2, map<int, string> *reverse_label_map,
		int cluster_id, int split_id, int &total_split_k) {
	if (normal_choice == NULL || f1 == NULL || f2 == NULL)
		return false;

	vector<pair<int, Node *> > ranked;
	wc2_collect_split_candidates_hlpr(
			f1->get_component(0), f1, f2,
			SPLIT_APPROX_THRESHOLD * 2, ranked);
	sort(ranked.begin(), ranked.end(),
			[](const pair<int, Node *> &a, const pair<int, Node *> &b) {
				if (a.first != b.first) return a.first > b.first;
				int as = (int) a.second->find_leaves().size();
				int bs = (int) b.second->find_leaves().size();
				return as < bs;
			});

	vector<int> split_prenums;
	split_prenums.push_back(normal_choice->get_preorder_number());
	for (size_t idx = 0;
			idx < ranked.size() &&
			(int)split_prenums.size() < WC4_SPLIT_TRIAL_LIMIT_VALUE;
			idx++) {
		if (ranked[idx].second == NULL)
			continue;
		int pre = ranked[idx].second->get_preorder_number();
		if (find(split_prenums.begin(), split_prenums.end(), pre) ==
				split_prenums.end()) {
			split_prenums.push_back(pre);
		}
	}

	list<pair<Forest,Forest> > split_candidates;
	int normal_k = 0;
	int normal_split_k = 0;
	int solved = 0;
	for (int idx = 0; idx < (int)split_prenums.size(); idx++) {
		int trial_k = 0;
		int trial_split_k = 0;
		bool full_window = idx == 0;
		if (wc4_collect_split_solution_by_pre(
				f1, f2, split_prenums[idx], full_window,
				split_candidates, &trial_k, &trial_split_k)) {
			solved++;
			if (idx == 0) {
				normal_k = trial_k;
				normal_split_k = trial_split_k;
			}
		}
	}
	if (split_candidates.empty())
		return false;

	Forest *chosen_first = NULL;
	Forest *chosen_second = NULL;
	bool ilp_ok = wc4_try_split_ilp(
			f1, f2, &split_candidates, reverse_label_map,
			cluster_id, split_id, &chosen_first, &chosen_second);
	int old_components = f1->num_components();
	int normal_score = wc4_split_commit_score(
			&split_candidates.front().first,
			&split_candidates.front().second,
			old_components);
	int ilp_score = ilp_ok && chosen_first != NULL && chosen_second != NULL
		? wc4_split_commit_score(chosen_first, chosen_second, old_components)
		: INT_MAX / 4;
	if (ilp_ok && ilp_score > normal_score) {
		cerr << "wc4_split_ilp_gate_reject cluster=" << cluster_id
			 << " split=" << split_id
			 << " normal_score=" << normal_score
			 << " ilp_score=" << ilp_score
			 << " normal_components="
			 << split_candidates.front().first.num_components()
			 << " ilp_components=" << chosen_first->num_components()
			 << endl;
		delete chosen_first;
		delete chosen_second;
		chosen_first = NULL;
		chosen_second = NULL;
		ilp_ok = false;
	}
	if (WC4_SPLIT_FUTURE_GATE_VALUE &&
			ilp_ok && chosen_first != NULL && chosen_second != NULL) {
		int normal_future_score = wc4_project_split_finish_score(
				&split_candidates.front().first,
				&split_candidates.front().second);
		int ilp_future_score = wc4_project_split_finish_score(
				chosen_first, chosen_second);
		if (ilp_future_score > normal_future_score) {
			cerr << "wc4_split_ilp_future_gate_reject cluster=" << cluster_id
				 << " split=" << split_id
				 << " normal_future=" << normal_future_score
				 << " ilp_future=" << ilp_future_score
				 << " normal_components="
				 << split_candidates.front().first.num_components()
				 << " ilp_components=" << chosen_first->num_components()
				 << endl;
			delete chosen_first;
			delete chosen_second;
			chosen_first = NULL;
			chosen_second = NULL;
			ilp_ok = false;
		}
		else {
			cerr << "wc4_split_ilp_future_gate_accept cluster=" << cluster_id
				 << " split=" << split_id
				 << " normal_future=" << normal_future_score
				 << " ilp_future=" << ilp_future_score
				 << endl;
		}
	}
	if (!ilp_ok || chosen_first == NULL || chosen_second == NULL) {
		chosen_first = new Forest(split_candidates.front().first);
		chosen_second = new Forest(split_candidates.front().second);
	}

	int committed_components = chosen_first->num_components();
	chosen_first->swap(f1);
	chosen_second->swap(f2);
	delete chosen_first;
	delete chosen_second;

	f2->unprotect_edges();
	f1->get_component(0)->allow_siblings_subtree();
	if (ilp_ok)
		total_split_k += max(0, committed_components - old_components);
	else if (normal_k > 0 || normal_split_k > 0)
		total_split_k += max(0, normal_k - normal_split_k);
	else
		total_split_k += max(0, f1->num_components() - old_components);

#ifdef WC1_CAPTURE_WHIDDEN_COMPONENTS
	wc1_capture_whidden_forest(f1, reverse_label_map,
			ilp_ok ? "split_ilp" : "split_commit",
			cluster_id, split_id);
#endif
	cerr << "wc4_split_ilp cluster=" << cluster_id
		 << " split=" << split_id
		 << " trials=" << split_prenums.size()
		 << " solved=" << solved
		 << " ilp_ok=" << (ilp_ok ? 1 : 0)
		 << " components=" << f1->num_components()
		 << " total_split_k=" << total_split_k
		 << endl;
	return true;
}
#endif

class ProblemSolution {
public:
string T1;
string T2;
int k;

ProblemSolution(Forest *t1, Forest *t2, int new_k) {
	T1 = t1->str();
	T2 = t2->str();
	k = new_k;
}
	};

	map<string, ProblemSolution> memoized_clusters = map<string, ProblemSolution>();

/*******************************************************************************
	RSPR WORSE_3_MULT_APPROX
*******************************************************************************/

/* rSPR_worse_3_mult_approx
 * Calculate an approximate maximum agreement forest and SPR distance for two multifurcating trees
 * RETURN At most 3 times the rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 * T1 and T2 can be  multifurcating forests.
 */
int rSPR_worse_3_mult_approx(Forest *T1, Forest *T2) {
	return rSPR_worse_3_mult_approx(T1, T2, true);
}

int rSPR_worse_3_mult_approx(Forest *T1, Forest *T2, bool sync) {
	// match up nodes of T1 and T2
	if (sync) {
if (!sync_twins(T1, T2))
	return 0;
	}

	if (LEAF_REDUCTION) {
	  reduction_leaf_mult(T1, T2);
	}
//	cout << "T1: "; T1->print_components();
//	cout << "T2: "; T2->print_components();
	// find sibling pairs of T1
	list<Node *> *sibling_groups = T1->find_sibling_groups();


	/*
	list<Node *>::iterator c;
	list<Node *>::iterator ch;
	
	for (c = sibling_groups->begin(); c != sibling_groups->end(); c++) {
	  cout << "Sibling group parent: " << (*c)->str() << endl;
		list<Node *> children = (*c)->get_children();
		for (ch = children.begin(); ch != children.end(); ch++) {
		  cout << "\tChild: " << (*ch)->str() <<endl;
		}
		list<list<Node *>> *identical_sibling_groups = (*c)->find_identical_sibling_groups();
		if (identical_sibling_groups->size() > 0) {
		  cout << "Printing identical sibling group:" << endl;
		  for (auto g = identical_sibling_groups->begin(); g != identical_sibling_groups->end(); g++) {
			cout << "Sibling identical group:" << endl;
			for (auto ge = g->begin(); ge != g->end(); ge++) {
			  cout << "\tchild:" << (*ge)->get_name() << endl;
			}			  
		  }
		}
		else {		  
		  cout << "No identical sibling groups" <<endl;
		}
	}	
	//return 0;
	*/
	
	
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	//list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	Forest *F1;
	Forest *F2;

	T2->max_preorder = T2->components[0]->get_max_preorder_number(0);//preorder_number(0);
	int ans = rSPR_worse_3_mult_approx_hlpr(T1, T2, &singletons, sibling_groups, &F1, &F2, true);

	F1->swap(T1);
	F2->swap(T2);
	sync_twins(T1,T2);


	delete sibling_groups;
	delete F1;
	delete F2;
	return ans;
}


// rSPR_worse_3_mult_approx recursive helper function
int rSPR_worse_3_mult_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_groups, Forest **F1, Forest **F2, bool save_forests) {

  int num_cut = 0;
  Node* previous_group = NULL;
  while(!singletons->empty() || !sibling_groups->empty()) {
	  
    // Case 1 - Remove singletons
    while(!singletons->empty()) {
      
      Node *T2_a = singletons->back();
      singletons->pop_back();
      #ifdef DEBUG_APPROX
      cout << "Handling singleton: " << T2_a->str() << endl;
      #endif
      // find twin in T1
      Node *T1_a = T2_a->get_twin();
      // if this is in the first component of T_2 then
      // it is not really a singleton.
      // TODO: problem when we cluster and have a singleton as the
      //		first comp of T2
      //    NEED TO MODIFY CUTTING?
      // 		HERE AND IN BB?
      Node *T1_a_p = T1_a->parent();
      if (T1_a_p == NULL)
	continue;

      
      if (T2_a == T2->get_component(0)){
	if (!T1->contains_rho()) {
	  T1->add_rho();
	  T2->add_rho();
	  num_cut++;
	}
      }
	
	//continue;

      bool is_sibling_group = T1_a_p->is_sibling_group();
      // cut the edge above T1_a
      T1_a->cut_parent();      
      if (!T1_a->is_leaf()) {
	T1_a_p->decrement_non_leaf_children();//although would this ever be a non leaf?
      }
      T1->add_component(T1_a);

      //only contract if one node
      if (T1_a_p->get_children().size() == 1) {
	if (is_sibling_group) {
	  sibling_groups->remove(T1_a_p);
          #ifdef DEBUG_APPROX	  
	  cout << "Removed ";
	  for (list<Node*>::iterator i = T1_a_p->get_children().begin(); i != T1_a_p->get_children().end(); i++) {
	    cout << (*i)->str();
	  }
	  cout  << " from sibling groups" << endl;
	  #endif
	}
	Node *possible_previous_sibling = T1_a_p->get_children().front();
	bool was_sibling_group = possible_previous_sibling->is_sibling_group();
	Node *node = T1_a_p->contract(true);
	if (node != NULL) {
	  node->recalculate_non_leaf_children(); //can we tell what this would be instead of recalculating?

	  if (node->is_sibling_group()) {
	    if (was_sibling_group) {
	      list<Node*>::iterator i = find(sibling_groups->begin(), sibling_groups->end(), possible_previous_sibling);
	      *i = node;
	    }
	    else{
	      sibling_groups->push_front(node);
	    }
	  }
	}
      }
    }//!singletons->empty()

    
    if(!sibling_groups->empty()) {
      //Get the first group that has identical sibling groups, otherwise default to the group on the back
      list<Node*>::iterator i = sibling_groups->end();
      i--;
      Node *T1_sibling_group = sibling_groups->back();
      list<list<Node*>> identical_sibling_groups = list<list<Node*>>();
      T1_sibling_group->find_identical_sibling_groups(&identical_sibling_groups);
      for (; i != sibling_groups->begin(); i-- ){
	(*i)->find_identical_sibling_groups(&identical_sibling_groups);
	if (identical_sibling_groups.size() > 0) {
	  T1_sibling_group = (*i);
	  break;
	}
      }
      
      #ifdef DEBUG_APPROX
      cout << "F2: ";
      T2->print_components();
      cout << endl;
      cout << "F1: ";
      T1->print_components();
      cout << endl;
      #endif
      
      /* 
	 Case where a subset of the group have the same parent both in T1 and T2
	 Step 5 in paper
      */
      // Case 2 - Contract identical sibling pair
      if (identical_sibling_groups.size() > 0) {	  		
	list<list<Node *>>::iterator i;
	for (i = identical_sibling_groups.begin(); i != identical_sibling_groups.end(); i++) {
	  list<Node *> T2_group = (*i);
	  Node *T2_p = T2_group.front()->parent();
	  #ifdef DEBUG_APPROX
	  cout << "Contracting T1... " << endl;
	  #endif
	  Node *T1_group_new = T1_sibling_group->contract_twin_group(&T2_group);
	  #ifdef DEBUG_APPROX
	  cout << "Contracting T2... " << endl;
	  #endif
	  Node *T2_group_new = T2_p->contract_sibling_group(&T2_group);
	  
	  T1_group_new->set_twin(T2_group_new);
	  T2_group_new->set_twin(T1_group_new);			

	  // check if T2_p is a singleton after the contraction
	  if (T2_p->is_singleton() && T2_p != T2->get_component(0)) {
	    //(T2_p->is_singleton() && T1_sibling_group != T1->get_component(0) && T2_p != T2->get_component(0)) {
	    singletons->push_front(T2_p);
	  }
	  if (T1_sibling_group->parent() != NULL) {
	    //Check if the contraction made a new sibling group
	    T1_sibling_group->parent()->recalculate_non_leaf_children();
	    if (T1_sibling_group->parent()->is_sibling_group()) {
	      sibling_groups->push_front(T1_sibling_group->parent());
	    }
	  }
	  if (!T1_sibling_group->is_sibling_group()) {
	    sibling_groups->remove(T1_sibling_group);
	  }	  
	}
      }

      /*
	4 branching case
	Step 6-8 in paper
	Cut above a1, b1, a2, b2
	Part 1: Get the LCA, this is the numbering to the top part
	Part 2: Get subset of sibling group that is descendant of this LCA
	Part 3: Sort them based on depth
	Part 4: Since this is the approximation, we cut deepest 2
      */

      // Case 3
      else {
	#ifdef DEBUG_APPROX
	cout << "Sibling group to be cutting: " << T1_sibling_group->str_subtree() << endl;
	#endif
	//vector of ints describing how many of the siblings are in its descendants, indexed by preorder number
	vector<int> descendant_count = T1_sibling_group->find_pseudo_lca_descendant_count(T2->max_preorder + 1);
	Node* arbitrary_lca = T1_sibling_group->find_arbitrary_lca(T2->components, descendant_count);
	vector<Node *> deepest_siblings;
	
	//All siblings are in different components, ie no path between them
	//Get depth of siblings from root of each component
	if (arbitrary_lca == NULL) {	
	  vector<vector<Node *>> siblings_by_depth = vector<vector<Node *>>(10);
	  for (int i = 0; i != T2->components.size(); i++) {
	    T2->components[i]->get_deepest_siblings(descendant_count, siblings_by_depth);
	  }
	  deepest_siblings = contract_deepest_siblings(siblings_by_depth);
	}
	//Otherwise they share an LCA
	else {
	  vector<vector<Node *>> siblings_by_depth = arbitrary_lca->get_deepest_siblings(descendant_count);
	  deepest_siblings = contract_deepest_siblings(siblings_by_depth);
	}
      
	// Should assert here
	if (deepest_siblings.size() < 2) { cout << "improper length" << endl; }

	//Get the deepest two of the siblings
	Node *T2_a1 = deepest_siblings[0];
	Node *T2_a2 = deepest_siblings[1];
	#ifdef DEBUG_APPROX
	cout << "a1: " << T2_a1->str() << " a2: " << T2_a2->str() << endl;
	#endif
	
	bool cut_a1   = false;
	bool cut_a1_p = false;
	bool cut_a2   = false;
	bool cut_a2_p = false;
      
	if (T1_sibling_group->get_children().size() == 2) {
	  /*
	    7.1 case
	    Cut a1, pa1, a2 in F2, add 3 to num_cut
	    Consider adding 3 regardless if we actually cut 3, 
	  */	
	  cut_a1   = true;
	  cut_a1_p = true;
	  cut_a2   = true;
	  //num_cut += 3;
	  /*
	    7.2 case
	    Cut a1, a2, pa1, pa2, add 4 to num_cuts
	  */
	  if (previous_group == T1_sibling_group) {
	    cut_a2_p = true;
	    #ifdef DEBUG_APPROX
	    cout <<"Case 7.2" << endl;
	    #endif
	    //num_cut += 1;
	  }
	  else {
	    #if DEBUG_APPROX
	    cout << "Case 7.1" << endl;
	    #endif
	  }
	} // size == 2
      
	else if (T1_sibling_group->get_children().size() > 2) {
	  /*
	    7.3 case
	    If a2's parent's only sibling is part of the sibling group, 
	    and a1's parent is a root or has a sibling that is not part of the sibling group
	    then cut a2 and a2_p otherwise a1 and a1_p
	  */
	  if (previous_group != T1_sibling_group) {
	    Node* T2_a2_p = T2_a2->parent();
	    bool x_2 = false;
	    bool a2_p_one_sibling = (T2_a2_p != NULL) &&
	      (T2_a2_p->parent() != NULL) &&
	      (T2_a2_p->parent()->get_children().size() == 2);	  
	    if (a2_p_one_sibling) {
	      list<Node *> group = T1_sibling_group->get_children();
	      //get the other one
	      Node *a2_p_sibling = T2_a2_p->parent()->get_children().front() == T2_a2_p ?
		T2_a2_p->parent()->get_children().back() :
		T2_a2_p->parent()->get_children().front();	
	      //check if it is part of sibling group
	      bool a2_p_sibling_in_group = descendant_count[a2_p_sibling->get_preorder_number()] == -1;
	      if (a2_p_sibling_in_group) {
		Node* T2_a1_p = T2_a1->parent();
		bool a1_p_is_root = T2_a1_p->parent() == NULL;
		bool a1_p_sibling_not_in_group = false;
		if (!a1_p_is_root) {
		  list<Node*> a1_p_siblings = T2_a1_p->parent()->get_children();
		  for (list<Node*>::iterator i = a1_p_siblings.begin(); i != a1_p_siblings.end(); i++) {
		    if (descendant_count[(*i)->get_preorder_number()] != -1) {
		      a1_p_sibling_not_in_group = true;
		      break;
		    }
		  }
		}
		if (a2_p_one_sibling && a2_p_sibling_in_group && (a1_p_is_root || a1_p_sibling_not_in_group)) {
		  x_2 = true;

		}
	      }
	    }
	    #ifdef DEBUG_APPROX
	    cout << "Case 7.3" << endl;
	    #endif
	    //num_cut += 2;
	    if (x_2){
	    
	      cut_a2   = true;
	      cut_a2_p = true;
	    }
	    else {
	      cut_a1   = true;
	      cut_a1_p = true;
	    }
	  }
	  /*7.4 case
	    cut a1 and a1_p
	  */
	  else if (previous_group == T1_sibling_group) {
	    #ifdef DEBUG_APPROX
	    cout << "Case 7.4" << endl;
	    #endif
	    cut_a1   = true;
	    cut_a1_p = true;
	    //num_cut += 2;
	  }
	}
	/*
	  Cutting section
	*/
	Node* T2_a2_p = NULL;
	if (cut_a1) {
	  Node *T2_a1_p = T2_a1->parent();
	  //singletons? components?
	  if (T2_a1_p != NULL) {
	    //Cut connections
	    T2_a1->cut_parent();
	    num_cut++;
	    //add as components
	    T2->add_component(T2_a1);
	    //Just cut 1 of two children of a1 parent
	    if (T2_a1_p->get_children().size() == 1) {
	      if (T2_a1_p->parent() == NULL) {
		T2_a1_p->contract(true);
		if (T2_a1_p->is_singleton() && T2_a1_p != T2->get_component(0)) {
		  singletons->push_front(T2_a1_p);
		}
	      }
	      else {
		Node* T2_b1 = T2_a1_p->get_children().front();
		T2_a1_p->contract(true);
		T2_a1_p = T2_b1; // <------------------------------------
	      }
	    }
	    //Check for singletons
	    if (T2_a1->is_singleton()) // wont ever be C0?
	      singletons->push_front(T2_a1);
	    

	    if (cut_a1_p) {
	      if (T2_a1_p->parent() != NULL) {
		bool aborted_a2 = false;
		if (T2_a1_p->parent() == T2_a2->parent() &&
		    T2_a2->parent()->get_children().size() == 2)
		  {
		    cut_a2 = false; //cutting a1_p will cause a2 to get contracted up, so there is no more a2 to cut
		    cut_a2_p = true; //Instead we cut a2_p
		    aborted_a2 = true;
		  }
		Node *T2_a1_gp = T2_a1_p->parent();
		//Cut connections
		T2_a1_p->cut_parent();
		num_cut++;
		//add as components
		T2->add_component(T2_a1_p);
		//Just cut 1 of two children of a1 parent
		if (T2_a1_gp->get_children().size() == 1) {
		  if (T2_a1_gp->parent() == NULL) {
		    T2_a1_gp->contract(true);		    
		  }
		  else {
		    Node* T2_b1_p = T2_a1_gp->get_children().front();
		    T2_a1_gp->contract(true);
		    T2_a1_gp = T2_b1_p;
		  }
		  if (aborted_a2) { T2_a2_p = T2_a1_gp; }
		  if (T2_a1_gp != NULL && T2_a1_gp->is_singleton() && (T2_a1_gp != T2->get_component(0) || aborted_a2)) {
		    singletons->push_front(T2_a1_gp);
		  }
		}
		//Check for singletons
		if (T2_a1_p->is_singleton() && T2_a1_p != T2->get_component(0)) {
		  singletons->push_front(T2_a1_p);	      
		}	
	      }
	    }
	  }//T2_a1_p() != NULL
	}
	if (cut_a2){
	  T2_a2_p = T2_a2->parent();
	  //could have cut a2's parent in previous steps, so could be singleton now
	  if (T2_a2->is_leaf() && T2_a2_p == NULL) {
	    singletons->push_front(T2_a2);
	  }
	  else if (T2_a2_p != NULL) {
	    //Cut connections
	    T2_a2->cut_parent();
	    num_cut++;
	    //add as components
	    T2->add_component(T2_a2);
	    //Just cut 1 of two children of a2 parent
	    if (T2_a2_p->get_children().size() == 1) {
	      if (T2_a2_p->parent() == NULL) {
		T2_a2_p->contract(true);
		if (T2_a2_p->is_singleton() && T2_a2_p != T2->get_component(0)) {
		  singletons->push_front(T2_a2_p);
		}
	      }
	      else {
		Node* T2_b2 = T2_a2_p->get_children().front();
		T2_a2_p->contract(true);
		T2_a2_p = T2_b2;
	      }
	    }
	    //Check for singletons
	    if (T2_a2->is_singleton())
	      singletons->push_front(T2_a2);
	    
	  }
	}
	if (cut_a2_p) {
	  if (T2_a2_p != NULL && T2_a2_p->parent() != NULL) {
	    Node *T2_a2_gp = T2_a2_p->parent();
	    //Cut connections
	    T2_a2_p->cut_parent();
	    	    num_cut++;
	    //add as components
	    T2->add_component(T2_a2_p);
	    //Just cut 1 of two children of a2 parent
	    if (T2_a2_gp->get_children().size() == 1) {
	      if (T2_a2_gp->parent() == NULL) {
		T2_a2_gp->contract(true);
	      }
	      else {
		Node* T2_b2_p = T2_a2_gp->get_children().front();
		T2_a2_gp->contract(true);
		T2_a2_gp = T2_b2_p;
	      }

	      if (T2_a2_gp->is_singleton() && T2_a2_gp != T2->get_component(0)) {
		singletons->push_front(T2_a2_gp);
	      }
	    }
	    //Check for singletons
	    if (T2_a2_p->is_singleton()) {
	      singletons->push_front(T2_a2_p);	      
	    }

	  }	 	 
	}
      

      }//else
      //delete identical_sibling_groups;
      previous_group = T1_sibling_group;
    }//!sibling_groups->empty()

  } //while(!sibling_groups->empty() && !singletons->empty()
  // if the first component of the forests differ then we have cut p
  /*
  if (T1->get_component(0)->get_twin() != T2->get_component(0)) {
    if (!T1->contains_rho()) {
      T1->add_rho();
      T2->add_rho();
    }
    else
      // hack to ignore rho when it shouldn't be in a cluster
      num_cut -=3;
  }
  */
  if (save_forests) {
    *F1 = new Forest(T1);
    *F2 = new Forest(T2);
  }
  //if (num_cut) num_cut--;
  return num_cut;
}

#ifdef DEBUG_CASE_COUNTER
struct mult_case_counter {
  int case_71, case_72, case_73, case_74;
  int case_81, case_82, case_83, case_84;
  int case_85, case_86, case_87;
};
static mult_case_counter case_counter = {};
void print_mult_case_count() {
  cout << "Stats:" << endl;
  cout << "\tCase 7.1:" << case_counter.case_71 << endl;
  cout << "\tCase 7.2:" << case_counter.case_72 << endl;
  cout << "\tCase 7.3:" << case_counter.case_73 << endl;
  cout << "\tCase 7.4:" << case_counter.case_74 << endl;
  cout << "\tCase 8.1:" << case_counter.case_81 << endl;
  cout << "\tCase 8.2:" << case_counter.case_82 << endl;
  cout << "\tCase 8.3:" << case_counter.case_83 << endl;
  cout << "\tCase 8.4:" << case_counter.case_84 << endl;
  cout << "\tCase 8.5:" << case_counter.case_85 << endl;
  cout << "\tCase 8.6:" << case_counter.case_86 << endl;
  cout << "\tCase 8.7:" << case_counter.case_87 << endl;
}
  
#endif

int rSPR_branch_and_bound_mult_range(Forest *T1, Forest *T2, int start_k){
  return rSPR_branch_and_bound_mult_range(T1, T2, start_k, MAX_SPR);
}
int rSPR_branch_and_bound_mult_range(Forest *T1, Forest *T2, int start_k, int end_k){
  int exact_spr = -1;
  int k;
  for (k = start_k; k <= end_k; k++) {
    Forest *F1 = new Forest(T1);
    Forest *F2 = new Forest(T2);
    if (!sync_twins(F1,F2)) {
      exact_spr = 0;
      continue;
    }
    if (LEAF_REDUCTION) {
      reduction_leaf_mult(F1, F2);
    }
    #ifdef DEBUG
    cout << "Trying K = " << k << endl << "------------------" << endl;
    #else
    cout << k << " " << endl;
    #endif
    exact_spr = rSPR_branch_and_bound_mult(F1, F2, k);
    #ifdef DEBUG
    cout << "------------------" << endl;
    cout << "Finished K = " << k << " return value : " << exact_spr << endl;   
    #endif
    if (exact_spr >= 0) {
      F1->swap(T1);
      F2->swap(T2);
    }
    delete F1;
    delete F2;

    if (exact_spr >= 0) {    
      break;
    }
  }
  #ifdef DEBUG_CASE_COUNTER
  print_mult_case_count();
  #endif
  if (k > end_k) {
    k = -1;
  }
  return k;
}
int rSPR_branch_and_bound_mult(Forest *T1, Forest *T2, int k){

  if (!sync_twins(T1,T2)) {
    return 0;      
  }
  T2->max_preorder = T2->components[0]->get_max_preorder_number(0);//preorder_number(0);          
  list<Node *> *sibling_groups = T1->find_sibling_groups();
  //sibling_groups->push_front(sibling_groups->back());
  //sibling_groups->pop_back();
  list<Node *> singletons     = T1->find_singletons();
  
  list<pair<Forest,Forest>> AFs = list<pair<Forest,Forest>>();
  //list<Node *> protected_stack = list<Node*>();
  int num_ties = 2;
  int final_k = rSPR_branch_and_bound_mult_hlpr(T1, T2, k, sibling_groups, &singletons, NULL, &AFs, &num_ties);  

  //print AFs
  if (!AFs.empty() && final_k > -1) {
    if (ALL_MAFS) {
      cout << endl << endl << "FOUND ANSWER" << endl;
      // TODO: this is a cheap hack
      for (list<pair<Forest,Forest> >::iterator x = AFs.begin(); x != AFs.end(); x++) {
	cout << "\tT1: ";
	x->first.print_components();
	cout << "\tT2: ";
	x->second.print_components();
      }
    }
      AFs.front().first.swap(T1);
  AFs.front().second.swap(T2);
  sync_twins(T1,T2);

  }

  
  delete sibling_groups;
  if (final_k >= 0)
  return k - final_k;
  else
    return final_k;
    
}


class bb_mult_recurse_data {
public:
  Forest *T1;
  Forest *T2;
  list<Node*> *sibling_groups;
  list<Node*> *singletons;
  map<Node *, Node*> node_map;
};

/* Generates a copy of the data used to recurse on rSPR_branch_and_bound_mult_hlpr so 
   that original forests aren't clobbered,
   updates pointers accordingly */
//TODO: 1 new not 5
bb_mult_recurse_data *generate_mult_recurse_data(Forest *T1, Forest *T2, list<Node*> *sibling_groups, list<Node*> *singletons) {

  bb_mult_recurse_data *data = new bb_mult_recurse_data();

  data->node_map = map<Node*, Node*>();
  data->T1 = new Forest(T1, &data->node_map);
  data->T2 = new Forest(T2, &data->node_map);
  //TODO: smarter way of syncing
  //for (auto n = node_map.begin(); n != node_map.end(); n++) {
    //(*n).first->set_twin(node_map[(*n).first->get_twin()]);
    /*
      if ((*n).second->is_leaf()) {
      (*n).second->set_twin(node_map[(*n).first->get_twin()]);
      node_map[(*n).first->get_twin()]->set_twin((*n).second);
      cout << "Synced twin : " << (*n).second->str() << " -> " << (*n).second->get_twin()->str() << endl;}*/
  //}
  sync_twins(data->T1, data->T2);
  data->sibling_groups = new list<Node*>();
  for (auto n = sibling_groups->begin(); n != sibling_groups->end(); n++) {
    data->sibling_groups->push_back(data->node_map[*n]);
  }
  data->singletons = new list<Node*>();
  for (auto n = singletons->begin(); n != singletons->end(); n++) {
    data->singletons->push_back(data->node_map[*n]);
  }
  return data;
}
//Cuts to_cut, adds to components, conditionally adds to singletons
//Assumes parent is not null and no Null parameters
void mult_cut_and_cleanup(Node* to_cut, Forest *T2, list<Node*> *singletons) {
  //Node* T2_a1_next = next_data->node_map[T2_a1];
  Node* to_cut_p = to_cut->parent();
  //Cut connections  
  to_cut->cut_parent();
  //add as components
  T2->add_component(to_cut);
  //Just cut 1 of two children of a1 parent
  if (to_cut_p->get_children().size() == 1) {
    if (to_cut_p->parent() == NULL) {
      to_cut_p->contract(true);
      if (to_cut_p->is_singleton() && to_cut_p != T2->get_component(0)) {
      singletons->push_front(to_cut_p);
      }
    }
    else {
      Node* to_cut_b = to_cut_p->get_children().front();
      to_cut_p->contract(true);
      if (to_cut_b->is_singleton() && to_cut_b != T2->get_component(0)) {
	singletons->push_front(to_cut_b);
      }
    }
  }
  //Check for singletons
  if (to_cut->is_singleton())
    singletons->push_front(to_cut);  
}
//cuts everything except node, possibly expanding, adds to components, conditionally adds to singletons
//Assumes parent is not null and no Null parameters
//NOTE: potentially unsafe for preorder numbers
void mult_cut_all_except_and_cleanup(Node* T2_a1, Forest *T2, list<Node*> *singletons) {
  Node* T2_b1;
  Node* parent = T2_a1->parent();
  if (parent->get_children().size() == 2) {
    T2_b1 = parent->get_children().front() == T2_a1 ?
      parent->get_children().back() :
      parent->get_children().front();		
  }
  else {
    list<Node*> all_but_a1 = list<Node*>(parent->get_children());
    all_but_a1.remove(T2_a1);	    
    T2_b1 = parent->expand_children_out(all_but_a1);
    T2_b1->set_preorder_number(parent->get_preorder_number());
  }
  //hack for now
  //We immediately contract a1 up so we know the parent's preorder number is available
  //Cut connections

  T2_b1->cut_parent();

  //add as components
  T2->add_component(T2_b1);
	    
  //Just cut 1 of two children of a2 parent (will always be in this case?)
  if (parent->get_children().size() == 1) {
    if (parent->parent() == NULL) {
      parent->contract(true);
      if (parent->is_singleton() && parent != T2->get_component(0)) {
	singletons->push_front(parent);
      }
    }
    else {
      parent = parent->contract(true);
      if (T2_a1->is_singleton() && T2_a1 != T2->get_component(0))
	singletons->push_front(T2_a1);
    }
  }
  if (T2_b1->is_singleton())
    singletons->push_front(T2_b1);
}

//Adds rho and a certain singleton to cut. Basically its for realising when we have cut rho but it is already a singleton so we explicitly continue on 
#define MULT_RHO_CUT_AND_RESOLVE(rho_to_singleton, node_to_protect) {			\
  bb_mult_recurse_data *next_data = generate_mult_recurse_data(T1, T2, sibling_groups, singletons); \
  Node* protect = NULL;							\
  next_data->T1->add_rho();						\
  next_data->T2->add_rho();						\
  next_data->singletons->push_front(next_data->node_map[rho_to_singleton]);	\
  int result_k = rSPR_branch_and_bound_mult_hlpr(next_data->T1, next_data->T2, k - 1, next_data->sibling_groups, next_data->singletons, protect, AFs, num_ties); \
  delete next_data->T1;							\
  delete next_data->T2;						\
  delete next_data->sibling_groups;				\
  delete next_data->singletons;					\
  delete next_data;						\
  if (result_k > best_k) {					\
    best_k = result_k;						\
  }								\
  }								\

//TODO (Ben): try using a routine instead of a macro, I suspect it will slow down because of
//            stack and parameter passing though
#define MULT_BB_CUT_AND_RESOLVE(nodes_to_cut, nodes_to_exclude_cutting, node_to_protect) { \
  bb_mult_recurse_data *next_data = generate_mult_recurse_data(T1, T2, sibling_groups, singletons);\
  int num_cuts = 0;							\
  for (int i = 0; i < nodes_to_cut.size(); i++) {			\
    Node* T2_ax_next = next_data->node_map[nodes_to_cut[i]];		\
    mult_cut_and_cleanup(T2_ax_next, next_data->T2, next_data->singletons);\
    num_cuts++;								\
  }									\
  for (int i = 0; i < nodes_to_exclude_cutting.size(); i++) {		\
    Node* T2_ax_next = next_data->node_map[nodes_to_exclude_cutting[i]]; \
    mult_cut_all_except_and_cleanup(T2_ax_next, next_data->T2, next_data->singletons); \
    num_cuts++;								\
  }									\
  Node* protect = next_data->node_map[node_to_protect];			\
  int result_k = rSPR_branch_and_bound_mult_hlpr(next_data->T1, next_data->T2, k - num_cuts, next_data->sibling_groups, next_data->singletons, protect, AFs, num_ties); \
  delete next_data->T1;							\
  delete next_data->T2;							\
  delete next_data->sibling_groups;					\
  delete next_data->singletons;						\
  delete next_data;							\
  if (result_k > best_k) {						\
    best_k = result_k;							\
    if (!ALL_MAFS && best_k > -1) {					\
    }									\
  }									\
}

//TODO: UndoMachine, then cleanup all constructors, bb_mult_recurse_data relying on copies of trees
int rSPR_branch_and_bound_mult_hlpr(Forest *T1, Forest *T2,
				    int k,
				    list<Node*> *sibling_groups, list<Node*> *singletons,
				    Node *protected_node, list<pair<Forest,Forest>> *AFs,
				    int* num_ties) {
  //run out of cuts if k < 0
  if (k < 0) {
    return k;
  }
  Node* previous_group = sibling_groups->back();
  int best_k = -1;
  while(!singletons->empty() || !sibling_groups->empty()) {
	  
    // Case 1 - Remove singletons
    while(!singletons->empty()) {
      
      Node *T2_a = singletons->back();
      singletons->pop_back();
      #ifdef DEBUG
      cout << "Handling singleton: " << T2_a->str() << endl;
      #endif
      
      Node *T1_a = T2_a->get_twin();
      Node *T1_a_p = T1_a->parent();
      
      if (T1_a_p == NULL)
	continue;      

      // find twin in T1
      //If we have added component 0, then we have made a B cut and added rho
      if (T2_a == T2->get_component(0)){// && T1_a != T1->get_component(0)) {
	if (!T1->contains_rho()) {
	  T1->add_rho();
	  T2->add_rho();
	  k--;
	  //continue;
	}
      }

      bool is_sibling_group = T1_a_p->is_sibling_group();
      // cut the edge above T1_a
      T1_a->cut_parent();
      if (!T1_a->is_leaf()) {
	T1_a_p->decrement_non_leaf_children();
      }
      T1->add_component(T1_a);
      
      //only contract if one node
      if (T1_a_p->get_children().size() == 1) {
	//If we contract this node, and it used to be a sibling group
	//then it is not a sibling group anymore
	if (is_sibling_group) {
	  sibling_groups->remove(T1_a_p);
          #ifdef DEBUG
	  cout << "Removed " << T1_a_p->str() << " from sibling groups" << endl;
	  #endif
	}
	
	Node *possible_previous_sibling = T1_a_p->get_children().front();
	bool was_sibling_group = possible_previous_sibling->is_sibling_group();
	
	Node *T1_a_gp = T1_a_p->parent();
	Node *T1_new_a_p = T1_a_p;
	
	//If the child is a sibling group, there is a possibility contract()
	//will delete it, so we need to update it in the sibling groups
	if (T1_a_p->parent() == NULL) {
	  if (was_sibling_group){
	    list<Node*>::iterator i = find(sibling_groups->begin(), sibling_groups->end(), possible_previous_sibling);
	    *i = T1_new_a_p;	    
	  }
	}
	else {
	  T1_new_a_p = possible_previous_sibling;
	}
       
        T1_a_p->contract(true);
	T1_new_a_p->recalculate_non_leaf_children();
	//After contracting this, the grandparent may be a sibling group now
	if (T1_a_gp != NULL) {
	  T1_a_gp->recalculate_non_leaf_children(); //can we tell what this would be instead of recalculating?
	  if (T1_a_gp->is_sibling_group()) {	 
	    sibling_groups->push_front(T1_a_gp);	 
	  }
	}
      }

    }//!singletons->empty()

    //NOTE: we know there are no singletons left here
    if(!sibling_groups->empty()) {
      //Find identical sibling groups
      //Get the first group that has identical sibling groups, otherwise default to the group on the back
      list<Node*>::reverse_iterator i = sibling_groups->rbegin();
      Node *T1_sibling_group = sibling_groups->back();
      list<list<Node*>> identical_sibling_groups;
      T1_sibling_group->find_identical_sibling_groups(&identical_sibling_groups);
      bool found_identical = false;
      for (; i != sibling_groups->rend(); i++ ){
	(*i)->find_identical_sibling_groups(&identical_sibling_groups);
	if (identical_sibling_groups.size() > 0) {
	  T1_sibling_group = (*i);
	  found_identical = true;
	  break;
	}
      }
      #ifdef MULT_PICK_LARGEST_GROUP
      if (!found_identical) {
	int max_size = 0;
	Node* largest_group = NULL;
	for (auto i = sibling_groups->begin(); i != sibling_groups->end(); i++) {
	  if ((*i)->get_children().size() > max_size) {
	    max_size = (*i)->get_children().size();
	    largest_group = (*i);
	  }
	}
	T1_sibling_group = largest_group;
      }
      #endif
      #ifdef DEBUG
      cout << "K = " << k << endl;
      cout << "F2: ";
      T2->print_components();
      cout << endl;
      cout << "F1: ";
      T1->print_components();
      cout << endl;
      #endif

      /* 
	 Case where a subset of the group have the same parent both in T1 and T2
      */
      // Case 2 - Contract identical sibling pair
      if (identical_sibling_groups.size() > 0) {	  		
	list<list<Node *>>::iterator i;
	for (i = identical_sibling_groups.begin(); i != identical_sibling_groups.end(); i++) {
	  //Contract the groups
	  list<Node *> T2_group = (*i);
	  Node *T2_p = T2_group.front()->parent();
	  #ifdef DEBUG
	  //cout << "Contracting T1... " << endl;
	  #endif
	  Node *T1_group_new = T1_sibling_group->contract_twin_group(&T2_group);
	  #ifdef DEBUG
	  //cout << "Contracting T2... " << endl;
	  #endif
	  Node *T2_group_new = T2_p->contract_sibling_group(&T2_group);

	  //Maintain twins
	  T1_group_new->set_twin(T2_group_new);
	  T2_group_new->set_twin(T1_group_new);			

	  // check if T2_p is a singleton after the contraction
	  if (T2_p->is_singleton() && T2_p != T2->components[0]){
	    singletons->push_front(T2_p);
	  }
	  if (T1_sibling_group->parent() != NULL) {
	    //Check if the contraction made a new sibling group
	    T1_sibling_group->parent()->recalculate_non_leaf_children();
	    if (T1_sibling_group->parent()->is_sibling_group()) {
	      sibling_groups->push_front(T1_sibling_group->parent());
	      #ifdef DEBUG
	      //cout << "Added new sibling group after contraction: " << T1_sibling_group->parent()->str_subtree() << endl;
	      #endif
	    }
	  }
	  if (!T1_sibling_group->is_sibling_group()) {
	    sibling_groups->remove(T1_sibling_group);
	      #ifdef DEBUG
	    //cout << "Removed new sibling group after contraction: " << T1_sibling_group->str() << endl;
	      #endif

	  }	  
	}
      }

      /*
	4 branching case
	Step 6-8 in paper
	Cut above a1, b1, a2, b2
	Part 1: Get the LCA, this is the numbering to the top part
	Part 2: Get subset of sibling group that is descendant of this LCA
	Part 3: Sort them based on depth
	Part 4: Since this is the approximation, we cut deepest 2
      */

      // Case 3
      else {

	//Check branch and bound
	if (BB) {
	  //copies for the approx so we dont clobber this tree
	  map<Node*, Node*> approx_map = map<Node*, Node*>();
	  Forest T1_approx = Forest(T1, &approx_map);
	  Forest T2_approx = Forest(T2, &approx_map);
	  sync_twins(&T1_approx, &T2_approx);
	  list<Node*> sibling_group_approx = list<Node*>();
	  for (auto n = sibling_groups->begin(); n != sibling_groups->end(); n++) {
	    sibling_group_approx.push_back(approx_map[*n]);
	  }
	  list<Node*> singletons_approx = list<Node*>();
	  for (auto n = singletons->begin(); n != singletons->end(); n++) {
	    singletons_approx.push_back(approx_map[*n]);
	  }
	  int approx_spr = rSPR_worse_3_mult_approx_hlpr(&T1_approx, &T2_approx, &singletons_approx, &sibling_group_approx, NULL, NULL, false);
	  //TODO: Sometimes approx returns 1 over the right amount. For example 4 for a 1 cut tree or 16 for a 3 cut tree
	  //If approx_spr > 3k then we will not have enough cuts
	  if (approx_spr > 3*k + 1) {
#ifdef DEBUG
	    cout << "approx failed approx k = " << approx_spr  <<  endl;
#endif
	    return -1;
	  }
	}

	//Finding deepest siblings
	#ifdef DEBUG
	cout << "Sibling group to be cutting: " << T1_sibling_group->str_subtree() << endl;
	#endif
	//vector of ints describing how many of the siblings are in its descendants, indexed by preorder number
	vector<int> descendant_count = T1_sibling_group->find_pseudo_lca_descendant_count(T2->max_preorder + 1);
	Node* arbitrary_lca = T1_sibling_group->find_arbitrary_lca(T2->components, descendant_count);
	vector<Node *> deepest_siblings;
	vector<vector<Node*>> siblings_by_depth;
	//Get depth of siblings from root of each component
	//If the lca is null, all siblings are in different components, ie no path between them
	if (arbitrary_lca == NULL) {	
	  siblings_by_depth = vector<vector<Node *>>(10);
	  for (int i = 0; i != T2->components.size(); i++) {
	    //if preorder is -1 then this is rho
	    if (T2->components[i]->get_preorder_number() == -1) { continue; }
	    T2->components[i]->get_deepest_siblings(descendant_count, siblings_by_depth);
	  }
	}
	//Otherwise they share an LCA
	else {
	  siblings_by_depth = arbitrary_lca->get_deepest_siblings(descendant_count);
	}
	map<Node*, int> s_map = map<Node*, int>();
	//Get deepest siblings returns sparse vector, this compacts it
	deepest_siblings = contract_deepest_siblings(siblings_by_depth, &s_map);      
	// Should assert here
	if (deepest_siblings.size() < 2) { cout << "improper length" << endl; }

	//Get the deepest two of the siblings
	Node *T2_a1 = deepest_siblings[0];
	Node *T2_a2 = deepest_siblings[1];
	#ifdef DEBUG
	cout << "a1: " << T2_a1->str() << " a2: " << T2_a2->str() << endl;
	#endif
	best_k = -1;

	//MULT_4_BRANCH is naive approach of making 4 cuts every time
	if (!MULT_4_BRANCH) {
	//step 7
	/*
	  if a1 == prot, x = 2 else x = 1
	*/
	if (protected_node != NULL) {
	  Node* T2_ax;
	  if (T2_a1 != protected_node) {
	    T2_ax = T2_a1;
	  }
	  else {
	    T2_ax = T2_a2;
	  }
	  
	  bool a0_descendant_of_lca = false;
	  Node* stepper = protected_node;
	  while (stepper->parent() != NULL && stepper->parent() != arbitrary_lca) {
	    stepper = stepper->parent();
	  }
	  if (stepper->parent() == arbitrary_lca) {
	    a0_descendant_of_lca = true;
	  }

	  //TODO (Ben) do one iteration of pl parent for has_aj_child and has_non_aj_child,
	  //           as well as any other flags, before step 7 check
	  bool pl_has_aj_child = false;
	  if (arbitrary_lca != NULL) {
	    Node* pl = arbitrary_lca->parent();	  
	    if (pl != NULL) {
	      for (auto i = pl->get_children().begin(); i != pl->get_children().end(); i++) {
		if (*i != arbitrary_lca && descendant_count[(*i)->get_preorder_number()] == -1) {
		  pl_has_aj_child = true;
		  break;
		}
	      }
	    }
	  }
	  //step 7.1
	  if (arbitrary_lca == NULL) {
	    #ifdef DEBUG
	    cout << "Case 7.1 T2_ax: " << T2_ax->str() << endl;
	    #endif
	    #ifdef DEBUG_CASE_COUNTER
	    case_counter.case_71++;
	    #endif
	    /* 
	       recurse on cut ax prot protected node
	    */
	    //7.1 : cut ax
	    if (T2_ax->parent() != NULL) {
	      vector<Node*> to_cut = {T2_ax};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	    }
	  }

	  //step 7.2
	  else if (!a0_descendant_of_lca) {
	    /*
	      recurse on cut a1's B's up to LCA prot protected node

	    */
	    #ifdef DEBUG
	    cout << "Case 7.2a cut all B1's Protected Node: " << protected_node->str() <<  endl;
	    #endif
	    #ifdef DEBUG_CASE_COUNTER
	    case_counter.case_72++;
	    #endif

	    {
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a1;
		
	      while(stepper->parent() != arbitrary_lca) { 
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	    }
	    
	    if (s_map[T2_a1] > 1 ||
		s_map[deepest_siblings[deepest_siblings.size()-1]] > 0){
	      /*
		recurse on cut a1 prot protected node
	      */
#ifdef DEBUG
	    cout << "Case 7.2b Cut a1 Protected Node: " << protected_node->str() <<  endl;
#endif
	      //7.2b cut a1
	      if (T2_a1->parent() != NULL) {
		vector<Node*> to_cut = {T2_a1};
		vector<Node*> to_cut_except = {};
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	      }     
	    }
	     /*

	      Not part of the outlined special cases

	     */
	    //if (T2_a1->parent()->get_children().size() == 2)
	    {
	      #ifdef DEBUG
	      cout << "Case 7.2c Cut a2 T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	      #endif
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);	      
	    }
	  }	
	    

	  //TODO: Is component 0 considered a root?
	  //If we do get component 0 can we consider this if rho has been added?
	  //step 7.3
	  else if (deepest_siblings.size() == 2 &&
		   T2_ax->parent() == arbitrary_lca &&
		   a0_descendant_of_lca &&
		   (
		    (arbitrary_lca->parent() == NULL && arbitrary_lca != T2->get_component(0)) ||
		    pl_has_aj_child
		    )) {
	    /*
	      recurse on cut Bx prot protected node
	    */
	    //7.3 cut Bx

	    #ifdef DEBUG
	    cout << "Case 7.3 Cut T2_bx T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	    #endif
	    #ifdef DEBUG_CASE_COUNTER
	    case_counter.case_73++;
	    #endif

	    {
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_ax};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	    }
	    /*

	      This is NOT in the cases. Figure out whats happening with 8.2 with r = 2

	     */
	    /*
	    {
	      cout << "Case 7.3 Cut T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	      vector<Node*> to_cut = {T2_ax};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	      }
	    */
	  }


	  //step 7.4
	  else if (a0_descendant_of_lca) {
	    /*
	      recurse on cut ax prot protected_node
	                 if m > 2 and r > 2 recurse on 
			    cut all Bx's then B`x
	     */
	    //7.4 cut ax
	    #ifdef DEBUG
	    cout << "Case 7.4 T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	    #endif
	    #ifdef DEBUG_CASE_COUNTER
	    case_counter.case_74++;
	    #endif

	    {
	      vector<Node*> to_cut = {T2_ax};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	    }
	    if (T1_sibling_group->get_children().size() > 2 && deepest_siblings.size() > 2) {
	      #ifdef DEBUG
	      cout << "Case 7.4b T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	      #endif
	      //7.4 cut all Bx B`x
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_ax;
		
	      while(stepper->parent() != arbitrary_lca) { //Do we know ax is descendant of lca?
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());
	      to_cut_except.push_back(T2_ax);
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);	      
	    }
	    /*

	      Not part of the outlined special cases

	     */
	    {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);	     
	    }
	    {
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);	     
	    }

	    /*
	    if (T1_sibling_group->get_children().size() == 2 &&
		deepest_siblings.size() == 2 &&
		T2_a1 == T2_ax &&
		T2_a2 != protected_node) {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	      }*/
	  }
	}

	//step 8
	else {
	  bool all_but_ar_s1 = true;
	  for (int i = 0; i < deepest_siblings.size() - 1; i++) {
	    if (s_map[deepest_siblings[i]] != 1) {
	      all_but_ar_s1 = false;
	      break;
	    }
	  }
	  bool all_s1 = s_map[deepest_siblings[deepest_siblings.size()-1]] == 1 &&
	    all_but_ar_s1;

	  bool lca_p_contains_sibling = false;	    
	  if (arbitrary_lca != NULL && arbitrary_lca->parent() != NULL) {
	    for (list<Node*>::iterator i = arbitrary_lca->parent()->get_children().begin();
		 i != arbitrary_lca->parent()->get_children().end();
		 i++) {
	      if (descendant_count[(*i)->get_preorder_number()] == -1) {
		lca_p_contains_sibling = true;
		break;
	      }
	    }
	  }

	  //step 8.1
	  if (arbitrary_lca == NULL) {
	    /*
	      recurse on cut a1 no prot,
	      cut a2 no prot
	    */
#ifdef DEBUG
	    cout << "Case 8.1a cut a1" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_81++;
#endif

	    //8.1 : Cut a1
	    if (T2_a1->parent() != NULL) {
	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    else if (T2_a1 == T2->get_component(0)) {
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(T2_a1, NULL);
	      }
	    }
#ifdef DEBUG
	    cout << "Case 8.1b cut a2" << endl;
#endif
	    //8.1 : Cut a2
	    if (T2_a2->parent() != NULL) {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    else if (T2_a2 == T2->get_component(0)) {
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(T2_a2, NULL);
	      }
	    }
	  
	  }
	  //step 8.2
	  else if (all_but_ar_s1 &&
		   deepest_siblings[deepest_siblings.size()-1]->parent() == arbitrary_lca) {
	    /*
	      recurse on cut all B's except shallowest
	      for each sibling except for shallowest,
	      cut all other B's protect ai
		           
	    */
	    //8.2 B's
	    {
#ifdef DEBUG
	    cout << "Case 8.2a cut B1 - B(r-1)" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_82++;
#endif

	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      for (int i = 0; i < deepest_siblings.size() - 1; i++) {
		to_cut_except.push_back(deepest_siblings[i]);
	      }
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //all other B's except ai
	    //Doesnt explicitly say this in the paper, but
	    //this would only happen if r > 2


	    if (deepest_siblings.size() > 2){
	      for (int i = 0; i < deepest_siblings.size() - 1; i++) {
#ifdef DEBUG
	      cout << "Case 8.2b for all ai cut all other B" << endl;
#endif

		vector<Node*> to_cut = {};
		vector<Node*> to_cut_except = {};
		for (int j = 0; j < deepest_siblings.size() - 1; j++) {
		  if (i != j) {
		    to_cut_except.push_back(deepest_siblings[j]);
		  }
		}
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, deepest_siblings[i]);
	      }
	      }
	    /*
	    else if (T1_sibling_group->parent() != NULL && T1_sibling_group->parent()->get_children().size() == 2) {
	      Node* gp = T1_sibling_group->parent();
	      Node* aunt = gp->get_children().front() == T1_sibling_group ?
		gp->get_children().back() :
		gp->get_children().front();
	      if (aunt->is_leaf()) {
		Node* a1_p = T2_a1->parent();
		vector<Node*> a1_p_leaves = a1_p->find_leaves();
		Node* T2_aunt = aunt->get_twin();
		for (int i = 0; i < a1_p_leaves.size(); i++) {
		  if (a1_p_leaves[i] == T2_aunt) {
#ifdef DEBUG
		    cout << "Case 8.2c cut a2" << endl;
#endif
		    vector<Node*> to_cut = {T2_a2};
		    vector<Node*> to_cut_except = {};
		    MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
		    break;
		  }
		}

	      }
	      }*/


	      else if (deepest_siblings.size() == 2) {
		  vector<Node*> to_cut = {T2_a2};
		  vector<Node*> to_cut_except = {};
		  MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	      }

	  }
	  

	  //step 8.3
	  else if (T1_sibling_group->get_children().size() == 2 &&
		   s_map[T2_a1] + s_map[T2_a2] >= 2) {
	    /*recurse on cut a1 no prot,
	      cut a2 no prot,
	      cut all along a1 to a2 no prot
	    */
	    //8.3 : Cut a1
	    if (T2_a1->parent() != NULL) {
#ifdef DEBUG
	      cout << "Case 8.3a cut a1" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_83++;
#endif

	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.3 : Cut a2
	    if (T2_a2->parent() != NULL) {
#ifdef DEBUG
	    cout << "Case 8.3b cut a2" << endl;
#endif
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.3 : All along path from a1, a2
	    {
#ifdef DEBUG
	    cout << "Case 8.3c cut all B1's and B2's" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a1;
		
	      while(stepper->parent() != arbitrary_lca) { //no need to check for null! we know theres a path
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      stepper = T2_a2;
	      while(stepper->parent() != arbitrary_lca) { //no need to check for null! we know theres a path
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      /*
		Cut from top to bottom,
		Cutting from bottom to top causes the contraction to invalidate
		the node, since contraction is implemented by cutting the parent, then giving
		the child to the parent
	      */
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }	      
	  }

	  //step 8.4
	  else if (T1_sibling_group->get_children().size() > 2 &&
		   deepest_siblings.size() == 2 &&
		   s_map[T2_a1] == 1 &&
		   s_map[T2_a2] == 1) {
	    /* 
	       recurse on cut a1 and a2 no prot
	       cut b1 and b2 no prot
	       cut all except b1 off of lca, b1
	       cut all except b2 off of lca, b2
	    */
	    //8.4 : Cut a1 and a2
	    {
#ifdef DEBUG
	    cout << "Case 8.4a cut a1 & a2" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_84++;
#endif

	      vector<Node*> to_cut = {T2_a1, T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.4 : Cut a1 and a2's B's
	    {
#ifdef DEBUG
	    cout << "Case 8.4b cut B1 & B2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1, T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }

	    //8.4 : All except a1->p, then b1
	    {
#ifdef DEBUG
	    cout << "Case 8.4c cut B1 & B`1" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1->parent(), T2_a1};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a2);
	    }
	    //8.4 : All except a2->p, then b2
	    {
#ifdef DEBUG
	    cout << "Case 8.4d cut B2 & B`2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a2->parent(), T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a1);
	    }
	    Node* pl = arbitrary_lca->parent();
	    bool pl_has_non_aj_child = false;
	    if (pl != NULL) {
	      for (auto i = pl->get_children().begin(); i != pl->get_children().end(); i++) {
		if (*i != arbitrary_lca && descendant_count[(*i)->get_preorder_number()] != -1) {
		  pl_has_non_aj_child = true;
		  break;
		}
	      }
	    }
	    bool gpl_has_non_aj_child = false;
	    //place this conditional so we do not unnecessarily iterate
	    //If pl_has_non_aj_child then the if will evaluate to true later on anyways
	    if (!pl_has_non_aj_child) {
	      if (pl != NULL) {
		Node* gpl = pl->parent();		    		    
		if (gpl != NULL) {
		  for (auto i = gpl->get_children().begin(); i != gpl->get_children().end(); i++) {
		    if (*i != pl && descendant_count[(*i)->get_preorder_number()] != -1) {
		      gpl_has_non_aj_child = true;
		      break;
		    }
		  }
		}
	      }
	    }
	    //8.4 special case
	    if (arbitrary_lca->parent() == NULL ||
		pl_has_non_aj_child ||
		gpl_has_non_aj_child) {
	      //8.4 Cut a1 prot a2
	      {
#ifdef DEBUG
	      cout << "Case 8.4e cut a1" << endl;
#endif
		vector<Node*> to_cut = {T2_a1};
		vector<Node*> to_cut_except = {};
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a2);	
	      }
	      //8.4 Cut a2 prot a1
	      {
#ifdef DEBUG
	      cout << "Case 8.4f a2" << endl;
#endif
		vector<Node*> to_cut = {T2_a2};
		vector<Node*> to_cut_except = {};
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a1);	
	      }

	    }
	  }

	  //step 8.5
	  else if (T1_sibling_group->get_children().size() > 2 &&
		   deepest_siblings.size() > 2 &&
		   all_s1) {
	    /*
	      recurse on cut all a1 through ar no prot,
	      cut all b1 through br no prot,
	      for each ai, cut all a's except ai prot ai
	      for each ai, cut all B's except for ai's B prot ai
	    */
	    //8.5 all ai
	    {
#ifdef DEBUG
	    cout << "Case 8.5a cut all a1-ar" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_85++;
#endif	    
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      for (int i = 0; i < deepest_siblings.size(); i++) {
		to_cut.push_back(deepest_siblings[i]);
	      }
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.5 cut all bi
	    {
#ifdef DEBUG
	      cout << "Case 8.5b cut all B1-Br" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      for (int i = 0; i < deepest_siblings.size(); i++) {
		to_cut_except.push_back(deepest_siblings[i]);
	      }
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.5 for each ai cut all a's except ai
	    {	      
	      for (int i = 0; i < deepest_siblings.size(); i++) {
#ifdef DEBUG
		cout << "Case 8.5c cut all a1-ar except ai" << endl;
#endif
		vector<Node*> to_cut = {};
		vector<Node*> to_cut_except = {};
		for (int j = 0; j < deepest_siblings.size(); j++) {
		  if (i != j) {
		    to_cut.push_back(deepest_siblings[j]);
		  }
		}
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, deepest_siblings[i]);//TODO protect ai
	      }
	    }
	    //8.5 for each ai cut all b's except ai's
	    {	      
	      for (int i = 0; i < deepest_siblings.size(); i++) {
#ifdef DEBUG
		cout << "Case 8.5d cut all B1-Br except Bi" << endl;
#endif
		vector<Node*> to_cut = {};
		vector<Node*> to_cut_except = {};
		for (int j = 0; j < deepest_siblings.size(); j++) {
		  if (i != j) {
		    to_cut_except.push_back(deepest_siblings[j]);
		  }
		}
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, deepest_siblings[i]);//TODO protect ai
	      }
	    }
	  }

	  //step 8.6
	  else if (T1_sibling_group->get_children().size() > 2 &&
		   deepest_siblings.size() == 2 &&
		   s_map[T2_a1] >= 2 &&
		   T2_a2->parent() == arbitrary_lca && //equivalent to s_map[T2_a2] == 0
		   (
		    arbitrary_lca->parent() == NULL ||
		    lca_p_contains_sibling
		    )) {
	    /*
	      recurse on cut a1 no prot,
	      cut all B's leading up to LCA from a1, no prot,
	      cut B2 (essentially a1's whole branch) no prot
	    */
	    //if (T2_a1->parent() != NULL) {
	    //8.6 cut a1
	    {
#ifdef DEBUG	     
	      cout << "Case 8.6a cut a1" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_86++;
#endif
	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.6 cut all B1s
	    {
#ifdef DEBUG	     
	      cout << "Case 8.6b cut all B1's" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a1;
		
	      while(stepper->parent() != arbitrary_lca) { 
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());

	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.6 cut B2
	    {
#ifdef DEBUG	     
	      cout << "Case 8.6c cut B2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }	      
	  }

	  //step 8.7
	  else if (T1_sibling_group->get_children().size() > 2 &&
		   s_map[T2_a1] >= 2) {
	    /*
	      recurse on cut a1 no prot,
	      cut a2 prot a1,
	      cut all B1's leading up to LCA no prot,
			   
	    */
	    //8.7 cut a1
	    {
#ifdef DEBUG	     
	      cout << "Case 8.7a cut a1" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_87++;
#endif
	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }	      
	    //8.7 cut a2
	    {
#ifdef DEBUG	     
	      cout << "Case 8.7b cut a2" << endl;
#endif
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a1);
	    }
	    //8.7 cut all B1s
	    {
#ifdef DEBUG	     
	      cout << "Case 8.7c cut all B1's" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a1;
		
	      while(stepper->parent() != arbitrary_lca) { 
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());

	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);		
	    }
	    //8.7 cut all B2s
	    {
#ifdef DEBUG
	    cout << "Case 8.7d";
	    if (deepest_siblings.size() == 2) {
	      cout << "2" << endl;
	    }
	    else {
	      cout << "n" << endl;
	    }
#endif

	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a2;
		
	      while(stepper->parent() != arbitrary_lca) { 
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());
	      //8.7 cut all B2's, cut B`2, otherwise r > 2
	      if(deepest_siblings.size() == 2) {
		to_cut_except.push_back(T2_a2); // cut B`2 at the end
	      }
	      //This exception should be caught by 8.6 case.
	      /*
		if (to_cut_except.size() == 0) {
		cout << "\n\n\nto_cut_except empty in 8.7n. Parent is lca ERROR \n\n\n"; 
		}*/

	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a1);	
	    }
	  }

	}

	}
	else { //4_MULT_BRANCH

	  //To debug: set all to true. Make all cuts all the time
	  //Otherwise uses same logic as approximation algorithm
	  bool cut_a1 = true;
	  bool cut_b1 = true;
	  bool cut_a2 = true;
	  bool cut_b2 = true;

	  if (T1_sibling_group->get_children().size() == 2) {
	    /*
	      7.1 case
	      Cut a1, pa1, a2 in F2, add 3 to num_cut
	    */	
	    cut_a1   = true;
	    cut_b1 = true;
	    cut_a2   = true;
	    /*
	      7.2 case
	      Cut a1, a2, pa1, pa2, add 4 to num_cuts
	    */
	    if (previous_group == T1_sibling_group) {
	      cut_b2 = true;
	    }
	  } // size == 2
      
	  else if (T1_sibling_group->get_children().size() > 2) {
	    /*
	      7.3 case
	      If a2's parent's only sibling is part of the sibling group, 
	      and a1's parent is a root or has a sibling that is not part of the sibling group
	      then cut a2 and a2_p otherwise a1 and a1_p
	    */
	    if (previous_group != T1_sibling_group) {
	      Node* T2_a2_p = T2_a2->parent();
	      bool x_2 = false;
	      bool a2_p_one_sibling = (T2_a2_p != NULL) &&
		(T2_a2_p->parent() != NULL) &&
		(T2_a2_p->parent()->get_children().size() == 2);	  
	      if (a2_p_one_sibling) {
		list<Node *> group = T1_sibling_group->get_children();
		//get the other one
		Node *a2_p_sibling = T2_a2_p->parent()->get_children().front() == T2_a2_p ?
		  T2_a2_p->parent()->get_children().back() :
		  T2_a2_p->parent()->get_children().front();	
		//check if it is part of sibling group
		bool a2_p_sibling_in_group = descendant_count[a2_p_sibling->get_preorder_number()] == -1;
		if (a2_p_sibling_in_group) {
		  Node* T2_a1_p = T2_a1->parent();
		  bool a1_p_is_root = T2_a1_p->parent() == NULL;
		  bool a1_p_sibling_not_in_group = false;
		  if (!a1_p_is_root) {
		    list<Node*> a1_p_siblings = T2_a1_p->parent()->get_children();
		    for (list<Node*>::iterator i = a1_p_siblings.begin(); i != a1_p_siblings.end(); i++) {
		      if (descendant_count[(*i)->get_preorder_number()] != -1) {
			a1_p_sibling_not_in_group = true;
			break;
		      }
		    }
		  }
		  if (a2_p_one_sibling && a2_p_sibling_in_group && (a1_p_is_root || a1_p_sibling_not_in_group)) {
		    x_2 = true;
		  }
		}
	      }
	      if (x_2){
		cut_a2 = true;
		cut_b2 = true;
		//cut_a1 = true;
		//cut_b1 = true;
	      }
	      else {
		cut_a1 = true;
		cut_b1 = true;
		cut_a2 = true;
		//cut_b2 = true;
	      }
	    }
	    /*7.4 case
	      cut a1 and a1_p
	    */
	    else if (previous_group == T1_sibling_group) {
	      cut_a1   = true;
	      cut_b1 = true;
	    }
	  }
	  /*
	    Cutting section
	  */
	  Node *T2_a1_p = T2_a1->parent();

	  if (cut_a1) {	  
	    if (T2_a1_p != NULL) {
#ifdef DEBUG
	      cout << "Case Cut a1" << endl;
#endif
	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //mimics case 8.1, should add rho
	    else if (T2_a1 == T2->get_component(0) && arbitrary_lca == NULL) {
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(T2_a1, NULL);
	      }
	    }
	  }
	  if (cut_b1) {
	    if (T2_a1_p != NULL) {
#ifdef DEBUG
	      cout << "Case Cut b1" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	  }
	  Node *T2_a2_p = T2_a2->parent();
	  if (cut_a2){
#ifdef DEBUG
	      cout << "Case Cut a2" << endl;
#endif
	    if (T2_a2_p != NULL) {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    else if (T2_a2 == T2->get_component(0) && arbitrary_lca == NULL) {
	      
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(T2_a2, NULL);
	      }
	    }

	  }
	  if (cut_b2) {
	    if (T2_a2_p != NULL) {
#ifdef DEBUG
	      cout << "Case Cut b2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	  }
	}
	sibling_groups->pop_back();
	return best_k;
      }//else (cutting)
      previous_group = T1_sibling_group;
    }//!sibling_groups->empty()

  } //while(!sibling_groups->empty() && !singletons->empty()

  //Made it to end, so add to results
  if (k >= 0) {
    //if (PREFER_RHO && !AFs->empty() && !AFs->front().first.contains_rho() && T1->contains_rho()) {
    if (true && !AFs->empty() && !AFs->front().first.contains_rho() && T1->contains_rho()) {
      if (!ALL_MAFS)
	AFs->clear();
      AFs->push_front(make_pair(Forest(T1),Forest(T2)));
      *num_ties = 2;
    }
    else if (ALL_MAFS || AFs->empty()) {
      AFs->push_back(make_pair(Forest(T1),Forest(T2)));
    }
    //else if (!PREFER_RHO || AFs->front().first.contains_rho() == T1->contains_rho()) {
    else if (false || AFs->front().first.contains_rho() == T1->contains_rho()) {
      if (rand() < RAND_MAX/ *num_ties) {
	AFs->clear();
	AFs->push_back(make_pair(Forest(T1),Forest(T2)));
      }
      (*num_ties)++;
    }
  }
  return k;
}


/* rSPR_3_approx
 * Calculate an approximate maximum agreement forest and SPR distance
 * RETURN At most 3 times the rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 */
int rSPR_3_approx(Forest *T1, Forest *T2) {
	// find sibling pairs of T1
	// match up nodes of T1 and T2
	if (!sync_twins(T1, T2))
return 0;
	// find singletons of T2
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	list<Node *> singletons = T2->find_singletons();
	int ans = rSPR_3_approx_hlpr(T1, T2, &singletons, sibling_pairs);
	delete sibling_pairs;
	return ans;
}

// rSPR_3_approx recursive helper function
int rSPR_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons,
list<Node *> *sibling_pairs) {
	int num_cut = 0;
	while(!singletons->empty() || !sibling_pairs->empty()) {
// Case 1 - Remove singletons
while(!singletons->empty()) {


	Node *T2_a = singletons->back();
	singletons->pop_back();
  // find twin in T1
	Node *T1_a = T2_a->get_twin();
	// if this is in the first component of T_2 then
	// it is not really a singleton.
	if (T2_a == T2->get_component(0))
		continue;

	Node *T1_a_parent = T1_a->parent();
	if (T1_a_parent == NULL)
		continue;
	bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
	// cut the edge above T1_a
	T1_a->cut_parent();
	T1->add_component(T1_a);
	if (T1_a->get_sibling_pair_status() > 0)
		T1_a->clear_sibling_pair(sibling_pairs);
	//delete(T1_a);

	Node *node = T1_a_parent->contract();
	if (node != NULL && potential_new_sibling_pair && node->is_sibling_pair()){
		node->rchild()->add_to_front_sibling_pairs(sibling_pairs, 2);
		node->lchild()->add_to_front_sibling_pairs(sibling_pairs, 1);
	}

}
if(!sibling_pairs->empty()) {
	Node *T1_a = sibling_pairs->back();
	sibling_pairs->pop_back();
	Node *T1_c = sibling_pairs->back();
	sibling_pairs->pop_back();
	T1_a->clear_sibling_pair_status();
	T1_c->clear_sibling_pair_status();
	if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
		continue;
	}
	Node *T1_ac = T1_a->parent();
	// lookup in T2 and determine the case
	Node *T2_a = T1_a->get_twin();
	Node *T2_c = T1_c->get_twin();

	// Case 2 - Contract identical sibling pair
	if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
		Node *T2_ac = T2_a->parent();
		T1_ac->contract_sibling_pair();
		T2_ac->contract_sibling_pair();
		T1_ac->set_twin(T2_ac);
		T2_ac->set_twin(T1_ac);
		T1->add_deleted_node(T1_a);
		T1->add_deleted_node(T1_c);
		T2->add_deleted_node(T2_a);
		T2->add_deleted_node(T2_c);

		// check if T2_ac is a singleton
		if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
			singletons->push_back(T2_ac);
		// check if T1_ac is part of a sibling pair
		if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
			T1_ac->parent()->lchild()->add_to_sibling_pairs(sibling_pairs, 1);
			T1_ac->parent()->rchild()->add_to_sibling_pairs(sibling_pairs, 2);
		}
	}
	// Case 3
	else {
		
		//  ensure T2_a is below T2_c
		if (T2_a->get_depth() < T2_c->get_depth()) {
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);
		}
		else if (T2_a->get_depth() == T2_c->get_depth()) {
			if (T2_a->parent() && T2_c->parent() &&
					(T2_a->parent()->get_depth() <
					T2_c->parent()->get_depth())) {
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);
			}
		}

		// get T2_b
		Node *T2_ab = T2_a->parent();
		Node *T2_b = T2_ab->rchild();
		if (T2_b == T2_a)
			T2_b = T2_ab->lchild();
		// cut T1_a, T1_c, T2_a, T2_b, T2_c

		bool cut_b_only = false;
		if (T2_a->parent() != NULL && T2_a->parent()->parent() != NULL && T2_a->parent()->parent() == T2_c->parent()) {
			cut_b_only = true;
			T1_a->add_to_sibling_pairs(sibling_pairs,1);
			T1_c->add_to_sibling_pairs(sibling_pairs,2);
		}

		if (!cut_b_only) {
			T1_a->cut_parent();
			T1_c->cut_parent();
			// contract parents
			Node *node = T1_ac->contract();
			// check for T1_ac sibling pair
			if (node != NULL && node && node->is_sibling_pair()){
				node->lchild()->add_to_sibling_pairs(sibling_pairs,1);
				node->rchild()->add_to_sibling_pairs(sibling_pairs,2);
			}
		}

		bool same_component = true;
		if (APPROX_CHECK_COMPONENT)
			same_component = (T2_a->find_root() == T2_c->find_root());

		if (!cut_b_only) {
			T2_a->cut_parent();
			num_cut++;
		}
		bool cut_b = false;
		if (same_component && T2_ab->parent() != NULL) {
			T2_b->cut_parent();
			num_cut++;
			cut_b = true;
		}
		// T2_b will move up after contraction
		else {
			T2_b = T2_b->parent();
		}
		// check for T2 parents as singletons
		Node *node = T2_ab->contract();
		if (node != NULL && node->is_singleton()
				&& node != T2->get_component(0))
			singletons->push_back(node);

		// if T2_c is gone then its replacement is in singleton list
		// contract might delete old T2_c, see where it is
		bool add_T2_c = true;
		T2_c = T1_c->get_twin();
		// ignore T2_c if it is a singleton
		if (T2_c != node && T2_c->parent() != NULL && !cut_b_only) {

			Node *T2_c_parent = T2_c->parent();
			T2_c->cut_parent();
			num_cut++;
			node = T2_c_parent->contract();
			if (node != NULL && node->is_singleton()
					&& node != T2->get_component(0))
				singletons->push_back(node);
		}
		else {
			add_T2_c = false;
		}

		
		if (!cut_b_only)
			T1->add_component(T1_a);
		if (!cut_b_only)
			T1->add_component(T1_c);
		// put T2 cut parts into T2
		if (!cut_b_only) {
			T2->add_component(T2_a);
		}
		// may have already been added
		if (cut_b) {
			T2->add_component(T2_b);
		}
		// problem if c is deleted
		if (add_T2_c) {
			T2->add_component(T2_c);
		}

		// may have already been added
		if (T2_b->is_leaf())
			singletons->push_back(T2_b);

	}
}
	}
// if the first component of the forests differ then we have to cut p
if (T1->get_component(0)->get_twin() != T2->get_component(0)) {
	num_cut++;
	T1->add_rho();
	T2->add_rho();
}
return num_cut;
}
/*******************************************************************************
	RSPR WORSE_3_APPROX
*******************************************************************************/

/* rSPR_worse_3_approx
 * Calculate an approximate maximum agreement forest and SPR distance
 * RETURN At most 3 times the rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 * T1 must be a binary tree. T2 can be a multifurcating forest.
 */
int rSPR_worse_3_approx(Forest *T1, Forest *T2) {
	return rSPR_worse_3_approx(T1, T2, true);
}

int rSPR_worse_3_approx(Forest *T1, Forest *T2, bool sync) {
	// match up nodes of T1 and T2
	if (sync) {
if (!sync_twins(T1, T2))
	return 0;
	}
//	cout << "T1: "; T1->print_components();
//	cout << "T2: "; T2->print_components();
	// find sibling pairs of T1
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	Forest *F1;
	Forest *F2;

	int ans = rSPR_worse_3_approx_hlpr(T1, T2, &singletons, sibling_pairs, &F1, &F2, true);

	F1->swap(T1);
	F2->swap(T2);
	sync_twins(T1,T2);


	delete sibling_pairs;
	delete F1;
	delete F2;
	return ans;
}

int rSPR_worse_3_approx_distance_only(Forest *T1, Forest *T2) {
if (!sync_twins(T1, T2))
	return 0;
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	list<Node *> singletons = T2->find_singletons();
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	int ans = rSPR_worse_3_approx_hlpr(T1, T2, &singletons, sibling_pairs, NULL, NULL, false);

	delete sibling_pairs;
	return ans;
}

int rSPR_worse_3_approx(Node *subtree, Forest *T1, Forest *T2) {
	return rSPR_worse_3_approx(subtree, T1, T2, true);
}

int rSPR_worse_3_approx(Node *subtree, Forest *T1, Forest *T2, bool sync) {
	// match up nodes of T1 and T2
	if (sync) {
if (!sync_twins(T1, T2))
	return 0;
	}
//	cout << "T1: "; T1->print_components();
//	cout << "T2: "; T2->print_components();
	// find sibling pairs of T1
	list<Node *> *sibling_pairs = subtree->find_sibling_pairs();
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	Forest *F1;
	Forest *F2;

	int ans = rSPR_worse_3_approx_hlpr(T1, T2, &singletons, sibling_pairs, &F1, &F2, true);

	F1->swap(T1);
	F2->swap(T2);
	sync_twins(T1,T2);


	delete sibling_pairs;
	delete F1;
	delete F2;
	return ans;
}

// rSPR_worse_3_approx recursive helper function
int rSPR_worse_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests) {
	#ifdef DEBUG_APPROX
cout << "rSPR_worse_3_approx_hlpr" << endl;
			cout << "\tT1: ";
			T1->print_components_with_twins();
			cout << "\tT2: ";
			T2->print_components_with_twins();
			cout << "sibling pairs:";
			for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
				cout << "  ";
				(*i)->print_subtree_hlpr();
			}
			cout << endl;
	#endif
	int num_cut = 0;
	UndoMachine um = UndoMachine();
	while(!singletons->empty() || !sibling_pairs->empty()) {
// Case 1 - Remove singletons
while(!singletons->empty()) {
	#ifdef DEBUG_APPROX
		cout << "Case 1" << endl;
	#endif

	Node *T2_a = singletons->back();
	singletons->pop_back();
	// find twin in T1
	Node *T1_a = T2_a->get_twin();
	// if this is in the first component of T_2 then
	// it is not really a singleton.
	// TODO: problem when we cluster and have a singleton as the
	//		first comp of T2
	//    NEED TO MODIFY CUTTING?
	// 		HERE AND IN BB?
	if (T2_a == T2->get_component(0))
		continue;

	Node *T1_a_parent = T1_a->parent();
	if (T1_a_parent == NULL)
		continue;
	bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
	// cut the edge above T1_a
	um.add_event(new CutParent(T1_a));
	T1_a->cut_parent();
	um.add_event(new AddComponent(T1));
	T1->add_component(T1_a);
	//if (T1_a->get_sibling_pair_status() > 0)
	//	T1_a->clear_sibling_pair(sibling_pairs);
	//delete(T1_a);

	ContractEvent(&um, T1_a_parent);
	Node *node = T1_a_parent->contract();
	if (node != NULL && potential_new_sibling_pair &&
			node->is_sibling_pair()){
		um.add_event(new AddToFrontSiblingPairs(sibling_pairs));
		sibling_pairs->push_front(node->rchild());
		sibling_pairs->push_front(node->lchild());
	}

	#ifdef DEBUG_APPROX
			cout << "\tT1: ";
			T1->print_components();
			cout << "\tT2: ";
			T2->print_components();
	#endif
}
if(!sibling_pairs->empty()) {
	/*
	if (PREORDER_SIBLING_PAIRS) {
		T1->get_component(0)->preorder_number();
		list<Node *>::iterator c;
		list<Node *>::iterator best_sib = sibling_pairs->end();
		int best_prenum = INT_MAX;
		list<Node *>::iterator T1_a_i;
		list<Node *>::iterator T1_c_i;
		for(c = sibling_pairs->begin(); c != sibling_pairs->end(); ) {
				T1_c_i = c;
				T1_c = *c;
				c++;
				T1_a_i = c;
				T1_a = *c;
				c++;
				cout << T1_a->str_subtree() << endl;
				cout << T1_c->str_subtree() << endl;
//					if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
//						cout << "invalid" << endl;
//						sibling_pairs->erase(T1_c_i);
//						sibling_pairs->erase(T1_a_i);
//						um.add_event(new PopSiblingPair(T1_a, T1_c, sibling_pairs));
//						continue;
//					}
//					else {
				//int prenum = T1_a->parent()->get_preorder_number();
				int prenum = T1_a->get_preorder_number();
				cout << "prenum=" << prenum << endl;
				cout << "old_prenum=" << best_prenum << endl;
				if (prenum < best_prenum) {
					best_sib = T1_c_i; 
					best_prenum = prenum;
				}
				cout << "new_prenum=" << best_prenum << endl;
//					}
		}
		cout << endl;
		if (best_prenum == INT_MAX)
			continue;
		else {
			T1_c_i = best_sib;
			T1_c = *T1_c_i;
			best_sib++;
			T1_a_i = best_sib;
			T1_a = *T1_a_i;
			sibling_pairs->erase(T1_a_i);
			sibling_pairs->erase(T1_c_i);
		}
	}
	else {
	*/
	Node *T1_a = sibling_pairs->back();
	sibling_pairs->pop_back();
	Node *T1_c = sibling_pairs->back();
	sibling_pairs->pop_back();
	um.add_event(new PopSiblingPair(T1_a, T1_c, sibling_pairs));

	//if (T1_a->get_sibling_pair_status() == 0 ||
	//		T1_c->get_sibling_pair_status() == 0) {
	//	continue;
//			}

	//T1_a->clear_sibling_pair_status();
	//T1_c->clear_sibling_pair_status();
	if (T1_a->parent() == NULL || T1_c->parent() == NULL || T1_a->parent() != T1_c->parent()) {
		continue;
	}
	if (!T1_a->can_be_sibling() || !T1_c->can_be_sibling()
	|| num_cut >= INT_MAX - 3) {
		continue;
	}
	Node *T1_ac = T1_a->parent();
	// lookup in T2 and determine the case
	Node *T2_a = T1_a->get_twin();
	Node *T2_c = T1_c->get_twin();

	#ifdef DEBUG_APPROX
		cout << "Fetching sibling pair" << endl;
		T1_ac->print_subtree();
		cout << "T2_a" << ": ";
		cout << " d=" << T2_a->get_depth() << " ";
		T2_a->print_subtree();
		cout << "T1_c" << ": ";
		T1_c->print_subtree();
		cout << "T2_c" << ": ";
		cout << " d=" << T2_c->get_depth() << " ";
		T2_c->print_subtree();
	#endif

	// Case 2 - Contract identical sibling pair
	if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
		#ifdef DEBUG_APPROX
			cout << "Case 2" << endl;
			T1->print_components();
			T2->print_components();
		#endif
		Node *T2_ac = T2_a->parent();
		um.add_event(new ContractSiblingPair(T1_ac));
		T1_ac->contract_sibling_pair_undoable();
		um.add_event(new ContractSiblingPair(T2_ac, T2_a, T2_c, &um));
		Node *T2_ac_new = T2_ac->contract_sibling_pair_undoable(T2_a, T2_c);
		if (T2_ac_new != NULL && T2_ac_new != T2_ac) {
			T2_ac = T2_ac_new;
			um.add_event(new CreateNode(T2_ac));
			um.add_event(new ContractSiblingPair(T2_ac));
			T2_ac->contract_sibling_pair_undoable();
		}
		um.add_event(new SetTwin(T1_ac));
		um.add_event(new SetTwin(T2_ac));
		T1_ac->set_twin(T2_ac);
		T2_ac->set_twin(T1_ac);
		//T2_ac->fix_contracted_order();
		//T1->add_deleted_node(T1_a);
		//T1->add_deleted_node(T1_c);
		//T2->add_deleted_node(T2_a);
		//T2->add_deleted_node(T2_c);

		// check if T2_ac is a singleton
		//if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
		if (T2_ac->is_singleton() && T1_ac != T1->get_component(0) && T2_ac != T2->get_component(0))
			singletons->push_back(T2_ac);
		// check if T1_ac is part of a sibling pair
		if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
			um.add_event(new AddToSiblingPairs(sibling_pairs));
			sibling_pairs->push_back(T1_ac->parent()->lchild());
			sibling_pairs->push_back(T1_ac->parent()->rchild());
		}
	}
	// Case 3
	else {
		#ifdef DEBUG_APPROX
			cout << "Case 3" << endl;
		#endif
		
		//  ensure T2_a is below T2_c
		if ((T2_a->get_depth() < T2_c->get_depth()
				&& T2_c->parent() != NULL)
				|| T2_a->parent() == NULL) {
			#ifdef DEBUG_APPROX
				cout << "swapping" << endl;
			#endif
		
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);

		}
		else if (T2_a->get_depth() == T2_c->get_depth()) {
			if (T2_a->parent() && T2_c->parent() &&
					(T2_a->parent()->get_depth() <
					T2_c->parent()->get_depth())) {
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);
			}
		}

		// get T2_b
		bool multi_node = false;
		Node *T2_ab = T2_a->parent();
		Node *T2_b = T2_ab;
		if (T2_ab->get_children().size() > 2) {
			multi_node = true;
		}
		else {
			T2_b = T2_ab->rchild();
			if (T2_b == T2_a)
				T2_b = T2_ab->lchild();
		}

		#ifdef DEBUG_APPROX
		cout << "T2_b" << ": ";
		cout.flush();
		T2_b->print_subtree();
	#endif
		// cut T1_a, T1_c, T2_a, T2_b, T2_c

		bool cut_a_only = false;
		bool cut_b_only = false;
		bool cut_c_only = false;
		bool cut_b_only_if_not_a_or_c = false;
		if (APPROX_CUT_ONE_B && T2_a->parent() != NULL && T2_a->parent()->parent() != NULL && T2_a->parent()->parent() == T2_c->parent() && !multi_node
						&& (!APPROX_EDGE_PROTECTION || !T2_b->is_protected())) {
			cut_b_only = true;
			um.add_event(new AddToSiblingPairs(sibling_pairs));
			sibling_pairs->push_back(T1_c);
			sibling_pairs->push_back(T1_a);
		}
	if (APPROX_CUT_TWO_B && !cut_b_only && T1_ac->parent() != NULL
						&& (!APPROX_EDGE_PROTECTION || !T2_b->is_protected())) {
		Node *T1_s = T1_ac->get_sibling();
		if (T1_s->is_leaf()) {
			Node *T2_l = T2_a->parent()->parent();
			if (T2_l != NULL) {
				if (T2_c->parent() != NULL && T2_c->parent()->parent() == T2_l
						&& T2_a->parent()->get_children().size() > 2
						&& T2_c->parent()->get_children().size() > 2) {
					if (T2_l->get_sibling() == T1_s->get_twin()) {
						cut_b_only=true;
					}
					else if (T2_l->parent() == NULL &&
							(T2->contains_rho() ||
							 T2->get_component(0) != T2_l)) {
						cut_b_only_if_not_a_or_c=true;
					}
				}
				else if ((T2_l = T2_l->parent()) != NULL
						&& T2_c->parent() == T2_l
						&& T2_a->parent()->get_children().size() > 2
						&& T2_a->parent()->parent()->get_children().size() > 2) {
					if (T2_l->get_sibling() == T1_s->get_twin()) {
						cut_b_only=true;
					}
					else if (T2_l->parent() == NULL &&
							(T2->contains_rho() ||
							 T2->get_component(0) != T2_l)) {
						cut_b_only_if_not_a_or_c=true;
					}
				}
			}
		}
	}
	if (APPROX_REVERSE_CUT_ONE_B && !cut_b_only && T1_ac->parent() != NULL) {
		Node *T1_s = T1_ac->get_sibling();
		if (T1_s->is_leaf()) {
			if (T1_s->get_twin()->parent() == T2_a->parent()//) {
						&& (!APPROX_EDGE_PROTECTION || !T2_c->is_protected())) {
				cut_c_only=true;
			}
			else if (T1_s->get_twin()->parent() == T2_c->parent()//) {
						&& (!APPROX_EDGE_PROTECTION || !T2_a->is_protected())
							&& T2_c->parent()->get_children().size() <= 2) {
				cut_a_only=true;
			}
		}
		else if (APPROX_REVERSE_CUT_ONE_B_2) {
			if (T2_c->parent() != NULL
				&& chain_match(T1_s, T2_c->get_sibling(), T2_a) //)
						&& (!APPROX_EDGE_PROTECTION || !T2_a->is_protected()))
			cut_a_only = true;
		}
	}
	if (APPROX_CUT_TWO_B_ROOT && cut_a_only == false && cut_c_only == false
			&& cut_b_only_if_not_a_or_c == true) {
		cut_b_only = true;
	}
	/*
	if (CUT_LOST) {
		if (T1_a->num_lost_children() > 0
				|| T2_a->num_lost_children() > 0) {
			cut_a_only = true;
			cut_b_only = false;
			cut_c_only = false;
			num_cut-=3;
		}
		else if (T1_c->num_lost_children() > 0
				|| T2_c->num_lost_children() > 0) {
			cut_a_only = false;
			cut_b_only = false;
			cut_c_only = true;
			num_cut-=3;
		}
		else if (T2_b->is_leaf()) {
			if (T2_b->num_lost_children() > 0
				|| T2_b->get_twin()->num_lost_children() > 0) {
			cut_a_only = false;
			cut_b_only = true;
			cut_c_only = false;
			num_cut-=3;
			}
		}
	}
	*/


		Node *node;

		bool cut_a = false;
		bool cut_c = false;
		if (!cut_b_only || T2_a->parent()->get_children().size() > 2) {
			if (!cut_c_only &&
					(!APPROX_EDGE_PROTECTION
					 	|| (!T2_a->is_protected()
							&& (T2_a->parent()->parent() != NULL
								|| !T2_b->is_protected()
								|| T2_a->parent()->get_children().size() > 2)))) {
//					|| cut_a_only)) {
				um.add_event(new CutParent(T1_a));
				T1_a->cut_parent();
				cut_a = true;

				ContractEvent(&um, T1_ac);
				node = T1_ac->contract();
			}
			else
				node = T1_ac;
			if (!cut_a_only &&
					(!APPROX_EDGE_PROTECTION
					 	|| (!T2_c->is_protected()
						&& (T2_c->parent() == NULL
								|| T2_c->parent()->parent() != NULL
								|| !T2_c->get_sibling()->is_protected()
								|| T2_c->parent()->get_children().size() > 2)))) {// &&
//					|| cut_c_only)) {
				um.add_event(new CutParent(T1_c));
				T1_c->cut_parent();
				cut_c = true;

				if (node) {
					ContractEvent(&um, node);
					node = node->contract();
				}
			}

			// contract parents
			// check for T1_ac sibling pair
			if (node && node->is_sibling_pair()){
				um.add_event(new AddToSiblingPairs(sibling_pairs));
				sibling_pairs->push_back(node->lchild());
				sibling_pairs->push_back(node->rchild());
			}
		}

		bool same_component = true;
		if (APPROX_CHECK_COMPONENT && !cut_a_only && !cut_c_only)
			same_component = (T2_a->find_root() == T2_c->find_root());

		Node *T2_ab_parent = T2_ab->parent();
		node = T2_ab;
		if (cut_a) {
			um.add_event(new CutParent(T2_a));
			T2_a->cut_parent();

			//ContractEvent(&um, T2_ab);
			//node = T2_ab->contract();
		}
		bool cut_b = false;
		if (same_component && T2_ab_parent != NULL
				&& !cut_a_only && !cut_c_only
				&& (!APPROX_EDGE_PROTECTION
					|| (!T2_b->is_protected() ))) {
//							&& (T2_b->parentT2_a->parent()->parent() != NULL
//								|| !T2_a->is_protected())))) {
//					|| cut_b_only)) {
			if (multi_node) {
				T2_b = T2_ab;
				um.add_event(new CutParent(T2_ab));
				T2_ab->cut_parent();
				if (T2_a->parent() != NULL) {
					um.add_event(new CutParent(T2_a));
					T2_a->cut_parent();
					um.add_event(new AddChild(T2_a));
					T2_ab_parent->add_child(T2_a);
				}
				else
					node = T2_ab_parent;
			}
			else {
				um.add_event(new CutParent(T2_b));
				T2_b->cut_parent();
				//ContractEvent(&um, node);
				//node = node->contract();
			}
			cut_b = true;
		}
		// T2_b will move up after contraction
		else if (!multi_node) {
			T2_b = T2_b->parent();
		}
		if (node != NULL) {
			ContractEvent(&um, node);
			node = node->contract();
		// check for T2 parents as singletons
		if (node != NULL && node->is_singleton()
				&& node != T2->get_component(0))
			singletons->push_back(node);
		}

		// if T2_c is gone then its replacement is in singleton list
		// contract might delete old T2_c, see where it is
		bool add_T2_c = true;
		T2_c = T1_c->get_twin();
		// ignore T2_c if it is a singleton
		if (cut_c && T2_c != node && T2_c->parent() != NULL) {
			Node *T2_c_parent = T2_c->parent();
			um.add_event(new CutParent(T2_c));
			T2_c->cut_parent();
			ContractEvent(&um, T2_c_parent);
			node = T2_c_parent->contract();
			if (node != NULL && node->is_singleton()
					&& node != T2->get_component(0))
				singletons->push_back(node);
		}
		else {
			add_T2_c = false;
		}

		
		if (cut_a) {
			um.add_event(new AddComponent(T1));
			T1->add_component(T1_a);
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_a);
		}
		if (cut_c) {
			um.add_event(new AddComponent(T1));
			T1->add_component(T1_c);
		}
		if (cut_b) {
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_b);
		}
		// problem if c is deleted
		if (add_T2_c) {
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_c);
		}

		// may have already been added
		if (T2_b->is_leaf() && cut_b)
			singletons->push_back(T2_b);

		num_cut+=3;

		if (cut_a == false && cut_b == false && cut_c == false) {
			num_cut = INT_MAX-3;
		}

	}
}
	}
// if the first component of the forests differ then we have cut p
if (T1->get_component(0)->get_twin() != T2->get_component(0)) {
	if (!T1->contains_rho()) {
		um.add_event(new AddRho(T1));
		um.add_event(new AddRho(T2));
		T1->add_rho();
		T2->add_rho();
	}
	else
		// hack to ignore rho when it shouldn't be in a cluster
		num_cut -=3;
}
if (save_forests) {
	*F1 = new Forest(T1);
	*F2 = new Forest(T2);
}
 
#ifdef DEBUG_APPROX
#ifdef DEBUG_UNDO
 while(um.num_events() > 0) {
		cout << "Undo step " << um.num_events() << endl;
		cout << "T1: ";
		T1->print_components();
		cout << "T2: ";
		T2->print_components();
			cout << "sibling pairs:";
			for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
				cout << "  ";
				(*i)->print_subtree_hlpr();
			}
			cout << endl;
	 um.undo();
	cout << endl;
 }
#else
 um.undo_all();
#endif
#else
 um.undo_all();
#endif

 
//		 for(int i = 0; i < T1->num_components(); i++)
//		 	T1->get_component(i)->fix_parents();
//		 for(int i = 0; i < T2->num_components(); i++)
//		 	T2->get_component(i)->fix_parents();
return num_cut;
}

/*******************************************************************************
	RSPR WORSE_3_APPROX_BINARY
*******************************************************************************/

int rSPR_worse_3_approx_binary(Forest *T1, Forest *T2, bool sync) {
	// match up nodes of T1 and T2
	if (sync) {
if (!sync_twins(T1, T2))
	return 0;
	}
	// find sibling pairs of T1
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	Forest *F1;
	Forest *F2;

	int ans = rSPR_worse_3_approx_binary_hlpr(T1, T2, &singletons, sibling_pairs, &F1, &F2, true);

	F1->swap(T1);
	F2->swap(T2);
	sync_twins(T1,T2);


	delete sibling_pairs;
	delete F1;
	delete F2;
	return ans;
}

// rSPR_worse_3_approx_binary recursive helper function
int rSPR_worse_3_approx_binary_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests) {
	#ifdef DEBUG_APPROX
cout << "rSPR_worse_3_approx_binary_hlpr" << endl;
			cout << "\tT1: ";
			T1->print_components_with_twins();
			cout << "\tT2: ";
			T2->print_components_with_twins();
			cout << "sibling pairs:";
			for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
				cout << "  ";
				(*i)->print_subtree_hlpr();
			}
			cout << endl;
	#endif
	int num_cut = 0;
	UndoMachine um = UndoMachine();
	while(!singletons->empty() || !sibling_pairs->empty()) {
// Case 1 - Remove singletons
while(!singletons->empty()) {
	#ifdef DEBUG_APPROX
		cout << "Case 1" << endl;
	#endif

	Node *T2_a = singletons->back();
	singletons->pop_back();
	// find twin in T1
	Node *T1_a = T2_a->get_twin();
	// if this is in the first component of T_2 then
	// it is not really a singleton.
	// TODO: problem when we cluster and have a singleton as the
	//		first comp of T2
	//    NEED TO MODIFY CUTTING?
	// 		HERE AND IN BB?
	if (T2_a == T2->get_component(0))
		continue;

	Node *T1_a_parent = T1_a->parent();
	if (T1_a_parent == NULL)
		continue;
	bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
	// cut the edge above T1_a
	um.add_event(new CutParent(T1_a));
	T1_a->cut_parent();
	um.add_event(new AddComponent(T1));
	T1->add_component(T1_a);
	//if (T1_a->get_sibling_pair_status() > 0)
	//	T1_a->clear_sibling_pair(sibling_pairs);
	//delete(T1_a);

	ContractEvent(&um, T1_a_parent);
	Node *node = T1_a_parent->contract();
	if (node != NULL && potential_new_sibling_pair && node->is_sibling_pair()){
		um.add_event(new AddToFrontSiblingPairs(sibling_pairs));
		sibling_pairs->push_front(node->rchild());
		sibling_pairs->push_front(node->lchild());
	}

	#ifdef DEBUG_APPROX
			cout << "\tT1: ";
			T1->print_components();
			cout << "\tT2: ";
			T2->print_components();
	#endif
}
if(!sibling_pairs->empty()) {
	Node *T1_a = sibling_pairs->back();
	sibling_pairs->pop_back();
	Node *T1_c = sibling_pairs->back();
	sibling_pairs->pop_back();
	um.add_event(new PopSiblingPair(T1_a, T1_c, sibling_pairs));

	//if (T1_a->get_sibling_pair_status() == 0 ||
	//		T1_c->get_sibling_pair_status() == 0) {
	//	continue;
	//}

	//T1_a->clear_sibling_pair_status();
	//T1_c->clear_sibling_pair_status();
	if (T1_a->parent() == NULL || T1_c->parent() == NULL || T1_a->parent() != T1_c->parent()) {
		continue;
	}
	Node *T1_ac = T1_a->parent();
	// lookup in T2 and determine the case
	Node *T2_a = T1_a->get_twin();
	Node *T2_c = T1_c->get_twin();

	#ifdef DEBUG_APPROX
		cout << "Fetching sibling pair" << endl;
		T1_ac->print_subtree();
		cout << "T2_a" << ": ";
		cout << " d=" << T2_a->get_depth() << " ";
		T2_a->print_subtree();
		cout << "T1_c" << ": ";
		T1_c->print_subtree();
		cout << "T2_c" << ": ";
		cout << " d=" << T2_c->get_depth() << " ";
		T2_c->print_subtree();
	#endif

	// Case 2 - Contract identical sibling pair
	if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
		#ifdef DEBUG_APPROX
			cout << "Case 2" << endl;
			T1->print_components();
			T2->print_components();
		#endif
		Node *T2_ac = T2_a->parent();
		um.add_event(new ContractSiblingPair(T1_ac));
		um.add_event(new ContractSiblingPair(T2_ac));
		T1_ac->contract_sibling_pair_undoable();
		T2_ac->contract_sibling_pair_undoable();
		um.add_event(new SetTwin(T1_ac));
		um.add_event(new SetTwin(T2_ac));
		T1_ac->set_twin(T2_ac);
		T2_ac->set_twin(T1_ac);
		//T1->add_deleted_node(T1_a);
		//T1->add_deleted_node(T1_c);
		//T2->add_deleted_node(T2_a);
		//T2->add_deleted_node(T2_c);

		// check if T2_ac is a singleton
		if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
			singletons->push_back(T2_ac);
		// check if T1_ac is part of a sibling pair
		if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
			um.add_event(new AddToSiblingPairs(sibling_pairs));
			sibling_pairs->push_back(T1_ac->parent()->lchild());
			sibling_pairs->push_back(T1_ac->parent()->rchild());
		}
	}
	// Case 3
	else {
		#ifdef DEBUG_APPROX
			cout << "Case 3" << endl;
		#endif
		
		//  ensure T2_a is below T2_c
		if ((T2_a->get_depth() < T2_c->get_depth()
				&& T2_c->parent() != NULL)
				|| T2_a->parent() == NULL) {
			#ifdef DEBUG_APPROX
				cout << "swapping" << endl;
			#endif
		
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);

		}
		else if (T2_a->get_depth() == T2_c->get_depth()) {
			if (T2_a->parent() && T2_c->parent() &&
					(T2_a->parent()->get_depth() <
					T2_c->parent()->get_depth())) {
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);
			}
		}

		// get T2_b
		Node *T2_ab = T2_a->parent();
		Node *T2_b = T2_ab->rchild();
		if (T2_b == T2_a)
			T2_b = T2_ab->lchild();

		#ifdef DEBUG_APPROX
		cout << "T2_b" << ": ";
		cout.flush();
		T2_b->print_subtree();
	#endif
		// cut T1_a, T1_c, T2_a, T2_b, T2_c

		bool cut_b_only = false;
		if (T2_a->parent() != NULL && T2_a->parent()->parent() != NULL && T2_a->parent()->parent() == T2_c->parent()) {
			cut_b_only = true;
			um.add_event(new AddToSiblingPairs(sibling_pairs));
			sibling_pairs->push_back(T1_c);
			sibling_pairs->push_back(T1_a);
		}

		Node *node;

		if (!cut_b_only) {
			um.add_event(new CutParent(T1_a));
			T1_a->cut_parent();

			ContractEvent(&um, T1_ac);
			node = T1_ac->contract();

			um.add_event(new CutParent(T1_c));
			T1_c->cut_parent();


			ContractEvent(&um, node);
			node = node->contract();

			// contract parents
			// check for T1_ac sibling pair
			if (node && node->is_sibling_pair()){
				um.add_event(new AddToSiblingPairs(sibling_pairs));
				sibling_pairs->push_back(node->lchild());
				sibling_pairs->push_back(node->rchild());
			}
		}

		bool same_component = true;
		if (APPROX_CHECK_COMPONENT)
			same_component = (T2_a->find_root() == T2_c->find_root());

		Node *T2_ab_parent = T2_ab->parent();
		node = T2_ab;
		if (!cut_b_only) {
			um.add_event(new CutParent(T2_a));
			T2_a->cut_parent();

			//ContractEvent(&um, T2_ab);
			//node = T2_ab->contract();
		}
		bool cut_b = false;
		if (same_component && T2_ab_parent != NULL) {
			um.add_event(new CutParent(T2_b));
			T2_b->cut_parent();
			//ContractEvent(&um, node);
			//node = node->contract();
			cut_b = true;
		}
		// T2_b will move up after contraction
		else {
			T2_b = T2_b->parent();
		}
			ContractEvent(&um, node);
			node = node->contract();
		// check for T2 parents as singletons
		if (node != NULL && node->is_singleton()
				&& node != T2->get_component(0))
			singletons->push_back(node);

		// if T2_c is gone then its replacement is in singleton list
		// contract might delete old T2_c, see where it is
		bool add_T2_c = true;
		T2_c = T1_c->get_twin();
		// ignore T2_c if it is a singleton
		if (T2_c != node && T2_c->parent() != NULL && !cut_b_only) {

			Node *T2_c_parent = T2_c->parent();
			um.add_event(new CutParent(T2_c));
			T2_c->cut_parent();
			ContractEvent(&um, T2_c_parent);
			node = T2_c_parent->contract();
			if (node != NULL && node->is_singleton()
					&& node != T2->get_component(0))
				singletons->push_back(node);
		}
		else {
			add_T2_c = false;
		}

		
		if (!cut_b_only) {
			um.add_event(new AddComponent(T1));
			T1->add_component(T1_a);
			um.add_event(new AddComponent(T1));
			T1->add_component(T1_c);
			// put T2 cut parts into T2
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_a);
			// may have already been added
		}
		if (cut_b) {
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_b);
		}
		// problem if c is deleted
		if (add_T2_c) {
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_c);
		}

		// may have already been added
		if (T2_b->is_leaf())
			singletons->push_back(T2_b);

		num_cut+=3;

	}
}
	}
// if the first component of the forests differ then we have cut p
if (T1->get_component(0)->get_twin() != T2->get_component(0)) {
	if (!T1->contains_rho()) {
		um.add_event(new AddRho(T1));
		um.add_event(new AddRho(T2));
		T1->add_rho();
		T2->add_rho();
	}
	else
		// hack to ignore rho when it shouldn't be in a cluster
		num_cut -=3;
}
if (save_forests) {
	*F1 = new Forest(T1);
	*F2 = new Forest(T2);
}
 

#ifdef DEBUG_UNDO
 while(um.num_events() > 0) {
		cout << "Undo step " << um.num_events() << endl;
		cout << "T1: ";
		T1->print_components();
		cout << "T2: ";
		T2->print_components();
			cout << "sibling pairs:";
			for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
				cout << "  ";
				(*i)->print_subtree_hlpr();
			}
			cout << endl;
	 um.undo();
	cout << endl;
 }
#else
 um.undo_all();
#endif

 
//		 for(int i = 0; i < T1->num_components(); i++)
//		 	T1->get_component(i)->fix_parents();
//		 for(int i = 0; i < T2->num_components(); i++)
//		 	T2->get_component(i)->fix_parents();
return num_cut;
}


int rSPR_branch_and_bound(Forest *T1, Forest *T2) {
	return rSPR_branch_and_bound_range(T1, T2, MAX_SPR);
}


int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int end_k) {
	string problem_key;
	map<string,ProblemSolution>::iterator i;

	if (MEMOIZE) {
problem_key = T1->str() + ":" + T2->str();
i = memoized_clusters.find(problem_key);
if (i != memoized_clusters.end()) {
	//cout << "already solved: " << endl;
	//cout << problem_key << endl;
	//cout << i->second.T2 << endl;
	//cout << "start" << endl;
	Forest *new_T1 = build_finished_forest(i->second.T1);
	//cout << "middle" << endl;
	Forest *new_T2 = build_finished_forest(i->second.T2);
	//cout << "end" << endl;
	T1->swap(new_T1);
	T2->swap(new_T2);
	sync_twins(T1, T2);
	delete new_T1;
	delete new_T2;
	return i->second.k;
}
	}
	Forest F1 = Forest(T1);
	Forest F2 = Forest(T2);
	int approx_spr = rSPR_worse_3_approx(&F1, &F2);
	int min_spr = approx_spr / 3;
	int exact_spr = rSPR_branch_and_bound_range(T1, T2, min_spr, end_k);
	if (MEMOIZE && exact_spr >= 0 && i == memoized_clusters.end()) {
//string solution_key = T1->str() + ":" + T2->str();
memoized_clusters.insert(make_pair(problem_key,
		ProblemSolution(T1,T2,exact_spr)));
	}

	return exact_spr;
}
	
int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int start_k,
int end_k) {
	int exact_spr = -1;
	bool in_main = MAIN_CALL;
	MAIN_CALL = false;
	int k;
	for(k = start_k; k <= end_k; k++) {
if (in_main) {
	cout << " " << k;
	cout.flush();
}
//Forest F1 = Forest(T1);
//Forest F2 = Forest(T2);
//exact_spr = rSPR_branch_and_bound(&F1, &F2, k);
exact_spr = rSPR_branch_and_bound(T1,T2, k);
//if (exact_spr >= 0 || k == end_k) {
if (exact_spr >= 0) {
//			F1.swap(T1);
//			F2.swap(T2);
	break;
}
	}
	if (in_main)
cout << endl;
	if (k > end_k)
k = -1;
	return k;
}

int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k) {
	return rSPR_branch_and_bound(T1, T2, k, NULL, NULL);
}

/* rSPR_branch_and_bound
 * Calculate a maximum agreement forest and SPR distance
 * Uses a branch and bound optimization to not explore paths
 * guaranteed to be incorrect based on rspr_3_approx
 * RETURN The rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 */
int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k,
		map<string, int> *label_map,
		map<int, string> *reverse_label_map) {
	// find sibling pairs of T1
//	cout << "foo1" << endl;
	if (!sync_twins(T1, T2))
return 0;
	if (PREORDER_SIBLING_PAIRS &&
			T1->get_component(0)->get_preorder_number() == -1) {
		T1->get_component(0)->preorder_number();
		T2->get_component(0)->preorder_number();
}
	if (DEEPEST_PROTECTED_ORDER
			&& T1->get_component(0)->get_edge_pre_start() == -1) {
		T1->get_component(0)->edge_preorder_interval();
		T2->get_component(0)->edge_preorder_interval();
	}

	set<SiblingPair> *sibling_pairs;
	list<Node *> singletons;
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();
	sibling_pairs = find_sibling_pairs_set(T1);
	singletons = T2->find_singletons();
	list<Node *> protected_stack = list<Node *>();
	int num_ties = 2;


	int final_k = 
rSPR_branch_and_bound_hlpr(T1, T2, k, sibling_pairs, &singletons, false, &AFs, &protected_stack, &num_ties);

//		cout << "foo" << endl;
	// TODO: this is a cheap hack
	if (!AFs.empty()) {
if (ALL_MAFS
#ifdef DEBUG
		|| true
#endif
		) {
	cout << endl << endl << "FOUND ANSWERS" << endl;
	// TODO: this is a cheap hack
	for (list<pair<Forest,Forest> >::iterator x = AFs.begin(); x != AFs.end(); x++) {
		if (label_map != NULL && reverse_label_map != NULL) {
			x->first.numbers_to_labels(reverse_label_map);
			x->second.numbers_to_labels(reverse_label_map);
		}
		cout << "\tT1: ";
		x->first.print_components();
		cout << "\tT2: ";
		x->second.print_components();
		if (label_map != NULL && reverse_label_map != NULL) {
			x->first.labels_to_numbers(label_map, reverse_label_map);
			x->second.labels_to_numbers(label_map, reverse_label_map);
		}
	}
}
AFs.front().first.swap(T1);
AFs.front().second.swap(T2);
sync_twins(T1,T2);
	}
	if (final_k >= 0)
final_k = k - final_k;
	delete sibling_pairs;
	return final_k;
}

void add_sibling_pair(set<SiblingPair> *sibling_pairs, Node *a, Node *c, UndoMachine *um) {
	SiblingPair sp = SiblingPair(a,c);
	pair< set<SiblingPair>::iterator, bool> ins = 
	sibling_pairs->insert(sp);
	if (ins.second == false) {
um->add_event(new RemoveSetSiblingPairs(sibling_pairs, *(ins.first)));
sibling_pairs->erase(ins.first);
ins = sibling_pairs->insert(sp);
	}
	um->add_event(new AddToSetSiblingPairs(sibling_pairs, *(ins.first)));
}

SiblingPair pop_sibling_pair(set<SiblingPair> *sibling_pairs, UndoMachine *um) {
	set<SiblingPair>::iterator s = sibling_pairs->begin();
	SiblingPair spair = SiblingPair(*s); 
	um->add_event(new RemoveSetSiblingPairs(sibling_pairs, spair));
	sibling_pairs->erase(s);
	return spair;
}

SiblingPair pop_sibling_pair(set<SiblingPair>::iterator s, set<SiblingPair> *sibling_pairs, UndoMachine *um) {
	SiblingPair spair = SiblingPair(*s); 
	um->add_event(new RemoveSetSiblingPairs(sibling_pairs, spair));
	sibling_pairs->erase(s);
	return spair;
}

inline int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
set<SiblingPair> *sibling_pairs, list<Node *> *singletons,
bool cut_b_only, list<pair<Forest,Forest> > *AFs,
list<Node *> *protected_stack, int *num_ties) {
	return rSPR_branch_and_bound_hlpr(T1, T2, k, sibling_pairs,
			singletons, cut_b_only, AFs, protected_stack, num_ties, NULL, NULL);
}

// rSPR_branch_and_bound recursive helper function
int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
set<SiblingPair> *sibling_pairs, list<Node *> *singletons,
bool cut_b_only, list<pair<Forest,Forest> > *AFs,
list<Node *> *protected_stack, int *num_ties, Node *prev_T1_a, Node *prev_T1_c) {
	#ifdef DEBUG
	cout << "rSPR_branch_and_bound_hlpr()" << endl;
	cout << "\tT1: ";
	T1->print_components();
	cout << "\tT2: ";
	T2->print_components();
	cout << "K=" << k << endl;
	cout << "sibling pairs:";
	for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
cout << "  ";
(*i).a->print_subtree_hlpr();
cout << ",";
(*i).c->print_subtree_hlpr();
	}
	cout << endl;
	cout << "protected_stack:";
	for (list<Node *>::iterator i = protected_stack->begin(); i != protected_stack->end(); i++) {
cout << "  ";
(*i)->print_subtree_hlpr();
	}
	cout << endl;
	#endif

	UndoMachine um = UndoMachine();


	while(!singletons->empty() || !sibling_pairs->empty()) {
		// Case 1 - Remove singletons
		while(!singletons->empty()) {
			Node *T2_a = singletons->back();
			#ifdef DEBUG
				cout << "Case 1" << endl;
				cout << "a " << T2_a->str_subtree() << endl;
			#endif
		
			singletons->pop_back();
			// find twin in T1
			Node *T1_a = T2_a->get_twin();
		
			//if (T1_a->get_sibling_pair_status() > 0)
			//	cout << T1_a->get_sibling(sibling_pairs) << endl;
			// if this is in the first component of T_2 then
			// it is not really a singleton.
			Node *T1_a_parent = T1_a->parent();
			if (T1_a_parent == NULL)
				continue;
			bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
			if (T2_a == T2->get_component(0)) {
				// TODO: should we do this when it happens?
				if (!T1->contains_rho()) {
					um.add_event(new AddRho(T1));
					um.add_event(new AddRho(T2));
					T1->add_rho();
					T2->add_rho();
					k--;
					#ifdef DEBUG
					cout << "adding p element, k=" << k << endl;
					#endif
				}
			}
		
			// cut the edge above T1_a
			um.add_event(new CutParent(T1_a));
			T1_a->cut_parent();
		
			um.add_event(new AddComponent(T1));
			T1->add_component(T1_a);
			ContractEvent(&um, T1_a_parent);
		
		
			Node *node = T1_a_parent->contract();
		
			if (node != NULL && potential_new_sibling_pair && node->is_sibling_pair()){
				add_sibling_pair(sibling_pairs, node->lchild(), node->rchild(),
						&um);
			}
			#ifdef DEBUG
				cout << "\tT1: ";
				T1->print_components();
				cout << "\tT2: ";
				T2->print_components();
			#endif
		//			sp_i = sibling_pairs->begin();
		}
		if(!sibling_pairs->empty()) {
			Node *T1_a;
			Node *T1_c;
			set<SiblingPair>::iterator deepest_valid = sibling_pairs->end();
			int deepest_depth = INT_MAX;
			int deepest_depth_2 = INT_MAX;
			Node *best_a = NULL;
			Node *best_c = NULL;
			/* pop protected_stack when out of order sibling pairs
			 have already contracted it */
			if(!protected_stack->empty()
					&& protected_stack->back()->get_twin()->parent() == NULL
					&& protected_stack->back()->get_twin() != T1->get_component(0)) {
				um.undo_all();
				return -1;
			}
			while(!protected_stack->empty()
					&& (protected_stack->back()->is_contracted()
					// this shouldn't happen
						|| protected_stack->back()->get_twin()->parent() == NULL)) {
				um.add_event(new ListPopBack(protected_stack));
				protected_stack->pop_back();
			}
			if (LEAF_REDUCTION && !cut_b_only) {
				bool found = false;
				set<SiblingPair>::iterator sp_i = sibling_pairs->begin();
				// correct in case sibling pair involves previous
		/*				if (sp_i != sibling_pairs->begin()) {
					if (check_all_pairs)
						sp_i = sibling_pairs->begin();
					else
						sp_i--;
				}
				*/
				while (sp_i != sibling_pairs->end()) {
					T1_a = (*sp_i).a;
					T1_c = (*sp_i).c;
					if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
						um.add_event(new RemoveSetSiblingPairs(sibling_pairs,
									SiblingPair(T1_a, T1_c)));
						set<SiblingPair>::iterator rem = sp_i;
						sp_i++;
						sibling_pairs->erase(rem);
						continue;
					}
					Node *T2_a = T1_a->get_twin();
					Node *T2_c = T1_c->get_twin();
					// select if this is a Case 2 or, optionally, nonbranching
					if ((T2_a->parent() != NULL && T2_a->parent() == T2_c->parent())
							|| (!cut_b_only && PREFER_NONBRANCHING
									&& is_nonbranching(T1, T2, T1_a, T1_c, T2_a, T2_c))) {
						um.add_event(new RemoveSetSiblingPairs(sibling_pairs,
									SiblingPair(T1_a, T1_c)));
						set<SiblingPair>::iterator rem = sp_i;
						sp_i++;
						sibling_pairs->erase(rem);
						found = true;
						break;
					}
					if (DEEPEST_ORDER) {
						int depth;
						int depth2;
						if (T1_a->get_depth() < T1_c->get_depth())
							depth = T1_c->get_depth();
						else
							depth = T1_a->get_depth();
						if (T2_a->get_depth() < T2_c->get_depth())
							depth2 = T2_c->get_depth();
						else
							depth2 = T2_a->get_depth();
						if (deepest_valid == sibling_pairs->end()
								|| deepest_depth < depth
								|| (deepest_depth == depth && deepest_depth_2 < depth2) ) {
							// TODO: this crashes on bigtest2
							// Why can we end up cutting the protected node?
							if (!DEEPEST_PROTECTED_ORDER
									|| protected_stack->empty()
									|| (protected_stack->back()->get_twin()->parent()->get_edge_pre_start()
											<= T1_a->get_preorder_number()
										&& protected_stack->back()->get_twin()->parent()->get_edge_pre_end()
											>= T1_a->get_preorder_number())) {
								deepest_valid = sp_i;
								deepest_depth = depth;
								deepest_depth_2 = depth2;
							}
						}
					}
					/* TODO: create a stack of protected nodes (intervals?) and only
						 accept the deepest sibling pair within the interval
					*/

					/* TODO: remember to pop the stack when we include the protected
					   node */
					sp_i++;
				}
				if (!found) {
					if (sibling_pairs->empty())
						continue;
					else {
						SiblingPair spair;
//						cout << "depth: " << deepest_depth << endl;
						if (DEEPEST_ORDER && deepest_valid != sibling_pairs->end())
							spair = pop_sibling_pair(deepest_valid, sibling_pairs, &um);
						else
							spair = pop_sibling_pair(sibling_pairs, &um);
						T1_a = spair.a;
						T1_c = spair.c;
					}
				}
			}
			else {
				if (prev_T1_a != NULL && prev_T1_c != NULL) {
					T1_a = prev_T1_a;
					T1_c = prev_T1_c;
					prev_T1_a = NULL;
					prev_T1_c = NULL;
				}
				else {
					SiblingPair spair = pop_sibling_pair(sibling_pairs, &um);
					T1_a = spair.a;
					T1_c = spair.c;
				}
			}
//			if (T1_a->parent() != NULL)
//				cout << "a_p: " << T1_a->parent()->str_subtree() << endl;
//			if (T1_c->parent() != NULL)
//				cout << "c_p: " << T1_c->parent()->str_subtree() << endl;
			if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
				continue;
			}
			if (!T1_a->can_be_sibling() || !T1_c->can_be_sibling()) {
				continue;
			}
			Node *T1_ac = T1_a->parent();
			// lookup in T2 and determine the case
			Node *T2_a = T1_a->get_twin();
			Node *T2_c = T1_c->get_twin();

			if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
				#ifdef DEBUG
					cout << "Case 2" << endl;
					T1_ac->print_subtree();
				#endif
				Node *T2_ac = T2_a->parent();

				if (CHECK_MERGE_DEPTH &&
						(T2_a->get_max_merge_depth() > T2_ac->get_depth()
							|| T2_c->get_max_merge_depth() > T2_ac->get_depth())) {
					um.undo_all();
					return -1;
				}

				if (!protected_stack->empty() &&
						(T2_a == protected_stack->back()
						 	|| T2_c == protected_stack->back())) {
					um.add_event(new ListPopBack(protected_stack));
					protected_stack->pop_back();
				}
				// CAN THIS HAPPEN TWICE?
				if (!protected_stack->empty() &&
						(T2_a == protected_stack->back()
						 	|| T2_c == protected_stack->back())) {
					um.add_event(new ListPopBack(protected_stack));
					protected_stack->pop_back();
				}


				um.add_event(new ContractSiblingPair(T1_ac));
				T1_ac->contract_sibling_pair_undoable();
				um.add_event(new ContractSiblingPair(T2_ac, T2_a, T2_c, &um));
				Node *T2_ac_new = T2_ac->contract_sibling_pair_undoable(T2_a, T2_c);
				if (T2_ac_new != NULL && T2_ac_new != T2_ac) {
					T2_ac = T2_ac_new;
					um.add_event(new CreateNode(T2_ac));
					um.add_event(new ContractSiblingPair(T2_ac));
					T2_ac->contract_sibling_pair_undoable();
				}

				um.add_event(new SetTwin(T1_ac));
				um.add_event(new SetTwin(T2_ac));
				T1_ac->set_twin(T2_ac);
				T2_ac->set_twin(T1_ac);
				//T1->add_deleted_node(T1_a);
				//T1->add_deleted_node(T1_c);
				//T2->add_deleted_node(T2_a);
				//T2->add_deleted_node(T2_c);

				// check if T2_ac is a singleton
				if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
					singletons->push_back(T2_ac);
				// check if T1_ac is part of a sibling pair
				if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
				add_sibling_pair(sibling_pairs, T1_ac->parent()->lchild(), T1_ac->parent()->rchild(),
						&um);
				}
				#ifdef DEBUG
					cout << "\tT1: ";
					T1->print_components();
					cout << "\tT2: ";
					T2->print_components();
				#endif
			}
			/* need to copy trees and lists for branching
			 * use forest copy constructor for T1 and T2 giving T1' and T2'
			 * T1' twins are in T2, and same for T2' and T1.
			 * singleton list will be empty except for maybe above the cut,
			 * so this can be created.
			 * fix one set of twins (T2->T1' or T1->T2' not sure)
			 * exploit chained twin relationship to copy sibling pair list
			 * fix other set of twins
			 * swap T2 and T2' root nodes
			 * now do the cut
			 *
			 * note: don't copy for 3rd cut, is a waste
			 */

			// Case 3
			// note: guaranteed that singleton list is empty
			else {
				if (k <= 0) {
					if ((!CUT_LOST || k < 0 ||
								(T1_a->num_lost_children() == 0 &&
								 T1_c->num_lost_children() == 0))
							&& ((T2_c->parent() != NULL && T2_a->parent() != NULL)|| !T2->contains_rho())) {
						singletons->clear();
						um.undo_all();
						return k-1;
					}
				}
				Forest *best_T1;
				Forest *best_T2;
				int best_k = -1;
				int answer_a = -1;
				int answer_b = -1;
				int answer_c = -1;
				bool cut_ab_only = false;
				bool cut_a_only = false;
				bool cut_c_only = false;
				bool cut_a_or_merge_ac = false;
				bool same_component = true;
				int lca_depth = -1;
				int path_length = -1;
				bool cut_b_only_if_not_a_or_c = false;
				bool cob = false;
				int undo_state = um.num_events();
				//  ensure T2_a is below T2_c
				if ((T2_a->get_depth() < T2_c->get_depth()
						&& T2_c->parent() != NULL)
						|| T2_a->parent() == NULL) {
					swap(&T1_a, &T1_c);
					swap(&T2_a, &T2_c);
				}
				else if (T2_a->get_depth() == T2_c->get_depth()) {
					if (T2_a->parent() && T2_c->parent() &&
							(T2_a->parent()->get_depth() <
							T2_c->parent()->get_depth()
							//|| (T2_c->parent()->parent()
							//&& T2_c->parent()->parent() == T2_a->parent())
							)) {
					swap(&T1_a, &T1_c);
					swap(&T2_a, &T2_c);
					}
				}
				Node *T2_b = T2_a->parent()->rchild();
				if (T2_b == T2_a)
					T2_b = T2_a->parent()->lchild();
				bool multi_node = false;
				if (T2_a->parent()->get_children().size() > 2)
					multi_node = true;

			if (CUT_ONE_B) {
				if (T2_a->parent()->parent() == T2_c->parent()
					&& T2_c->parent() != NULL && !cut_b_only)
					cut_b_only=true;
					cob = true;
			}
			else if (CUT_ONE_AB) {
				if (T2_a->parent()->parent() == T2_c->parent()
					&& T2_c->parent() != NULL)
					cut_ab_only=true;
			}
			if (CUT_TWO_B && !cut_b_only && T1_ac->parent() != NULL) {
				Node *T1_s = T1_ac->get_sibling();
				if (T1_s->is_leaf()) {
					Node *T2_l = T2_a->parent()->parent();
					// Note: is this too harsh? If T2_l is nonbinary then can we do cut_b_only_if_not_a_or_c ?
					if (T2_l != NULL && T2_l->get_children().size() <= 2) {
						if (T2_c->parent() != NULL && T2_c->parent()->parent() == T2_l
								&& ((T2_a->parent()->get_children().size() <= 2
										&& T2_c->parent()->get_children().size() <= 2)
									|| T1_s->get_twin()->is_protected())) {
							if (T2_l->get_sibling() == T1_s->get_twin()) {
								cut_b_only=true;
							}
							else if (T2_l->parent() == NULL &&
									(T2->contains_rho() ||
									 T2->get_component(0) != T2_l)) {
								cut_b_only_if_not_a_or_c=true;
							}
						}
						else if ((T2_l = T2_l->parent()) != NULL
								&& T2_c->parent() == T2_l
								&& ((T2_a->parent()->get_children().size() <= 2
										&& T2_a->parent()->parent()->get_children().size() <= 2
										&& T2_l->get_children().size() <= 2)
									|| T1_s->get_twin()->is_protected())) {
							if (T2_l->get_sibling() == T1_s->get_twin()) {
								cut_b_only=true;
							}
							else if (T2_l->parent() == NULL &&
									(T2->contains_rho() ||
									 T2->get_component(0) != T2_l)) {
								cut_b_only_if_not_a_or_c=true;
							}
						}
					}
				}
			}
			if (REVERSE_CUT_ONE_B && (!cut_b_only || (cob && multi_node)) &&
					T1_ac->parent() != NULL) {
				Node *T1_s = T1_ac->get_sibling();
				if (T1_s->is_leaf()) {
					Node *T2_s = T1_s->get_twin();
					if (T2_s->parent() == T2_a->parent()) {
						cut_c_only=true;
						cut_b_only=false;
						cob=false;
					}
					else if (T2_s->parent() == T2_c->parent()) {
						if (T2_c->parent()->get_children().size() <= 2) {
							cut_a_only=true;
							cut_b_only=false;
							cob=false;
						}
						else {
							cut_a_or_merge_ac=true;
							cut_b_only=false;
							cob=false;
						}
					}
					else if (REVERSE_CUT_ONE_B_3
							// TODO: there is a chance for an additional optimization
							// here. If T2_s is not protected then we can cut c or
							// have to cut a (and s?)
//							&& T2_c->parent()->parent() != NULL
							// TODO: buggy? Do we also have to consider cutting b_1 through the other b's?
							&& T2_s->is_protected()
							&& T2_s->parent() != NULL
							&& T2_s->parent()->parent() == T2_a->parent()
							&& T2_s->parent()->get_children().size() <= 2) {
						//cut_c_only = true;
						cut_b_only=false;
						cob=false;
						if (!T2_a->is_protected()) {
							um.add_event(new ProtectEdge(T2_a));
							T2_a->protect_edge();
						}
					}
				}
				else if (REVERSE_CUT_ONE_B_2 && T2_c->parent() != NULL
						&& chain_match(T1_s, T2_c->get_sibling(), T2_a)) {
					cut_a_only = true;
					cut_b_only=false;
					cob=false;
				}
			}
			if (REVERSE_CUT_ONE_B_3 && (!cut_b_only || (cob && multi_node)) &&
					T1_ac->parent() != NULL && T1_ac->parent()->parent() != NULL) {
				Node *T1_s = T1_ac->parent()->get_sibling();
				Node *T2_s = T1_s->get_sibling();
				if (T1_s->is_leaf()
						&& T2_s->is_protected()
						&& T2_s->parent() != NULL
						&& T2_s->parent()->parent() == T2_a->parent()
						&& T2_s->parent()->get_children().size() <= 2) {
					//cut_c_only = true;
					cut_b_only=false;
					cob=false;
					if (!T2_a->is_protected()) {
						um.add_event(new ProtectEdge(T2_a));
						T2_a->protect_edge();
					}
				}
			}
			if (CUT_TWO_B_ROOT && cut_a_only == false && cut_c_only == false
					&& cut_b_only_if_not_a_or_c == true) {
				cut_b_only = true;
			}
/*			if (CUT_LOST) {
				if (T1_a->num_lost_children() > 0
						|| T2_a->num_lost_children() > 0) {
					cut_a_only = true;
					cut_b_only = false;
					cut_c_only = false;
					k++;
				}
				else if (T1_c->num_lost_children() > 0
						|| T2_c->num_lost_children() > 0) {
					cut_a_only = false;
					cut_b_only = false;
					cut_c_only = true;
					k++;
				}
				else if (T2_b->is_leaf()) {
					if (T2_b->num_lost_children() > 0
						|| T2_b->get_twin()->num_lost_children() > 0) {
					cut_a_only = false;
					cut_b_only = true;
					cut_c_only = false;
					k++;
					}
				}
			}
*/

			#ifdef DEBUG
					cout << "Case 3" << endl;
					cout << "\tT1: ";
					T1->print_components();
					cout << "\tT2: ";
					T2->print_components();
					cout << "\tK=" << k << endl;
					cout << "\tsibling pairs:";
					for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						cout << "  ";
						(*i).a->print_subtree_hlpr();
						cout << ",";
						(*i).c->print_subtree_hlpr();
					}
					cout << endl;
					cout << "\tprotected_stack:";
					for (list<Node *>::iterator i = protected_stack->begin(); i != protected_stack->end(); i++) {
						cout << "  ";
						(*i)->print_subtree_hlpr();
					}
					cout << endl;
					cout << "\tcut_a_only=" << cut_a_only << endl;
					cout << "\tcut_b_only=" << cut_b_only << endl;
					cout << "\tcut_c_only=" << cut_c_only << endl;
					cout << "\tT2_a " << T2_a->str() << " "
						<< T2_a->get_depth() << endl;
					cout << "\tT2_c " << T2_c->str() << " "
						<< T2_c->get_depth() << endl;
					cout << "\tT2_b " << T2_b->str_subtree() << " "
						<< T2_b->get_depth() << endl;
				#endif
			


				// copy elements
					/*
				Forest *T1_copy;
				Forest *T2_copy;
				list<Node *> *sibling_pairs_copy;
				Node *T1_a_copy;
				Node *T1_c_copy;
				Node *T2_a_copy;
				Node *T2_c_copy;
				*/
				//list<Node *> *singletons_copy = new list<Node *>();

				// make copies for the approx
				// be careful we do not kill real T1 and T2
				// ie use the copies
				if (BB && !cut_a_only && !cut_b_only && !cut_c_only) {
					list<Node *> *spairs;
					spairs = new list<Node *>();
					spairs->push_back(T1_c);
					spairs->push_back(T1_a);
					for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						spairs->push_back((*i).a);
						spairs->push_back((*i).c);
					}
					int approx_spr = rSPR_worse_3_approx_hlpr(T1, T2,
							singletons, spairs, NULL, NULL, false);
					delete spairs;
					#ifdef DEBUG
						cout << "\tT1: ";
						T1->print_components();
						cout << "\tT2: ";
						T2->print_components();
						cout << "approx =" << approx_spr << endl;
					#endif
					if (approx_spr  >  3*k){
						#ifdef DEBUG
							cout << "approx failed" << endl;
						#endif
						um.undo_all();
						return -1;
					}
					um.undo_to(undo_state);
				}

				if (!cut_a_only && !cut_c_only && !cut_b_only)
					same_component = T2_a->same_component(T2_c, lca_depth, path_length);
				if (cut_b_only)
					same_component = true;
				
			if (CLUSTER_REDUCTION && (MAX_CLUSTERS < 0 || NUM_CLUSTERS < MAX_CLUSTERS)) {
				// clean up singletons
				// TODO: this is duplication
				/*
			while(!singletons->empty()) {
				Node *T2_a = singletons->back();
				singletons->pop_back();
				// find twin in T1
				Node *T1_a = T2_a->get_twin();
				// if this is in the first component of T_2 then
				// it is not really a singleton.
				Node *T1_a_parent = T1_a->parent();
				if (T1_a_parent == NULL)
					continue;
				bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
				if (T2_a == T2->get_component(0)) {
					T1->add_rho();
					T2->add_rho();
					k--;
				}
	
				// cut the edge above T1_a
				T1_a->cut_parent();
				T1->add_component(T1_a);
				Node *node = T1_a_parent->contract();
				if (node != NULL && potential_new_sibling_pair && node->is_sibling_pair()){
					sibling_pairs->push_front(node->lchild());
					sibling_pairs->push_front(node->rchild());
				}
			}
			*/
	//			cout << "foo" << endl;
	//			cout << "foo2" << endl;
//				cout << "\tT1: ";
//				T1->print_components();
//				cout << "\tT2: ";
//				T2->print_components();
				sync_interior_twins_real(T1, T2);
				list<Node *> *cluster_points = find_cluster_points(T1, T2);
				//cluster_points->erase(++cluster_points->begin(),cluster_points->end());
	
				// TODO: could this be faster by using the approx to allocate
				// a certain amount of the k to different clusters?
				// TODO: write pseudocode for what we need
				// TODO: then implement it
				// NOTE: need to make a list of ClusterInstances and then
				// solve each.
				// TODO: where should we do this? Just before we would
				// normally branch?
	
//				cout << "k=" << k << endl;
//				cout << "cp=" << cluster_points->size() << endl;
				if (!cluster_points->empty()) {
					NUM_CLUSTERS++;
					sibling_pairs->clear();
#ifdef DEBUG_CLUSTERS
					cout << "CLUSTERS" << endl;
					for(int j = 0; j < 70; j++) {
						cout << "*";
					}
					cout << endl;
					for(list<Node *>::iterator i = cluster_points->begin();
							i != cluster_points->end(); i++) {
						cout << (*i)->str_subtree() << endl;
						cout << (*i)->get_twin()->str_subtree() << endl;
						for(int j = 0; j < 70; j++) {
							cout << "*";
						}
						cout << endl;
					}
					cout << endl;
#endif
	
					list<ClusterInstance> clusters =
						cluster_reduction(T1, T2, cluster_points);
	
					// TODO: make it so we don't need this?
					T1->unsync_interior();
					T2->unsync_interior();
					while(!clusters.empty()) {
						ClusterInstance cluster = clusters.front();
						clusters.pop_front();
						cluster.F1->unsync_interior();
						cluster.F2->unsync_interior();
#ifdef DEBUG_CLUSTERS
						cout << "CLUSTER_START" << endl;
						cout << &(*cluster.F1->get_component(0)) << endl;
						if (!clusters.empty())
							cout << &(*clusters.front().F1->get_component(0)) << endl;
						cout << "\tF1: ";
						cluster.F1->print_components();
						cout << "\tF2: ";
						cluster.F2->print_components();
						cout << "K=" << k << endl;
#endif
						int cluster_spr = -1;
						if (k >= 0) {
							// hack for clusters with no rho
//							cout << __LINE__ << endl;
//							cout << cluster.F2_cluster_node << endl;
//							cout << cluster.F2_has_component_zero << endl;
							if ((cluster.F2_cluster_node == NULL
										|| (cluster.F2_cluster_node->is_leaf()
												&& cluster.F2_cluster_node->parent() == NULL
												&& cluster.F2_cluster_node->
												get_num_clustered_children() <= 1
												&& (cluster.F2_cluster_node !=
														cluster.F2_cluster_node->get_forest()->
															get_component(0))))
										&& cluster.F2_has_component_zero == false) {
//							cout << __LINE__ << endl;
								cluster.F1->add_rho();
								cluster.F2->add_rho();
							}
							cluster_spr = rSPR_branch_and_bound_range(cluster.F1,
									cluster.F2, k);
							if (cluster_spr >= 0) {
//							cout << "cluster k=" << cluster_spr << endl;
//							cout << "\tF1: ";
//							cluster.F1->print_components();
//							cout << "\tF2: ";
//							cluster.F2->print_components();
								k -= cluster_spr;
							}
							else {
								k = -1;
							}
						}
						if (k > -1) {
//							cout << "\tF1: ";
//							T1->print_components();
//							cout << "\tF2: ";
//							T2->print_components();
								if (!cluster.is_original()) {
									int adjustment = cluster.join_cluster(T1, T2);
									k -= adjustment;
									delete cluster.F1;
									delete cluster.F2;
	
	//						cout << cluster.F1_cluster_node->str_subtree() << endl;
	//						cout << cluster.F2_cluster_node->str_subtree() << endl;
	//						Node *p = cluster.F1_cluster_node;
	//						while (p->parent() != NULL)
	//							p = p->parent();
	//						cout << &(*p) << endl;
	//						cout << p->str_subtree() << endl;
	
	
							//cout << "\tF1: ";
							//T1->print_components();
							//cout << "\tF2: ";
							//T2->print_components();
								}
//						else {
//							cout << "original" << endl;
//						}
							}
							else {
								if (!cluster.is_original()) {
									//if (cluster.F1_cluster_node != NULL)
									//	cluster.F1_cluster_node->contract();
									//if (cluster.F2_cluster_node != NULL)
									//	cluster.F2_cluster_node->contract();
									delete cluster.F1;
									delete cluster.F2;
								}
							}
					}
					delete cluster_points;
//					cout << "returning k=" << k << endl;
					NUM_CLUSTERS--;
					return k;
				}
				else {
					T1->unsync_interior();
					T2->unsync_interior();
					// HACK to allow only initial clusters
					// TODO: use UndoMachine in this clustering section
					// and update this clustering to not require copying
					NUM_CLUSTERS++;
				}
				delete cluster_points;
	
	//			cout << "done" << endl;
			}
				 // make copies for the branching
			/*
				copy_trees(&T1, &T2, &sibling_pairs, &T1_a, &T1_c, &T2_a, &T2_c,
						&T1_copy, &T2_copy, &sibling_pairs_copy,
					&T1_a_copy, &T1_c_copy, &T2_a_copy, &T2_c_copy);
					*/

				Node *node;
				Node *T2_ab = T2_a->parent();

				bool balanced = false;
				bool multi_b1 = T2_ab->get_children().size() > 2;
				bool multi_b2 = false;
				Node *T2_d = NULL;
				if (path_length == 4) {
					if (T2_a->parent()->parent() == T2_c->parent()->parent())
						balanced = true;
					if (balanced)
						multi_b2 = T2_c->parent()->get_children().size() > 2;
					else
						multi_b2 = T2_ab->parent()->get_children().size() > 2;
					if (balanced)
						T2_d = T2_c->get_sibling();
				}

				// cut T2_a
//				if (cut_b_only == false && T2_a->is_protected())
//					cout << "protected k=" << k << endl;
				if (cut_b_only == false && cut_c_only == false &&
						!T2_a->is_protected()
						&& (T2_a->parent()->parent() != NULL
								|| (T2_a->parent() == T2->get_component(0)
										&& !T2->contains_rho())
								|| !T2_b->is_protected()
								|| T2_a->parent()->get_children().size() > 2)) {// &&
//						(!T2_a->parent()->is_protected() ||
//							T2_a->parent()->get_children().size() > 2)) { }
					um.add_event(new CutParent(T2_a));
					T2_a->cut_parent();
					ContractEvent(&um, T2_ab);
					node = T2_ab->contract();
					if (node != NULL && node->is_singleton() &&
							node != T2->get_component(0))
						singletons->push_back(node);
					um.add_event(new AddComponent(T2));
					T2->add_component(T2_a);
					singletons->push_back(T2_a);

					// also require !cut_a_only ?
//					if (EDGE_PROTECTION_TWO_B && T2_c->is_protected() && !cut_a_only) {
//				}

					if (EDGE_PROTECTION_TWO_B && T2_c->is_protected() && !cut_a_only){
						if (path_length == 4) {
							if (!multi_b1 && !multi_b2 && !T2_b->is_protected()) {
								um.add_event(new ProtectEdge(T2_b));
								T2_b->protect_edge();
							}
							if (!multi_b2 && !multi_b1) {
								Node *T2_b2 = T2_b->parent()->get_sibling();
								if (balanced)
									T2_b2 = T2_d;
								if (!T2_b2->is_protected()) {
									um.add_event(new ProtectEdge(T2_b2));
									T2_b2->protect_edge();
								}
						}
						}
					}
					if (cut_a_only) {
						answer_a =
							rSPR_branch_and_bound_hlpr(T1, T2, k-1, sibling_pairs,
									singletons, false, AFs, protected_stack, num_ties, T1_c, T1_c->get_sibling());
					}
					else {
						answer_a =
							rSPR_branch_and_bound_hlpr(T1, T2, k-1, sibling_pairs,
									singletons, false, AFs, protected_stack, num_ties);
					}
				}
				best_k = answer_a;
				best_T1 = T1;
				best_T2 = T2;

				um.undo_to(undo_state);
/*				#ifdef DEBUG
					cout << "Case 3 CHECK" << endl;
					cout << "\tT1: ";
					T1->print_components();
					cout << "\tT2: ";
					T2->print_components();
					cout << "\tK=" << k << endl;
					cout << "\tsibling pairs:";
					for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						cout << "  ";
						(*i)->print_subtree_hlpr();
					}
					cout << endl;
					cout << "\tcut_b_only=" << cut_b_only << endl;
					cout << "\tT2_a " << T2_a->str() << " "
						<< T2_a->get_depth() << endl;
					cout << "\tT2_c " << T2_c->str() << " "
						<< T2_c->get_depth() << endl;
					cout << "\tT2_b " << T2_b->str_subtree() << " "
						<< T2_b->get_depth() << endl;
				#endif
*/

				//load the copy
				/*
				T1 = T1_copy;
				T2 = T2_copy;
				T1_a = T1_a_copy;
				T1_c = T1_c_copy;
				T2_a = T2_a_copy;
				T2_c = T2_c_copy;
				sibling_pairs = sibling_pairs_copy;
				singletons = new list<Node *>();
				*/


				// make copies for the branching
				/*
				copy_trees(&T1, &T2, &sibling_pairs, &T1_a, &T1_c, &T2_a, &T2_c,
						&T1_copy, &T2_copy, &sibling_pairs_copy,
						&T1_a_copy, &T1_c_copy, &T2_a_copy, &T2_c_copy);
						*/


				// get T2_b
				T2_ab = T2_a->parent();
				T2_b = T2_ab->rchild();

				if (T2_b == T2_a)
					T2_b = T2_ab->lchild();


//				if ((!CUT_AC_SEPARATE_COMPONENTS ||
//							T2_a->find_root() == T2_c->find_root())
//						&& (!multi_node && T2_b->is_protected()))
//					cout << "protected k=" << k << endl;

				// BUG: we need to cut b or c when b is a multifurcating COB
				// TODO: if the LCA of B_1's leaves in T1 is not an ancestor of 
				// a then this is safe
				if (multi_node && cob) {
					vector<Node *> B_1 = T2_a->parent()->find_leaves();
					LCA T1_LCA = LCA(T1->get_component(0));
					Node *B_1_lca = NULL;
					for(int i = 0; i < B_1.size(); i++) {
						if (B_1[i] == T2_a) {
							continue;
						}
						if (B_1_lca == NULL) {
							B_1_lca = B_1[i]->get_twin();
						}
						else {
							B_1_lca = T1_LCA.get_lca(B_1_lca, B_1[i]->get_twin());
						}
					}
					Node *T1_a_ancestor = T1_a;
					while(T1_a_ancestor != NULL && T1_a_ancestor != B_1_lca) {
						T1_a_ancestor = T1_a_ancestor->parent();
					}
					if (T1_a_ancestor != NULL) {
						#ifdef DEBUG
							cout << "\tmultifurcating cut_b_only so will also cut c" << endl;
						#endif
						cut_b_only=false;
					}
				}

				// cut T2_b
				if ((!CUT_AC_SEPARATE_COMPONENTS || same_component)
//						&& ((!T2_b->parent()->is_protected()
								&& (((multi_node || !T2_b->is_protected())))
						&& (!ABORT_AT_FIRST_SOLUTION || best_k < 0
							|| !PREFER_RHO || !AFs->front().first.contains_rho() )
						&& !cut_a_only && !cut_c_only
						&& (T2_a->parent()->parent() != NULL
								|| !T2_a->is_protected()
								|| (T2_a->parent() == T2->get_component(0)
										&& !T2->contains_rho()))) {
					if (multi_node) {
						um.add_event(new ChangeEdgePreInterval(T2_a));
						T2_a->copy_edge_pre_interval(T2_ab);
						um.add_event(new CutParent(T2_a));
						T2_a->cut_parent();
						um.add_event(new ChangeEdgePreInterval(T2_ab));
						T2_ab->set_edge_pre_start(-1);
						T2_ab->set_edge_pre_end(-1);
						Node *T2_ab_parent = T2_ab->parent();
						if (T2_ab_parent != NULL) {
							um.add_event(new CutParent(T2_ab));
							T2_ab->cut_parent();
							um.add_event(new AddChild(T2_a));
							T2_ab_parent->add_child(T2_a);
							um.add_event(new AddComponent(T2));
							T2->add_component(T2_ab);
						}
						else {
							if (T2->get_component(0) == T2_ab) {
								um.add_event(new AddComponentToFront(T2));
								T2->add_component(0, T2_a);
							}
							else {
								um.add_event(new AddComponent(T2));
								T2->add_component(T2_a);
								singletons->push_back(T2_a);
							}
						}
					}
					else {
						um.add_event(new CutParent(T2_b));
						T2_b->cut_parent();
						ContractEvent(&um, T2_ab);
						node = T2_ab->contract();
						if (node != NULL && node->is_singleton()
								&& node != T2->get_component(0))
								singletons->push_back(node);
						um.add_event(new AddComponent(T2));
						T2->add_component(T2_b);
						if (T2_b->is_leaf())
							singletons->push_back(T2_b);
					}
				add_sibling_pair(sibling_pairs, T1_a, T1_c,
						&um);

					// TODO: check carefully

					if (cut_a_or_merge_ac) {
						if (!T2_a->is_protected()) {
							um.add_event(new ProtectEdge(T2_a));
							T2_a->protect_edge();
							um.add_event(new ListPushBack(protected_stack));
							protected_stack->push_back(T2_a);
						}
						if (!T2_c->is_protected()) {
							um.add_event(new ProtectEdge(T2_c));
							T2_c->protect_edge();
						}
					}

					if (CUT_ALL_B) {
						answer_b =
							rSPR_branch_and_bound_hlpr(T1, T2, k-1,
									sibling_pairs, singletons, true, AFs, protected_stack,
									num_ties, T1_a, T1_c);
					}
					else {
						answer_b =
							rSPR_branch_and_bound_hlpr(T1, T2, k-1,
									sibling_pairs, singletons, false, AFs, protected_stack,
									num_ties, T1_a, T1_c);
					}
				}
				if (answer_b > best_k
						|| (answer_b == best_k
							&& PREFER_RHO
							&& T2->contains_rho() )) {
					best_k = answer_b;
					//swap(&best_T1, &T1);
					//swap(&best_T2, &T2);
				}

				um.undo_to(undo_state);

				/*
				delete T1;
				delete T2;
				delete sibling_pairs;
				delete singletons;
				*/


				// load the copy
				/*
				T1 = T1_copy;
				T2 = T2_copy;
				T1_a = T1_a_copy;
				T1_c = T1_c_copy;
				T2_a = T2_a_copy;
				T2_c = T2_c_copy;
				sibling_pairs = sibling_pairs_copy;
				singletons = new list<Node *>();
				*/
//				if (T2_c->is_protected())
//					cout << "protected k=" << k << endl;
				if (!T2_c->is_protected() &&
						!cut_a_or_merge_ac &&
	//					(T2_c->parent() == NULL || !T2_c->parent()->is_protected() ||
	//						T2_c->parent()->get_children().size() > 2) &&
						(!ABORT_AT_FIRST_SOLUTION || best_k < 0
							|| !PREFER_RHO || !AFs->front().first.contains_rho() )
						&& cut_b_only == false && cut_ab_only == false
						&& cut_a_only == false
						// TODO: do we allow this if T2_c has no parent?
						// it has to be under rho, right?
						&& (T2_c->parent() == NULL
								|| T2_c->parent()->parent() != NULL
								|| (T2_c->parent() == T2->get_component(0)
										&& !T2->contains_rho())
								|| !T2_c->get_sibling()->is_protected()
								|| T2_c->parent()->get_children().size() > 2)) {// &&


					if (T2_c->parent() != NULL) {
						Node *T2_c_parent = T2_c->parent();
						um.add_event(new CutParent(T2_c));
						T2_c->cut_parent();
						ContractEvent(&um, T2_c_parent);
						node = T2_c_parent->contract();
						if (node != NULL && node->is_singleton()
								&& node != T2->get_component(0))
							singletons->push_back(node);
						um.add_event(new AddComponent(T2));
						T2->add_component(T2_c);
					}
					else {
						// don't decrease k
						k++;
					}
					if (EDGE_PROTECTION && !cut_c_only) {
						if (!T2_a->is_protected()) {
							um.add_event(new ProtectEdge(T2_a));
							T2_a->protect_edge();
//							if (DEEPEST_PROTECTED_ORDER && !cut_c_only) {
							if (DEEPEST_PROTECTED_ORDER) {
								um.add_event(new ListPushBack(protected_stack));
								protected_stack->push_back(T2_a);
							}
							// TODO: add to protected list
						}
//						if (EDGE_PROTECTION_TWO_B && !cut_c_only) {
//					}
						// TODO: problem here :(
						if (EDGE_PROTECTION_TWO_B) {
							if (path_length == 4) {
								if (!multi_b1 && !multi_b2 && !T2_b->is_protected()) {
									um.add_event(new ProtectEdge(T2_b));
									T2_b->protect_edge();
								}
								if (!multi_b2 && !multi_b1) {
									Node *T2_b2 = T2_b->parent()->get_sibling();
									if (balanced)
										T2_b2 = T2_d;
									if (!T2_b2->is_protected()) {
										um.add_event(new ProtectEdge(T2_b2));
										T2_b2->protect_edge();
									}
								}
							}
						}
						if (path_length == 5) 
							T2_a->set_max_merge_depth(lca_depth);
					}
						singletons->push_back(T2_c);
						if (cut_c_only) {
							answer_c =
								rSPR_branch_and_bound_hlpr(T1, T2, k-1, sibling_pairs,
										singletons, false, AFs, protected_stack, num_ties, T1_a, T1_a->get_sibling());
						}
						else {
							answer_c =
								rSPR_branch_and_bound_hlpr(T1, T2, k-1, sibling_pairs,
										singletons, false, AFs, protected_stack, num_ties);
						}
						if (answer_c > best_k
									|| (answer_c == best_k
									&& PREFER_RHO
									&& T2->contains_rho() )) {
							best_k = answer_c;
							//swap(&best_T1, &T1);
							//swap(&best_T2, &T2);
						}
				}
				/*
				delete T1;
				delete T2;
				delete sibling_pairs;
				delete singletons;
				*/

				um.undo_to(undo_state);

				//T1 = best_T1;
				//T2 = best_T2;

#ifdef DEBUG_UNDO
		 while(um.num_events() > 0) {
				cout << "Undo step " << um.num_events() << endl;
				cout << "T1: ";
				T1->print_components();
				cout << "T2: ";
				T2->print_components();
					cout << "sibling pairs:";
					for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						cout << "  ";
						(*i).a->print_subtree_hlpr();
						cout << ",";
						(*i).c->print_subtree_hlpr();
					}
					cout << endl;
			 um.undo();
			cout << endl;
		 }
#else
		 um.undo_all();
#endif
				singletons->clear();
				return best_k;
			}
			cut_b_only = false;
		}
	}

	if (k >= 0) {
		if (PREFER_RHO && !AFs->empty() && !AFs->front().first.contains_rho() && T1->contains_rho()) {
			if (!ALL_MAFS)
				AFs->clear();
			AFs->push_front(make_pair(Forest(T1),Forest(T2)));
			*num_ties = 2;
		}
		else if (ALL_MAFS || AFs->empty()) {
			AFs->push_back(make_pair(Forest(T1),Forest(T2)));
		}
		else if (!PREFER_RHO || AFs->front().first.contains_rho() == T1->contains_rho()) {
			if (rand() < RAND_MAX/ *num_ties) {
				AFs->clear();
				AFs->push_back(make_pair(Forest(T1),Forest(T2)));
			}
			(*num_ties)++;
		}
	}

#ifdef DEBUG_UNDO
		 while(um.num_events() > 0) {
				cout << "Undo step " << um.num_events() << endl;
				cout << "T1: ";
				T1->print_components();
				cout << "T2: ";
				T2->print_components();
					cout << "sibling pairs:";
					for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						cout << "  ";
						(*i).a->print_subtree_hlpr();
						cout << ",";
						(*i).c->print_subtree_hlpr();
					}
					cout << endl;
			 um.undo();
			cout << endl;
		 }
#else
		 um.undo_all();
#endif

	return k;
}

int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map) {
	return rSPR_branch_and_bound_simple_clustering(T1,T2, verbose, label_map, reverse_label_map, -1, -1, NULL, NULL);
}

int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map, int min_k, int max_k) {
	return rSPR_branch_and_bound_simple_clustering(T1,T2, verbose, label_map, reverse_label_map, min_k, max_k, NULL, NULL);
}

int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map, int min_k, int max_k, Forest **out_F1, Forest **out_F2) {
	bool do_cluster = true;
	if (max_k > MAX_SPR)
		max_k = MAX_SPR;
	else if (max_k == -1)
		max_k = INT_MAX;

	if (T2->get_preorder_number() == -1) {
	  T2->preorder_number();
	}

	ClusterForest F1 = ClusterForest(T1);
	ClusterForest F2 = ClusterForest(T2);
	Forest F3 = Forest(F1);
	Forest F4 = Forest(F2);


//	bool old_rho = PREFER_RHO;
	PREFER_RHO = true;
	if (verbose) {
		cout << "T1: ";
		F1.print_components();
		cout << "T2: ";
		F2.print_components();
	}
	// h25: use the local singleton-greedy distance to decide whether the
	// complete instance is large enough to benefit from Whidden clustering.
	// Fall back to Whidden's approximation only if the adapter cannot parse it.
	int full_approx_spr = H25_GLOBAL_GREEDY_DISTANCE >= 0
		? H25_GLOBAL_GREEDY_DISTANCE
		: h25_local_greedy_distance(&F3, &F4);
	if (full_approx_spr < 0) {
		if (MULTIFURCATING) {
		  full_approx_spr = rSPR_worse_3_mult_approx(&F3, &F4);
		}
		else {
		  full_approx_spr = rSPR_worse_3_approx(&F3, &F4);
		}
	}
#ifdef H35_WHIDDEN_TRACE
	cerr << "h25 global_greedy_distance=" << full_approx_spr
		 << " cluster_tune=" << CLUSTER_TUNE << endl;
#endif
	if (full_approx_spr < CLUSTER_TUNE) {
		do_cluster = false;
	}
	if (verbose) {

		cout << "approx F1: ";
		F3.print_components();
		cout << "approx F2: ";
		F4.print_components();
		// what the AF shows
		cout << "approx drSPR=" << F4.num_components()-1 << endl;
		/* what we use to get the lower bound: 3 * the number of cutting rounds in
			 the approx algorithm
		*/
		//cout << "approx drSPR=" << full_approx_spr << endl;
		cout << "\n";
	}
//	if (F1.get_component(0)->get_preorder_number() == -1)
//		F1.get_component(0)->preorder_number();
//	if (F2.get_component(0)->get_preorder_number() == -1)
//		F2.get_component(0)->preorder_number();

	if (!sync_twins(&F1, &F2))
		return 0;
	if (F1.get_component(0)->is_leaf())
		return 0;
	if (F1.get_component(0)->get_preorder_number() == -1) {
		F1.get_component(0)->preorder_number();
		F2.get_component(0)->preorder_number();
	}
	int loss = 0;
	list<Node *> *cluster_points;
	if (F1.get_component(0)->get_edge_pre_start() == -1) {
		F1.get_component(0)->edge_preorder_interval();
		F2.get_component(0)->edge_preorder_interval();
	}
	if (LEAF_REDUCTION2) {
	  
	  if (MULTIFURCATING) {
	    reduction_leaf_mult(&F1, &F2);
	  }
	  else {	    
	    reduction_leaf(&F1, &F2);
	    }
//		F1.get_component(0)->preorder_number();
//		F2.get_component(0)->preorder_number();
//		F1.get_component(0)->edge_preorder_interval();
//		F2.get_component(0)->edge_preorder_interval();
	}
	if (COUNT_LOSSES) {
		loss += F1.get_component(0)->count_lost_subtree();
		loss += F2.get_component(0)->count_lost_subtree();
	}
		//F1.print_components();
		//F2.print_components();
	if (do_cluster) {
		sync_interior_twins(&F1, &F2);
		cluster_points = find_cluster_points(&F1, &F2);
		//	list<Node *> *cluster_points = new list<Node *>();
		for(list<Node *>::iterator i = cluster_points->begin();
				i != cluster_points->end(); i++) {
			string cluster_name = "X";
			/*
			if (verbose) {
					stringstream ss;
					ss << F1.size();
					cluster_name += ss.str();
					//int num_labels = label_map.size();
					//label_map.insert(make_pair(cluster_name,num_labels));
					//reverse_label_map.insert(
					//		make_pair(num_labels,cluster_name));
					//ss.str("");
					//ss << num_labels;
					//cluster_name = ss.str();
			}
			*/
	
			Node *n = *i;
			if (n->parent()->parent() == NULL
					&& n->get_sibling() != NULL &&
					n->get_sibling()->get_name() == "X")
				continue;
			Node *n_parent = n->parent();
			Node *twin = n->get_twin();
			Node *twin_parent = twin->parent();
	
			if (twin_parent == NULL)
				continue;
	
			F1.add_cluster(n,cluster_name);
	
			F2.add_cluster(twin,cluster_name);
	
			Node *n_cluster =
					F1.get_cluster_node(F1.num_clusters()-1);
			Node *twin_cluster =
					F2.get_cluster_node(F2.num_clusters()-1);
			n_cluster->set_twin(twin_cluster);
			twin_cluster->set_twin(n_cluster);
	
		}
		if (verbose)
			cout << endl << "CLUSTERS" << endl;

	}

	// component 0 needs to be done last
	F1.add_component(F1.get_component(0));
	F2.add_component(F2.get_component(0));

	int k;
	int num_clusters = F1.num_components();
	int total_k = 0;

	if(SHOW_CLUSTERS){
		cout << "Clusters start" << endl;
		for(int i = 1; i < num_clusters; i++) {
			Forest f1 = Forest(F1.get_component(i));
			f1.print_components();
		}
		cout << "Clusters end" << endl;
	}


	for(int i = 1; i < num_clusters; i++) {
		if (i == num_clusters - 1) {
			PREFER_RHO = false;
		}
		int exact_spr = -1;
		//vector<Node *> comps = vector<Node *>();
		//comps.push_back(F1.get_component(i));
		Forest f1 = Forest(F1.get_component(i));
		//Forest f1 = Forest(comps);
		//comps.clear();

		//comps.push_back(F1.get_component(i));
		Forest f2 = Forest(F2.get_component(i));
		//Forest f2 = Forest(comps);
		//comps.clear();
		Forest f1a = Forest(f1);
		Forest f2a = Forest(f2);
		Forest *f1_cluster;
		Forest *f2_cluster;

		if (verbose) {
			cout << "C" << i << "_1: ";
			f1.print_components();
			cout << "C" << i << "_2: ";
			f2.print_components();
		}
		// h25: replace Whidden's initial per-cluster 3-approximation with our
		// singleton greedy. Its result is a heuristic UB, so start 12% below
		// it and search upward. Whidden's restricted split-approx remains
		// unchanged later in this loop.
#ifdef WC1_CAPTURE_WHIDDEN_COMPONENTS
		Forest *wc1_greedy_first = NULL;
		Forest *wc1_greedy_second = NULL;
		int greedy_ub = h25_local_greedy_distance(
				&f1, &f2, &wc1_greedy_first, &wc1_greedy_second);
		if (wc1_greedy_first != NULL) {
			wc1_capture_whidden_forest(wc1_greedy_first,
					reverse_label_map, "cluster_greedy", i, -1);
			delete wc1_greedy_first;
			delete wc1_greedy_second;
		}
#else
		int greedy_ub = h25_local_greedy_distance(&f1, &f2);
#endif
		int approx_spr = greedy_ub;
		if (greedy_ub < 0) {
			if (MULTIFURCATING) {
			  approx_spr = rSPR_worse_3_mult_approx(&f1a, &f2a);
			}
			else {
			  approx_spr = rSPR_worse_3_approx(&f1a, &f2a);
			}
		}
		if (verbose) {
			cout << "cluster approx drSPR=" << f2a.num_components()-1 << endl;
			//cout << "cluster approx drSPR=" << approx_spr << endl;

			cout << endl;
		}

		int min_spr =
			greedy_ub >= 0 ? (greedy_ub * 88) / 100 : approx_spr / 3;
#ifdef H35_WHIDDEN_TRACE
		cerr << "h25 cluster=" << i
			 << " greedy_ub=" << greedy_ub
			 << " start=" << min_spr
			 << " split_threshold=" << SPLIT_APPROX_THRESHOLD << endl;
#endif
		if (min_spr < MIN_SPR - total_k)
			min_spr = MIN_SPR - total_k;
		int total_split_k = 0;

		bool done_cluster = false;
		bool done_split = false;

		double tree_fraction = INITIAL_TREE_FRACTION;

		if (min_spr < min_k)
			min_spr = min_k;

		while(!done_cluster) {
			done_cluster = true;

			for(k = min_spr - total_split_k; true; k++) {
				if (k < 0)
					k = 0;
				if (SPLIT_APPROX && !done_split && k >= SPLIT_APPROX_THRESHOLD) {
					done_cluster = false;
					break;
				}
				Forest f1t = Forest(f1);
//				Forest f1t = f1;
				Forest f2t = Forest(f2);
//				Forest f2t = f2;
				f1t.unsync();
				f2t.unsync();
				exact_spr = -1;
				if (verbose) {
					cout << k << " ";
  				cout.flush();
				}
				if (k + total_k <= max_k && k <= CLUSTER_MAX_SPR) {
					if (f1t.get_component(0)->get_name() == DEAD_COMPONENT) {
						f1t.add_rho();
						f2t.add_rho();
					}
					if (MULTIFURCATING) {
					  exact_spr = rSPR_branch_and_bound_mult(&f1t, &f2t, k);
					}
					else {
					  exact_spr = rSPR_branch_and_bound(&f1t, &f2t, k);
					}
				}
				if (exact_spr >= 0 || k + total_k > max_k ||
						k > CLUSTER_MAX_SPR) {
					if (k > CLUSTER_MAX_SPR) {
						f1t.swap(&f1a);
						f2t.swap(&f2a);
//						cout << "foo" << endl;
					}
						if (exact_spr >= 0) {
							exact_spr += total_split_k;
						if (verbose) {
	  					cout << endl;
	  					cout << "F" << i << "_1: ";
	  					f1t.print_components();
	  					cout << "F" << i << "_2: ";
	  					f2t.print_components();
	  					cout << "cluster exact drSPR=" << exact_spr << endl;
	  					cout << endl;
						}
#ifdef WC1_CAPTURE_WHIDDEN_COMPONENTS
						wc1_capture_whidden_forest(&f1t, reverse_label_map,
								"cluster_exact", i, -1);
#endif
						total_k += exact_spr;
					}
					else {
						// TODO: don't just the MAX_SPR here
						// incorporate extra information
						// toggle?
						if (verbose) {
							cout << "cluster exact drSPR=?  " << "k=" << k << " too large"
								<< endl;
							cout << "\n";
						}
						if (false && k > CLUSTER_MAX_SPR) {
							// TODO: this should be an approx of the remaining forest
//							total_k += approx_spr;
						}
						else if (CLAMP) {
							total_k = max_k;
						}
						else {
							Forest f1a = Forest(f1);
							Forest f2a = Forest(f2);
							
							int approx_spr;
							if (MULTIFURCATING) {
							  approx_spr = rSPR_worse_3_mult_approx(&f1a, &f2a);
							}
							else{
							  approx_spr = rSPR_worse_3_approx(&f1a, &f2a);
							}
								//total_k += min_spr;
								total_k += approx_spr / 3;
						}
						}
#ifdef WC4_ENABLE_PREJOIN_CLUSTER_ILP
						{
							Forest *wc4_cluster_first = NULL;
							Forest *wc4_cluster_second = NULL;
							int wc4_before_components = f1t.num_components();
							if (wc4_try_cluster_ilp(&f1t, &f2t,
									reverse_label_map, i, &F1, &F2, i,
									i >= num_clusters - 1,
									&wc4_cluster_first, &wc4_cluster_second)) {
								int wc4_after_components =
									wc4_cluster_first->num_components();
								f1t.swap(wc4_cluster_first);
								f2t.swap(wc4_cluster_second);
								delete wc4_cluster_first;
								delete wc4_cluster_second;
								int wc4_gain =
									wc4_before_components - wc4_after_components;
								if (wc4_gain > 0) {
									total_k -= wc4_gain;
									if (total_k < 0)
										total_k = 0;
								}
								cerr << "wc4_cluster_prejoin_commit cluster=" << i
									 << " before=" << wc4_before_components
									 << " after=" << wc4_after_components
									 << " total_k=" << total_k << endl;
							}
						}
#endif
						if ( i < num_clusters - 1) {
							F1.join_cluster(i,&f1t);
							F2.join_cluster(i,&f2t);
					}
					else {
							F1.join_cluster(&f1t);
							F2.join_cluster(&f2t);
						}
						break;
				}
			}
			done_split = done_cluster;
			int num_splits = 0;
			while (SPLIT_APPROX && !done_split) {
				//IN_SPLIT_APPROX = true;
				Node *original_split_node = find_subtree_of_approx_distance(
						f1.get_component(0), &f1, &f2, SPLIT_APPROX_THRESHOLD*2);
#ifdef WC2_ALT_SPLIT_SELECT
				original_split_node = wc2_select_split_candidate(
						original_split_node, &f1, &f2,
						SPLIT_APPROX_THRESHOLD * 2, i, num_splits);
#endif
#ifdef WC3_CAPTURE_ALT_SPLIT_TRIALS
				wc3_capture_alt_split_trials(original_split_node, &f1, &f2,
						SPLIT_APPROX_THRESHOLD * 2, reverse_label_map,
						i, num_splits);
#endif
#ifdef WC4_CLUSTER_ILP_SPLIT
				bool wc4_committed = false;
				if (!(original_split_node == f1.get_component(0) &&
						num_splits > 0)) {
					wc4_committed = wc4_try_split_ilp_commit(
							original_split_node, &f1, &f2,
							reverse_label_map, i, num_splits,
							total_split_k);
				}
				if (wc4_committed) {
					if (verbose)
						cout << "split_ilp" << endl;
				}
				else
#endif
				if (original_split_node == f1.get_component(0) &&
						num_splits > 0)
					done_split = true;
				else {
					Forest f1a = Forest(f1);
					Forest f2a = Forest(f2);
					Node *a_split_node =
					f1a.find_by_prenum(original_split_node->get_preorder_number());
					f1a.get_component(0)->disallow_siblings_subtree();
						a_split_node->allow_siblings_subtree();
//					if (a_split_node->lchild() != NULL)
//						a_split_node->lchild()->allow_siblings_subtree();
//					if (a_split_node->rchild() != NULL)
//						a_split_node->rchild()->allow_siblings_subtree();
					// something odd going on here
					int start = rSPR_worse_3_approx(a_split_node, &f1a, &f2a);
					if (start == INT_MAX)
						start = 0;
					start /= 3;
					int end = f1.get_component(0)->size();
#ifdef RR4_SPLIT_TRACE
					{
						vector<Node *> rr4_leaves =
							original_split_node->find_leaves();
						vector<int> rr4_leaf_ids;
						for (size_t rr4_i = 0; rr4_i < rr4_leaves.size();
								rr4_i++) {
							string rr4_name = rr4_leaves[rr4_i]->get_name();
							if (!rr4_name.empty() &&
									isdigit((unsigned char)rr4_name[0])) {
								rr4_leaf_ids.push_back(atoi(rr4_name.c_str()));
							}
						}
						sort(rr4_leaf_ids.begin(), rr4_leaf_ids.end());
#ifdef RR4_CAPTURE_SPLIT_REGIONS
						if (!rr4_leaf_ids.empty())
							RR4_CAPTURED_SPLIT_REGIONS.push_back(rr4_leaf_ids);
#endif
						cerr << "rr4_split_region cluster=" << i
							 << " split_no=" << num_splits
							 << " pre="
							 << original_split_node->get_preorder_number()
							 << " leaves=" << rr4_leaves.size()
							 << " numeric_leaves=" << rr4_leaf_ids.size()
							 << " start=" << start
							 << " end=" << end
							 << " ids=";
						for (size_t rr4_i = 0; rr4_i < rr4_leaf_ids.size();
								rr4_i++) {
							if (rr4_i) cerr << ",";
							cerr << rr4_leaf_ids[rr4_i];
						}
						cerr << endl;
					}
#endif
					for(k = start; true; k++) {
						// TODO: figure out the bug here
						if (k > end) {
							k = 0;
							done_split = true;
							break;
						}
				/*	if (k > SPLIT_APPROX_THRESHOLD) {
						k = 0;
						tree_fraction *= 0.75;
						if (verbose)
							cout << "tree_fraction: " << tree_fraction << endl;
						continue;
					}*/
						Forest f1s = Forest(f1);
						Forest f2s = Forest(f2);
						if (!sync_twins(&f1s, &f2s)) {
							k = 0;
							done_split = true;
							break;
						}
						if (verbose) {
							cout << k << " ";
		  				cout.flush();
						}
						Node *split_node = f1s.find_by_prenum(original_split_node->get_preorder_number());
						f1s.get_component(0)->disallow_siblings_subtree();
							split_node->allow_siblings_subtree();
//						if (split_node->lchild() != NULL)
//							split_node->lchild()->allow_siblings_subtree();
//						if (split_node->rchild() != NULL)
//							split_node->rchild()->allow_siblings_subtree();
							//f1s.get_component(0)->find_subtree_of_size(tree_fraction);
							set<SiblingPair > *sibling_pairs =
								find_sibling_pairs_set(split_node);
							list<Node *> singletons = f2s.find_singletons();
							list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();
							list<Node *> protected_stack = list<Node *>();

							int num_ties = 2;

							int split_k = rSPR_branch_and_bound_hlpr(&f1s, &f2s, k,
									sibling_pairs, &singletons, false, &AFs,
									&protected_stack, &num_ties);
							delete sibling_pairs;
							if (!AFs.empty()) {
								AFs.front().first.swap(&f1);
								AFs.front().second.swap(&f2);
#ifdef WC1_CAPTURE_WHIDDEN_COMPONENTS
								wc1_capture_whidden_forest(&f1, reverse_label_map,
										"split_commit", i, num_splits);
#endif
								f2.unprotect_edges();
								f1.get_component(0)->allow_siblings_subtree();
								AFs.clear();
								total_split_k += k - split_k;
#ifdef RR4_SPLIT_TRACE
								cerr << "rr4_split_commit cluster=" << i
									 << " split_no=" << num_splits
									 << " pre="
									 << original_split_node->get_preorder_number()
									 << " k=" << k
									 << " split_k=" << split_k
									 << " gained=" << (k - split_k)
									 << " total_split_k=" << total_split_k
									 << " f1_components=" << f1.num_components()
									 << " f2_components=" << f2.num_components()
									 << endl;
#endif
		//						if (k < SPLIT_APPROX_THRESHOLD * 0.75) {
		//							tree_fraction *= 2;
		//							if (tree_fraction > INITIAL_TREE_FRACTION)
		//								tree_fraction = INITIAL_TREE_FRACTION;
		//						}
								if (verbose)
									cout << "split_k: " << k << endl;
								break;
							}
					}
				}
				//IN_SPLIT_APPROX = false;
				num_splits++;
			}
	
			// TODO: approx again? seperate approxes ?
		}
	}

		if (F1.contains_rho()) {
			F1.get_component(0)->delete_tree();
			F2.get_component(0)->delete_tree();
			F1.erase_components(0, num_clusters);
			F2.erase_components(0, num_clusters);
		}
		else {
			F1.get_component(num_clusters)->delete_tree();
			F2.get_component(num_clusters)->delete_tree();
			F1.erase_components(1, num_clusters+1);
			F2.erase_components(1, num_clusters+1);
		}
		// fix hanging roots
		for(int i = 0; i < F1.num_components(); i++) {
			F1.get_component(i)->contract(true);
			F2.get_component(i)->contract(true);
		}
		if (verbose) {
			F1.numbers_to_labels(reverse_label_map);
			F2.numbers_to_labels(reverse_label_map);
			cout << "F1: ";
			F1.print_components();
			cout << "F2: ";
			F2.print_components();
			cout << "total exact drSPR=" << total_k << endl;
		}
		if (out_F1 != NULL)
			*out_F1 = new Forest(&F1);
			//F1.swap(out_F1);
		if (out_F2 != NULL)
			*out_F2 = new Forest(&F2);
//			F2.swap(out_F2);
		if (out_F1 != NULL && out_F2 != NULL) {
//			out_F1->resync();
			sync_twins(*out_F1, *out_F2);
	//		sync_interior_twins_real(out_F1, out_F2);
		}

	if (do_cluster) {
		delete cluster_points;
	}
//	PREFER_RHO = old_rho;
	total_k += loss;
/*	cout << "F1: ";
	for (int i = 0; i < F1.num_components(); i++) {
		if (i > 0)
			cout << " ";
		F1.get_component(i)->expand_contracted_nodes();
		cout << F1.get_component(i)->str_edge_pre_interval_subtree();
	}
	cout << endl;
	cout << "F2: ";
	for (int i = 0; i < F2.num_components(); i++) {
		if (i > 0)
			cout << " ";
		F2.get_component(i)->expand_contracted_nodes();
		cout << F2.get_component(i)->str_edge_pre_interval_subtree();
	}
	cout << endl << endl;
*/
	#ifdef DEBUG_CASE_COUNTER
	print_mult_case_count();
        #endif

	return total_k;
}

int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map) {
	Forest F1 = *T1;//Forest(T1);
	Forest F2 = *T2;//Forest(T2);
	Forest F3 = Forest(F1);
	Forest F4 = Forest(F2);

	bool do_cluster = true;

//	bool old_rho = PREFER_RHO;
	PREFER_RHO = true;
	if (verbose) {
		cout << "T1: ";
		F1.print_components();
		cout << "T2: ";
		F2.print_components();
	}

	int full_approx_spr;
	if (MULTIFURCATING) {
	  full_approx_spr = rSPR_worse_3_mult_approx(&F3, &F4);
	}
	else {
	  full_approx_spr = rSPR_worse_3_approx(&F3, &F4);
	}
	if (full_approx_spr <= CLUSTER_TUNE) {
		do_cluster = false;
	}
	if (verbose) {

		cout << "approx F1: ";
		F3.print_components();
		cout << "approx F2: ";
		F4.print_components();
		// what the AF shows
		cout << "approx drSPR=" << F4.num_components()-1 << endl;
		/* what we use to get the lower bound: 3 * the number of cutting rounds in
			 the approx algorithm
		*/
		//cout << "approx drSPR=" << full_approx_spr << endl;
		cout << "\n";
	}

	if (!sync_twins(&F1, &F2))
		return 0;
	if (F1.get_component(0)->is_leaf())
		return 0;
	list<Node *> *cluster_points;
	sync_interior_twins_real(&F1, &F2);
	if (do_cluster) {
		cluster_points = find_cluster_points(&F1, &F2);
	}

		int total_k = 0;

		if (do_cluster && !cluster_points->empty()) {
	
			list<ClusterInstance> clusters =
				cluster_reduction(&F1, &F2, cluster_points);
	
			F1.unsync_interior();
			F2.unsync_interior();
		int k = 0;
		int i = 0;
		while(!clusters.empty()) {
			i++;
			ClusterInstance cluster = clusters.front();
			clusters.pop_front();
			cluster.F1->unsync_interior();
			cluster.F2->unsync_interior();
			if (verbose) {
				cout << "C" << i << "_1: ";
				cluster.F1->print_components();
				cout << "C" << i << "_2: ";
				cluster.F2->print_components();
			}
			Forest f1 = Forest(cluster.F1);
			Forest f2 = Forest(cluster.F2);
			
			int min_spr;
			if (MULTIFURCATING) {
			  min_spr = rSPR_worse_3_mult_approx(&f1, &f2);
			}
			else {
			  min_spr = rSPR_worse_3_approx(&f1, &f2);
			}
			min_spr /= 3;

			if (verbose) {
				cout << "cluster approx drSPR=" << f1.num_components()-1 << endl;
				//cout << "cluster approx drSPR=" << min_spr << endl;
				cout << endl;
			}

			int cluster_spr = -1;
			k = MAX_SPR - total_k;
			if (k >= 0) {
				// hack for clusters with no rho
				if ((cluster.F2_cluster_node == NULL
							|| (cluster.F2_cluster_node->is_leaf()
									&& cluster.F2_cluster_node->parent() == NULL
									&& cluster.F2_cluster_node->
									get_num_clustered_children() <= 1
									&& (cluster.F2_cluster_node !=
											cluster.F2_cluster_node->get_forest()->
												get_component(0))))
							&& cluster.F2_has_component_zero == false) {
					cluster.F1->add_rho();
					cluster.F2->add_rho();
				}
				if (MULTIFURCATING) {
				  cluster_spr = rSPR_branch_and_bound_mult_range(cluster.F1,
						cluster.F2, min_spr, MAX_SPR - total_k);
				}
				else {
				  cluster_spr = rSPR_branch_and_bound_range(cluster.F1,
						cluster.F2, min_spr, MAX_SPR - total_k);
				}
				if (cluster_spr >= 0) {
					if (verbose) {
	  				cout << endl;
	  				cout << "F" << i << "_1: ";
	  				cluster.F1->print_components();
	  				cout << "F" << i << "_2: ";
	  				cluster.F2->print_components();
	  				cout << "cluster exact drSPR=" << cluster_spr << endl;
	  				cout << endl;
					}
					total_k += cluster_spr;
					k -= cluster_spr;
				}
				else {
					k = -1;
					if (verbose) {
						cout << "cluster exact drSPR=?  " << "k=" << k << " too large"
							<< endl;
						cout << "\n";
					}
					if (CLAMP) {
						total_k = MAX_SPR + 1;
					}
					else {
							total_k += min_spr;
					}
				}
			}
			if (k > -1) {
				if (!cluster.is_original()) {
					int adjustment = cluster.join_cluster(&F1, &F2);
					total_k += adjustment;
					delete cluster.F1;
					delete cluster.F2;
				}
			}
			else {
				if (!cluster.is_original()) {
					//if (cluster.F1_cluster_node != NULL)
					//	cluster.F1_cluster_node->contract();
					//if (cluster.F2_cluster_node != NULL)
					//	cluster.F2_cluster_node->contract();
					delete cluster.F1;
					delete cluster.F2;
				}
			}
		}
		delete cluster_points;
	}
	else {
		if (do_cluster) {
			delete cluster_points;
		}
		full_approx_spr /= 3;
		if (MULTIFURCATING) {
		  total_k = rSPR_branch_and_bound_mult_range(&F1, &F2, full_approx_spr, MAX_SPR);
		}
		else {
		  total_k = rSPR_branch_and_bound_range(&F1, &F2, full_approx_spr, MAX_SPR);
		}
		int i = 1;
		if (total_k < 0) {
			if (CLAMP) {
				total_k = MAX_SPR;
			}
			else {
		 	 total_k = full_approx_spr;
			}
		}

	}

	if (verbose) {
		F1.numbers_to_labels(reverse_label_map);
		F2.numbers_to_labels(reverse_label_map);
	 	cout << endl;
	 	cout << "F1: ";
	 	F1.print_components();
	 	cout << "F2: ";
	 	F2.print_components();
	 	cout << "total exact drSPR=" << total_k << endl;
	 	cout << endl;
	}
	return total_k;
}

/*Joel's part*/
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2, bool verbose){
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false, NULL, NULL);
}

int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2){
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false, NULL, NULL);
}

int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose) {
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false, NULL, NULL, -1, -1);
}

int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, int min_k, int max_k) {
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false, NULL, NULL, min_k, max_k);
}

int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2) {
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false);
}
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, Forest **out_F1, Forest **out_F2) {
	return rSPR_branch_and_bound_simple_clustering(T1,T2, false, NULL, NULL, -1, -1, out_F1, out_F2);
}

// T1 and T2 are assumed to already be synced
void reduction_leaf_mult(Forest *T1, Forest* T2) {
  
  list<Node *> *sibling_groups = T1->find_sibling_groups();
  while (!sibling_groups->empty()) {
    //Get a sibling group with identical siblings
    list<Node*>::reverse_iterator i = sibling_groups->rbegin();
    Node *T1_sibling_group = sibling_groups->back();
    list<list<Node*>> identical_sibling_groups;
    T1_sibling_group->find_identical_sibling_groups(&identical_sibling_groups);
    for (; i != sibling_groups->rend(); i++ ){
      (*i)->find_identical_sibling_groups(&identical_sibling_groups);
      if (identical_sibling_groups.size() > 0) {
	T1_sibling_group = (*i);
	break;
      }
    }

    // Checked all sibling groups, found none with identical groups in T2
    // Therefore there are no more contractions to be made

    if (identical_sibling_groups.size() == 0) {
      delete sibling_groups;
      return;
      }
    // Contract them
    else {
      list<list<Node *>>::iterator i;
      for (i = identical_sibling_groups.begin(); i != identical_sibling_groups.end(); i++) {
	list<Node *> T2_group = (*i);
	Node *T2_p = T2_group.front()->parent();
	Node *T1_group_new = T1_sibling_group->contract_twin_group(&T2_group);
	Node *T2_group_new = T2_p->contract_sibling_group(&T2_group);
	  
	T1_group_new->set_twin(T2_group_new);
	T2_group_new->set_twin(T1_group_new);			

	if (T1_sibling_group->parent() != NULL) {
	  //Check if the contraction made a new sibling group
	  T1_sibling_group->parent()->recalculate_non_leaf_children();
	  if (T1_sibling_group->parent()->is_sibling_group()) {
	    sibling_groups->push_front(T1_sibling_group->parent());
	  }
	}
	//Check if this contraction removed a sibling group
	if (!T1_sibling_group->is_sibling_group()) {
	  sibling_groups->remove(T1_sibling_group);
	}	  
      }
    }
  }
  delete sibling_groups;
}
  
// T1 and T2 are assumed to already be synced
void reduction_leaf(Forest *T1, Forest *T2) {
	reduction_leaf(T1, T2, NULL);
}

void reduction_leaf(Forest *T1, Forest *T2, UndoMachine *um) {
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	Node *T1_a;
	Node *T1_c;
	while (!sibling_pairs->empty()) {
		T1_a = sibling_pairs->front();
		sibling_pairs->pop_front();
		T1_c = sibling_pairs->front();
		sibling_pairs->pop_front();
		// shouldn't happen here
		if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
				continue;
		}
		Node *T2_a = T1_a->get_twin();
		Node *T2_c = T1_c->get_twin();
		if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
			Node *T1_ac = T1_a->parent();
			Node *T2_ac = T2_a->parent();
			T1_ac->contract_sibling_pair_undoable();
			Node *T2_ac_new = T2_ac->contract_sibling_pair_undoable(T2_a, T2_c);
			if (T2_ac_new != NULL && T2_ac_new != T2_ac) {
				T2_ac = T2_ac_new;
				T2_ac->contract_sibling_pair_undoable();
			}

			T1_ac->set_twin(T2_ac);
			T2_ac->set_twin(T1_ac);

			// check if T2_ac is a singleton
			// also shouldn't happen
//				if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))

			// check if T1_ac is part of a sibling pair
			if (T1_ac->parent() != NULL &&
					T1_ac->parent()->is_sibling_pair()) {
				sibling_pairs->push_back(T1_ac->parent()->lchild());
				sibling_pairs->push_back(T1_ac->parent()->rchild());
			}
			#ifdef DEBUG
				cout << "\tT1: ";
				T1->print_components();
				cout << "\tT2: ";
				T2->print_components();
			#endif
		}
	}
	delete sibling_pairs;
}

/* return true if T1_node matches the chain between T2_node and
	 T2_node_end
*/
bool chain_match(Node *T1_node, Node *T2_node, Node *T2_node_end) {
	Node *T1_pendant;
	Node *T2_pendant;
	bool pendant_found = false;
	if (T2_node->is_leaf())
		return false;
	// T1_node is a leaf
	if (T1_node->is_leaf()) {
		T1_pendant = T1_node;
		if (T1_pendant->get_twin() == T2_node->lchild()) {
			if (T2_node->rchild() == T2_node_end)
				return true;
		}
		else if (T1_pendant->get_twin() == T2_node->rchild()) {
			if (T2_node->lchild() == T2_node_end)
				return true;
		}
		return false;
	}
	// T1_pendant is T1_node->lchild()
	T1_pendant = T1_node->lchild();
	if (T1_pendant->is_leaf()) {
		T2_pendant = T2_node->lchild();
		if (T2_pendant->is_leaf() && T1_pendant->get_twin() == T2_pendant) {
			return chain_match(T1_pendant->get_sibling(),
					T2_pendant->get_sibling(), T2_node_end);
		}
		T2_pendant = T2_node->rchild();
		if (T2_pendant->is_leaf() && T1_pendant->get_twin() == T2_pendant) {
			return chain_match(T1_pendant->get_sibling(),
					T2_pendant->get_sibling(), T2_node_end);
		}
	}
	// T1_pendant is T1_node->rchild()
	if (T1_pendant->is_leaf()) {
		T2_pendant = T2_node->lchild();
		if (T2_pendant->is_leaf() && T1_pendant->get_twin() == T2_pendant) {
			return chain_match(T1_pendant->get_sibling(),
					T2_pendant->get_sibling(), T2_node_end);
		}
		T2_pendant = T2_node->rchild();
		if (T2_pendant->is_leaf() && T1_pendant->get_twin() == T2_pendant) {
			return chain_match(T1_pendant->get_sibling(),
					T2_pendant->get_sibling(), T2_node_end);
		}
	}
	return false;
}

int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees) {
	return rSPR_total_distance(T1, gene_trees, NULL);
}

int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees,
		vector<int> *original_scores) {
	int total = 0;
	MAIN_CALL = false;
	int end = gene_trees.size();
//	T1->preorder_number();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
//	for(int j = 0; j < 10; j++)
//	cout << "T1: " << T1->str_subtree() << endl;
	for(int i = 0; i < end; i++) {
			//		cout << i << endl;
	  cout << "Trying tree #" << i << " : " << gene_trees[i]->str_subtree() << endl;
	  int k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i], VERBOSE);

	  MULTIFURCATING = true;
		//MULT_4_BRANCH = true;
		int mult_k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i], VERBOSE);
		MULTIFURCATING = false;
		//MULT_4_BRANCH = false;
		
		if (k != mult_k) {
		  cout << "BINARY DOES NOT MATCH MULT" << endl;
		  cout << "T1: " << T1->str_subtree() << endl;;
		  cout << "BINARY k = " << k << " mult_k = " << mult_k << endl;
		  break;
		  }
		else {
		  cout << "\tMATCHES: k = " << k << endl;
		}
		
		//cout << "\t: k = " << k << endl;
//		k *= mylog2(gene_trees[i]->size());

		if (original_scores != NULL)
			(*original_scores)[i] = k;
		total += k;
//		cout << "T2: " << gene_trees[i]->str_subtree() << endl;
//		cout << " k: " << k << endl;
		if (FIND_RATE) {
			if (k > 0) {
				int size = gene_trees[i]->find_leaves().size();
//				cout << k << endl;
//				cout << size << endl;
				cout << "rate=" << (float)k / size << endl;
//				cout << T1->str_subtree() << endl;
//				cout << gene_trees[i]->str_subtree() << endl;
			}
		}
//		Forest F1 = Forest(T1);
//		Forest F2 = Forest(gene_trees[i]);
//		total += rSPR_branch_and_bound(&F1, &F2);
	}
	return total;
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees) {
	rSPR_pairwise_distance(T1, gene_trees, 0, gene_trees.size());
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, bool APPROX) {
	rSPR_pairwise_distance(T1, gene_trees, 0, gene_trees.size(), APPROX);
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end) {
	rSPR_pairwise_distance(T1, gene_trees, start, end, false);
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end, bool approx) {
	MAIN_CALL = false;
//	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int k;
		if (approx) {
			Forest F1 = Forest(T1);
			Forest F2 = Forest(gene_trees[i]);
			k = rSPR_worse_3_approx_distance_only(&F1, &F2)/3;
		}
		else {
			k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i]);
		}
		distances[i-start] = k;
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}


void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int max_spr) {
	rSPR_pairwise_distance(T1, gene_trees, max_spr, 0, (int)gene_trees.size());
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int max_spr, int start, int end) {
	MAIN_CALL = false;
//	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		Forest F1 = Forest(T1);
		Forest F2 = Forest(gene_trees[i]);
		int k = rSPR_branch_and_bound_range(&F1, &F2, 0, max_spr);
		distances[i-start] = k;
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	rSPR_pairwise_distance_unrooted(T1, gene_trees, 0, gene_trees.size());
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, bool approx) {
	rSPR_pairwise_distance_unrooted(T1, gene_trees, 0, gene_trees.size(), approx);
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end) {
	rSPR_pairwise_distance_unrooted(T1, gene_trees, 0, gene_trees.size(), false);
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end, bool approx) {
	MAIN_CALL = false;
	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int best_k = INT_MAX;
		Node *T2_copy = new Node(*(gene_trees[i]));
		vector<Node *> descendants = 
				T2_copy->find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			T2_copy->reroot(descendants[j]);
			T2_copy->set_depth(0);
			T2_copy->fix_depths();
			T2_copy->preorder_number();
	//				cout << i << "," << j << endl;
	//				cout << T1->str_subtree() << endl;
	//				cout << gene_trees[i]->str_subtree() << endl;
			int k;
			if (approx) {
				Forest F1 = Forest(T1);
				Forest F2 = Forest(T2_copy);
				k = rSPR_worse_3_approx(&F1, &F2) / 3;
			}
			else {
				k = rSPR_branch_and_bound_simple_clustering(T1, T2_copy, false);
			}
			if (k < best_k) {
				best_k = k;
			}
		}
		distances[i-start] = best_k;
		T2_copy->delete_tree();
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int max_spr) {
	rSPR_pairwise_distance_unrooted(T1, gene_trees, max_spr, 0, (int)gene_trees.size());
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int max_spr, int start, int end) {
	MAIN_CALL = false;
	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int best_k = -1;
		Node *T2_copy = new Node(*(gene_trees[i]));
		vector<Node *> descendants = 
				T2_copy->find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			T2_copy->reroot(descendants[j]);
			T2_copy->set_depth(0);
			T2_copy->fix_depths();
			T2_copy->preorder_number();
	//				cout << i << "," << j << endl;
	//				cout << T1->str_subtree() << endl;
	//				cout << gene_trees[i]->str_subtree() << endl;
			Forest F1 = Forest(T1);
			Forest F2 = Forest(T2_copy);
			int k = rSPR_branch_and_bound_range(&F1, &F2, 0, max_spr);
			if ((best_k == -1) || (k < best_k && k >= 0)) {
				best_k = k;
			}
		}
		distances[i-start] = best_k;
		T2_copy->delete_tree();
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

int rSPR_total_distance_precomputed(Node *T1, vector<Node *> &gene_trees,
		vector<int> *original_scores, vector<int> *new_original_scores, Node *old_T1) {
	int total = 0;
	MAIN_CALL = false;
	int end = gene_trees.size();
//	T1->preorder_number();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
	for(int i = 0; i < end; i++) {
		// check that the SPR move affects the projection of T1
		Forest F1 = Forest(T1);
		Forest F2 = Forest(gene_trees[i]);
		Forest F1_old = Forest(old_T1);
		sync_twins(&F1, &F2);
		sync_twins(&F1, &F1_old);
		int k = 0;
		if (original_scores == NULL
				|| rSPR_worse_3_approx(&F1, &F1_old) > 0) {
			k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i], VERBOSE);
		}
		else {
			k = (*original_scores)[i];
		}
		if (new_original_scores != NULL) {
			(*new_original_scores)[i] = k;
		}

		total += k;
	}
	return total;
}


int rf_total_distance(Node *T1, vector<Node *> &gene_trees) {
	int total = 0;
	int end = gene_trees.size();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
	for(int i = 0; i < end; i++) {
			//		cout << i << endl;
		int k = rf_distance(T1, gene_trees[i]);
		total += k;
	}
	return total;
}

int rf_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	int total = 0;
	int end = gene_trees.size();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
	for(int i = 0; i < end; i++) {
		int best_k = INT_MAX;
		Node T2_copy = Node(*(gene_trees[i]));
		vector<Node *> descendants = 
				T2_copy.find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			T2_copy.reroot(descendants[j]);
			T2_copy.set_depth(0);
			T2_copy.fix_depths();
			T2_copy.preorder_number();
			int k = rf_distance(T1, &T2_copy);
			if (k < best_k) {
				best_k = k;
			}
		}
		total += best_k;
	}
	return total;
}



void rf_pairwise_distance(Node *T1, vector<Node *> &gene_trees) {
	rf_pairwise_distance(T1, gene_trees, 0, gene_trees.size());
}

void rf_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end) {
	MAIN_CALL = false;
//	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int k = rf_distance(T1, gene_trees[i]);
		distances[i-start] = k;
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

void rf_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	rf_pairwise_distance_unrooted(T1, gene_trees, 0, gene_trees.size());
}

void rf_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end) {
	MAIN_CALL = false;
	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int best_k = INT_MAX;
		Node T2_copy = Node(*(gene_trees[i]));
		vector<Node *> descendants = 
				T2_copy.find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			T2_copy.reroot(descendants[j]);
			T2_copy.set_depth(0);
			T2_copy.fix_depths();
			T2_copy.preorder_number();
	//				cout << i << "," << j << endl;
	//				cout << T1->str_subtree() << endl;
	//				cout << gene_trees[i]->str_subtree() << endl;
			int k = rf_distance(T1, &T2_copy);
			if (k < best_k) {
				best_k = k;
			}
		}
		distances[i-start] = best_k;
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees, int threshold) {
	int total = 0;
	MAIN_CALL = false;
	int end = gene_trees.size();
	T1->preorder_number();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
	for(int i = 0; i < end; i++) {
		int k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i], VERBOSE);
//		k *= mylog2(gene_trees[i]->size());
		total += k;
//		if (total > threshold) {
//			break;
//		}
	}
	return total;
}

/*Joel's part*/
int rSPR_total_distance(Forest *T1, vector<Node *> &gene_trees){
	int total = 0;
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)
	for(int i = 0; i < gene_trees.size(); i++) {
		Forest T2 = Forest(gene_trees[i]);
		total += rSPR_branch_and_bound_simple_clustering(&T2, T1, VERBOSE);
	}
	return total;
}

int rSPR_total_approx_distance(Forest *T1, vector<Node *> &gene_trees) {
	int total = 0;
	#pragma omp parallel for reduction(+ : total)
	for(int i = 0; i < gene_trees.size(); i++) {
		Forest F1 = Forest(T1);
		Forest F2 = Forest(gene_trees[i]);
//		cout << i << endl;
//		cout << T1->str_subtree() << endl;
//		cout << gene_trees[i]->str_subtree() << endl;
		//total += rSPR_worse_3_approx(&F2, &F1)/3;
		total += rSPR_worse_3_approx(&F2, &F1)/3;
	}
	return total;
}

int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	return rSPR_total_distance_unrooted(T1, gene_trees, INT_MAX, NULL);
}

int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees,
		int threshold) {
	return rSPR_total_distance_unrooted(T1, gene_trees, threshold, NULL);
}

int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int threshold, vector<int> *original_scores) {
	//cout << "rSPR_total_distance_unrooted" << endl;
	int total = 0;
	MAIN_CALL = false;
	T1->preorder_number();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO) firstprivate(MAX_SPR) firstprivate(MIN_SPR)
	for(int i = 0; i < gene_trees.size(); i++) {
//		cout << "T1: " << T1->str_subtree() << endl;
//		cout << "T2: " << gene_trees[i]->str_subtree() << endl;
		Forest f1 = Forest(T1);
		//f1.print_components();
		Forest f2 = Forest(gene_trees[i]);
		//f2.print_components();
		if (!sync_twins(&f1, &f2))
			continue;
		if (f2.get_component(0)->get_children().size() > 2) {
			f2.get_component(0)->fixroot();
			f2.get_component(0)->set_depth(0);
			f2.get_component(0)->fix_depths();
			f2.get_component(0)->preorder_number();
		}
		//f1.print_components();
		//f2.print_components();
		int size = f2.get_component(0)->size();
		int best_distance = INT_MAX;
		int old_max = MAX_SPR;
		bool done = false;
		int NO_CLUSTER_ROUNDS=15;
//		cout << "boo" << endl;
		if (!UNROOTED_MIN_APPROX) {
//			for(int k = 0; !done; k++) {
			int best_min_spr = INT_MAX;
			vector<Node *> descendants = 
				f2.get_component(0)->find_descendants();
/*
			for(int j = 0; j < descendants.size(); j++) {
				f2.get_component(0)->reroot(descendants[j]);
				f2.get_component(0)->set_depth(0);
				f2.get_component(0)->fix_depths();
				f2.get_component(0)->preorder_number();
				Forest F1 = Forest(f1);
				Forest F2 = Forest(f2);
				int distance = rSPR_worse_3_approx(&F1, &F2)/3;
				if (distance < best_min_spr)
					best_min_spr = distance;
			}
			*/
		
			int min_spr = 0;
			if (best_min_spr < INT_MAX)
				min_spr = best_min_spr;
			for(int k = min_spr; k <= NO_CLUSTER_ROUNDS; k++) {
//			for(int k = min_spr; !done; k++) {
////				cout << k << endl;
				MIN_SPR=k;
				MAX_SPR=k;
//				Node *original_lc = f2.get_component(0)->lchild();
////					f2.print_components();
////					cout << endl;
//				vector<Node *> descendants = 
//				f2.get_component(0)->find_descendants();
				for(int j = 0; j < descendants.size(); j++) {
//					cout << "J=" << j << endl;
//					cout << i << "," << k << "," << j << endl;
//					cout << "rooting at: " << descendants[j]->str_subtree() << endl;
					//f2.get_component(0)->reroot(original_lc);
					f2.get_component(0)->reroot(descendants[j]);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
////					f2.print_components();
////					cout << endl;
	//			cout << T1->str_subtree() << endl;
	//			cout << gene_trees[i]->str_subtree() << endl;
					int distance;
					Forest *F1 = new Forest(f1);
					Forest *F2 = new Forest(f2);
					if (k <= NO_CLUSTER_ROUNDS)
						distance = rSPR_branch_and_bound_range(F1, F2, MIN_SPR, MAX_SPR);
//					else
//						break;
//						distance = rSPR_branch_and_bound_simple_clustering(F1->get_component(0), F2->get_component(0), VERBOSE, k, k);
					if (distance < 0)
						distance = k+1;
					delete F1;
					delete F2;
					if (distance <= k) {
						best_distance = distance;
						k=NO_CLUSTER_ROUNDS+1;
						done = true;
						break;
					}
				}
////				cout << endl;
				//f2.get_component(0)->reroot(original_lc);
//				cout << endl;
			}
			MAX_SPR=old_max;
			MIN_SPR=0;
			if (!done) {
				vector<Node *> descendants = 
					f2.get_component(0)->find_descendants();
				for(int j = 0; j < descendants.size(); j++) {
					f2.get_component(0)->reroot(descendants[j]);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
	//				cout << i << "," << j << endl;
	//				cout << T1->str_subtree() << endl;
	//				cout << gene_trees[i]->str_subtree() << endl;
					int distance = rSPR_branch_and_bound_simple_clustering(f1.get_component(0), f2.get_component(0), VERBOSE);
					if (distance <= best_distance) {
							best_distance = distance;
					}
				}
			}
	//		cout << "best_distance: " << best_distance << endl;
			if (best_distance == INT_MAX)
				best_distance = 0;
			total += best_distance;
			if (original_scores != NULL)
				(*original_scores)[i] = best_distance;
	//		cout << "total: " << total << endl;
		}
		else {
			int best_approx = INT_MAX;
			Node *best_rooting = f2.get_component(0)->lchild();
			int num_ties = 2;
			vector<Node *> descendants = 
				f2.get_component(0)->find_descendants();
			int NUM_ROOTINGS = 0;
//			int NUM_ROOTINGS = 6;
//			if (descendants.size() > 70)
//				NUM_ROOTINGS = descendants.size() / 10;
//			NUM_ROOTINGS = sqrt(descendants.size());
			if (NUM_ROOTINGS > 0 && descendants.size() < NUM_ROOTINGS + 1) {
				vector<Node *> rand_descendants =
					random_select(descendants, NUM_ROOTINGS);
				descendants = rand_descendants;
				descendants.push_back(f2.get_component(0)->lchild());
			}
			for(int j = 0; j < descendants.size(); j++) {
				f2.get_component(0)->reroot(descendants[j]);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
				//Forest F1 = Forest(f1);
				//Forest F2 = Forest(f2);
				int distance = rSPR_worse_3_approx_distance_only(&f1, &f2)/3;
				if (distance < best_approx) {
					best_approx = distance;
					best_rooting = descendants[j];
					num_ties = 2;
				}
				else if (distance == best_approx) {
					int r = rand();
					if (r < RAND_MAX/num_ties) {
						best_approx = distance;
						best_rooting = descendants[j];
					}
					num_ties++;
				}
			}
			f2.get_component(0)->reroot(best_rooting);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
			int k;
			if (best_approx > 20)
				k = rSPR_branch_and_bound_simple_clustering(f1.get_component(0), f2.get_component(0), VERBOSE);
			else
					k = rSPR_branch_and_bound_range(&f1, &f2, best_approx/3, best_approx);
			total += k;
		if (original_scores != NULL)
			(*original_scores)[i] = k;
		}
//		if (total > threshold)
//			break;
	}
	return total;
}

int rSPR_total_approx_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	int total = 0;
	MAIN_CALL = false;
	#pragma omp parallel for reduction(+: total)
	for(int i = 0; i < gene_trees.size(); i++) {
		Forest f1 = Forest(T1);
		Forest f2 = Forest(gene_trees[i]);
		if (!sync_twins(&f1, &f2))
			continue;
		if (f2.get_component(0)->get_children().size() > 2) {
			f2.get_component(0)->fixroot();
			f2.get_component(0)->set_depth(0);
			f2.get_component(0)->fix_depths();
			f2.get_component(0)->preorder_number();
		}
		int size = f2.get_component(0)->size();
		int best_distance = INT_MAX;
		vector<Node *> descendants = 
			f2.get_component(0)->find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			f2.get_component(0)->reroot(descendants[j]);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
			Forest F1 = Forest(f1);
			Forest F2 = Forest(f2);

			int distance = rSPR_worse_3_approx(&F1, &F2)/3;
			if (distance < best_distance)
				best_distance = distance;
		}
		if (best_distance == INT_MAX)
			best_distance = 0;
		total += best_distance;
	}
	return total;
}

int rSPR_total_approx_distance(Node *T1, vector<Node *> &gene_trees) {
	return rSPR_total_approx_distance(T1, gene_trees, INT_MAX);
}

int rSPR_total_approx_distance(Node *T1, vector<Node *> &gene_trees,
		int threshold) {
	int total = 0;
	MAIN_CALL = false;
	#pragma omp parallel for reduction(+ : total)
	for(int i = 0; i < gene_trees.size(); i++) {
		Forest F1 = Forest(T1);
		Forest F2 = Forest(gene_trees[i]);
//		cout << i << endl;
//		cout << T1->str_subtree() << endl;
//		cout << gene_trees[i]->str_subtree() << endl;
		total += rSPR_worse_3_approx(&F1, &F2)/3;
//		if (total > threshold)
//			break;
	}
	return total;
}


string itos(int i) {
	stringstream ss;
	string a;
	ss << i;
	a = ss.str();
	return a;
}

Node *find_subtree_of_approx_distance_hlpr(Node *n, Forest *F1, Forest *F2, int target_size) {
	Node *best_child_subtree = NULL;
	Node *native_child_subtree = NULL;
	int best_native_size = 0;
	int best_greedy_size = -1;
	int best_combined_size = INT_MIN;
	int largest_native_size = 0;
	list<Node *>::iterator c;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		Forest f1 = Forest(F1);
		Forest f2 = Forest(F2);
		Node *subtree = f1.find_by_prenum((*c)->get_preorder_number());
		f1.get_component(0)->disallow_siblings_subtree();
		subtree->allow_siblings_subtree();
//		if (subtree->lchild() != NULL)
//			subtree->lchild()->allow_siblings_subtree();
//		if (subtree->rchild() != NULL)
//			subtree->rchild()->allow_siblings_subtree();

		int cs_size = rSPR_worse_3_approx(subtree, &f1, &f2);
		if (native_child_subtree == NULL || cs_size > largest_native_size) {
			native_child_subtree = *c;
			largest_native_size = cs_size;
		}
		// h25: our induced-subproblem greedy chooses the split region.
		// Whidden still performs every cut in its native forest format.
		int greedy_size = h25_greedy_subtree_score(*c, F1, F2);
		int combined_size = cs_size +
			H25_GREEDY_SELECTOR_WEIGHT * max(0, greedy_size);
		bool better = best_child_subtree == NULL ||
			combined_size > best_combined_size ||
			(combined_size == best_combined_size &&
			 cs_size > best_native_size);
		if (better) {
			best_child_subtree = *c;
			best_native_size = cs_size;
			best_greedy_size = greedy_size;
			best_combined_size = combined_size;
		}
	}
	if (best_child_subtree != NULL) {
#ifdef H35_WHIDDEN_TRACE
		cerr << "h25 split_select pre="
			 << best_child_subtree->get_preorder_number()
			 << " leaves=" << best_child_subtree->find_leaves().size()
			 << " native=" << best_native_size
			 << " greedy=" << best_greedy_size
			 << " combined=" << best_combined_size
			 << " target=" << target_size
			 << " native_pre="
			 << native_child_subtree->get_preorder_number()
			 << " native_best=" << largest_native_size
			 << (best_child_subtree == native_child_subtree
				 ? " same" : " DIFFERENT")
			 << endl;
#endif
	}
	if (best_native_size <= 0)
		return n;
	else if (best_native_size < target_size)
		return best_child_subtree;
	else
		return find_subtree_of_approx_distance_hlpr(best_child_subtree,
				F1, F2, target_size);
}

Node *find_subtree_of_approx_distance(Node *n, Forest *F1, Forest *F2, int target_size) {
		Forest f1 = Forest(F1);
		Forest f2 = Forest(F2);
		Node *subtree = f1.find_by_prenum(n->get_preorder_number());
		f1.get_component(0)->disallow_siblings_subtree();
		subtree->allow_siblings_subtree();
//		if (subtree->lchild() != NULL)
//			subtree->lchild()->allow_siblings_subtree();
//		if (subtree->rchild() != NULL)
//			subtree->rchild()->allow_siblings_subtree();
		int size = rSPR_worse_3_approx(subtree, &f1, &f2);
		if (size > target_size)
			return find_subtree_of_approx_distance_hlpr(n, F1, F2, target_size);
		else 
			return n;
}

Node *find_best_root(Node *T1, Node *T2, double *best_root_b_acc) {
	Forest F1 = Forest(T1);
	Forest F2 = Forest(T2);
	Node *t1 = F1.get_component(0);
	Node *t2 = F2.get_component(0);
	//F1->preorder_number();
	int lchild_pre = t1->lchild()->get_preorder_number();
	int rchild_pre = t1->rchild()->get_preorder_number();
	if (!sync_twins(&F1, &F2)) {
		return NULL;
	}
//	if (t1->lchild()->get_preorder_number() != lchild_pre ||
//			t1->rchild()->get_preorder_number() != rchild_pre)
//		return NULL;
	if (F2.get_component(0)->get_children().size() > 2)
		F2.get_component(0)->fixroot();
	// TODO: maybe stop if F2 is too small?
	int pre_separator = t1->lchild()->get_preorder_number();
	int group_1_total;
	int group_2_total;
	if (t1->rchild()->get_preorder_number() > pre_separator) {
		pre_separator = t1->rchild()->get_preorder_number();
		group_1_total = t1->lchild()->find_leaves().size();
		group_2_total = t1->rchild()->find_leaves().size();
	}
	else {
		group_1_total = t1->rchild()->find_leaves().size();
		group_2_total = t1->lchild()->find_leaves().size();
	}
//	cout << "g1_total: " << group_1_total << endl;
//	cout << "g2_total: " << group_2_total << endl;
	Node *best_root = t2->lchild();
	find_best_root_hlpr(t2, pre_separator, group_1_total, group_2_total,
			&best_root, best_root_b_acc);
//	cout << "t1: " << t1->str_subtree() << endl;
//	cout << "t2: " << t2->str_subtree() << endl;
//	cout << "best_root: " << best_root->str_subtree() << endl;
	best_root = T2->find_by_prenum(best_root->get_preorder_number());
//	cout << "T1: " << T1->str_subtree() << endl;
//	cout << "best_root: " << best_root->str_subtree() << endl;
//	T2->reroot(best_root);
//	cout << "T2: " << T2->str_subtree() << endl;
	return best_root;
}

Node *find_best_root(Node *T1, Node *T2) {
	double best_root_b_acc = 0;
	Forest f1 = Forest(T1);
	Forest f2 = Forest(T2);
	sync_twins(&f1, &f2);
	Node *new_root =
		find_best_root(f1.get_component(0), f2.get_component(0), &best_root_b_acc);
	if (new_root != NULL)
		new_root = T2->find_by_prenum(new_root->get_preorder_number());
	return new_root;
}

double find_best_root_acc(Node *T1, Node *T2) {
	double best_root_b_acc = -1;
	Forest f1 = Forest(T1);
	Forest f2 = Forest(T2);
	sync_twins(&f1, &f2);
	find_best_root(f1.get_component(0), f2.get_component(0), &best_root_b_acc);
	return best_root_b_acc;
}

void find_best_root_hlpr(Node *T2, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc) {
	list<Node*>::iterator c;
	int group_1_descendants = 0;
	int group_2_descendants = 0;
	int num_ties = 2;
	for(c = T2->get_children().begin(); c != T2->get_children().end(); c++) {
		find_best_root_hlpr(*c, pre_separator, group_1_total,
				group_2_total, best_root, best_root_b_acc,
				&group_1_descendants, &group_2_descendants, &num_ties);
	}
}

/*	class child_ba_comp {
		private:
			int g_1_total;
			int g_2_total;
			bool d;
		public:
			child_ba_comp(int group_1_total, int group_2_total, bool direction) {
				g_1_total = group_1_total;
				g_2_total = group_2_total;
				d = direction;
			}
			bool operator()(const pair<int,int> x,const pair<int,int> y) {
				if (d) {
					return ((x.first / x.second * (x.first + x.second)) - 
							(y.first / y.second * (y.first + y.second)));
				}
				else {
					return ((x.second / x.first * (x.first + x.second)) - 
							(y.second / y.first * (y.first + y.second)));
				}
			}
	};
	*/

void find_best_root_hlpr(Node *n, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc,
		int *p_group_1_descendants, int *p_group_2_descendants, int *num_ties) {
	list<Node*>::iterator c;
	int group_1_descendants = 0;
	int group_2_descendants = 0;
//	vector<pair<int,int> > children_splits = vector<pair<int,int>>();
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
//		int g_1_desc = 0;
//		int g_2_desc = 0;
		find_best_root_hlpr(*c, pre_separator, group_1_total,
				group_2_total, best_root, best_root_b_acc,
				&group_1_descendants, &group_2_descendants, num_ties);
//				&g_1_desc, &g_2_desc, num_ties);
//		children_splits.push_back(make_pair(g_1_desc,g_2_desc));
//		group_1_descendants += g_1_desc;
//		group_2_descendants += g_2_desc;
	}
	if (n->is_leaf()) {
		int pre = n->get_twin()->get_preorder_number();
		if (pre < pre_separator)
			group_1_descendants++;
		else
			group_2_descendants++;
	}
//	else if (n->get_children.size() > 2) {
//		if (group_1_descendants / (double) group_1_total >= group_2_descendants / (double) group_2_total) {
//			sort(children_splits.begin(),children_splits.end(), g_1_comp(group_1_total, group_2_total, true));
//		}
//	}
	// balanced accuracy
	// don't bother averaging since we only directly compare them
	double tpos = group_1_descendants;
	double fpos = group_2_descendants;
	double fneg = (group_1_total - group_1_descendants);
	double tneg = (group_2_total - group_2_descendants);
	// balanced accuracy for this bipartition
	double b_acc =  tpos / (tpos + fneg)
			+ tneg / (tneg + fpos);
	// balanced accuracy for the opposite bipartition
	double b_acc_opp = fpos / (fpos + tneg)
			+ fneg / (fneg + tpos);
	// use the better bipartition
//	cout << "n: " << n->str_subtree() << endl;
//	cout << "b_acc: " << b_acc << endl;
//	cout << "b_acc_opp: " << b_acc_opp << endl;
//	cout << n->str_subtree() << endl;
//	cout << b_acc << "\t" << b_acc_opp << endl;
	if (b_acc_opp > b_acc)
		b_acc = b_acc_opp;
	if (b_acc > *best_root_b_acc) {
		*best_root = n;
		*best_root_b_acc = b_acc;
		*num_ties = 2;
	}
	else if(b_acc == *best_root_b_acc) {
		int r = rand();
		if (r < RAND_MAX/ *num_ties) {
			*best_root = n;
			*best_root_b_acc = b_acc;
		}
	}
	*p_group_1_descendants += group_1_descendants;
	*p_group_2_descendants += group_2_descendants;
}

Node *find_random_root(Node *T1, Node *T2) {
	vector<Node *> rroots = T2->find_descendants();
	int r = rand() % rroots.size();
	return rroots[r];
}
Node *find_best_root_rspr(Node *T1, Node *T2) {
	Node *t1 = new Node(*T1);
//	t1->preorder_number();
	Node *t2  = new Node(*T2);
	int new_prenum = T2->lchild()->get_preorder_number();
	vector<Node *> roots = t2->find_descendants();
	vector<int> root_prenums = vector<int>(roots.size());
	for(int i = 0; i < roots.size(); i++) {
		root_prenums[i] = roots[i]->get_preorder_number();
	}
	int best_distance = INT_MAX;
	int num_ties = 2;
//	cout << endl;
//	cout << "find_best_root_rspr" << endl;
//	cout << "T1: " << T1->str_subtree() << endl;
//	cout << "T2: " << T2->str_subtree() << endl;
//	cout << roots.size() << " Rootings" << endl;
	for(int i = 0; i < roots.size(); i++) {
		Node *root = roots[i];
		t2->reroot(root);
//		cout << "\t" << t2->str_subtree() << endl;
		t2->set_depth(0);
		t2->fix_depths();
		t2->preorder_number();
		int distance = rSPR_branch_and_bound_simple_clustering(t1, t2);
//		int distance = rf_distance(t1, t2);
		if (distance < best_distance) {
			best_distance = distance;
			new_prenum = root_prenums[i];
			num_ties = 2;
		}
		else if (distance == best_distance) {
			int r = rand();
			if (r < RAND_MAX / num_ties) {
				best_distance = distance;
				new_prenum = root_prenums[i];
			}
			num_ties++;
		}
	}
	Node *new_root = T2->find_by_prenum(new_prenum);
	t1->delete_tree();
	t2->delete_tree();
	return new_root;
}

// assume already sync_twins and preorder numbered
bool contains_bipartition(Node *n, int pre_start, int pre_end,
		int group_1_total, int group_2_total, int *p_group_1_descendants,
		int *p_group_2_descendants) {
	list<Node*>::iterator c;
	int group_1_descendants = 0;
	int group_2_descendants = 0;
	bool found = false;
	bool proper_split = true;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		int c_group_1_descendants = 0;
		int c_group_2_descendants = 0;
		found = contains_bipartition(*c, pre_start, pre_end, group_1_total,
				group_2_total, &c_group_1_descendants, &c_group_2_descendants);
		if (found)
			return true;
		group_1_descendants += c_group_1_descendants;
		group_2_descendants += c_group_2_descendants;
		if (c_group_1_descendants > 0 && c_group_2_descendants > 0)
			proper_split = false;
	}

	if (n->is_leaf()) {
		int pre = n->get_twin()->get_preorder_number();
		if (pre >= pre_start && pre <= pre_end)
			group_1_descendants++;
		else
			group_2_descendants++;
	}
	else if (proper_split) {
		if (group_1_descendants == group_1_total
				|| group_2_descendants == group_2_total)
			return true;
	}
	if (p_group_1_descendants != NULL)
		*p_group_1_descendants += group_1_descendants;
	if (p_group_2_descendants != NULL)
		*p_group_2_descendants += group_2_descendants;
	return false;
}

// root the tree based on an outgroup
// returns false if the outgroup is not found or not a clade
bool outgroup_root(Node *T, set<string, StringCompare> outgroup) {
	vector<int> num_in = vector<int>();
	vector<int> num_out = vector<int>();
	count_in_out(T, num_in, num_out, outgroup);
	list<Node *>::iterator c;
	int pre = T->get_preorder_number();
	bool clean_split = true;
	if (num_out[pre] == 0)
		return false;
	for(c = T->get_children().begin(); c != T->get_children().end(); c++) {
		int c_pre = (*c)->get_preorder_number();
		if (num_out[pre] == num_out[c_pre]) {
			// split the out_group
			return outgroup_root(*c, num_in, num_out);
		}
		else if (num_in[pre] == num_in[c_pre]) {
			// split the in_group
			return outgroup_root(*c, num_out, num_in);
		}
		else if (num_out[c_pre] > 0 && num_in[c_pre] > 0)
			clean_split = false;
	}
	if (clean_split)
		return outgroup_reroot(T, num_in, num_out);
	else
		return false;
}

bool outgroup_root(Node *n, vector<int> &num_in, vector<int> &num_out) {
	list<Node *>::iterator c;
	int pre = n->get_preorder_number();
	bool clean_split = true;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		int c_pre = (*c)->get_preorder_number();
		if (num_out[pre] == num_out[c_pre]) {
			// split the out_group
			return outgroup_root(*c, num_in, num_out);
		}
		else if (num_out[c_pre] > 0 && num_in[c_pre] > 0)
			clean_split = false;
	}
	if (clean_split)
		return outgroup_reroot(n, num_in, num_out);
	else
		return false;
}

bool outgroup_reroot(Node *n, vector<int> &num_in, vector<int> &num_out) {
	Node *T = n->find_root();
	if (num_in[n->get_preorder_number()] == 0) {
		if (!n->is_leaf()) {
			T->reroot(n);
			T->set_depth(0);
			T->fix_depths();
			T->preorder_number();
		}
		return true;
	}
	Node *new_split = new Node("");
	n->add_child(new_split);
	list<Node *>::iterator c = n->get_children().begin();
	while(c != n->get_children().end()) {
		Node *child = *c;
		c++;
		int c_pre = child->get_preorder_number();
		if (c_pre >= 0 && num_in[c_pre] == 0) {
			//child->cut_parent();
			new_split->add_child(child);
		}
	}
	T->reroot(new_split);
	T->set_depth(0);
	T->fix_depths();
	T->preorder_number();
	return true;
}

void count_in_out(Node *n, vector<int> &num_in, vector<int> &num_out,
		set<string, StringCompare> &outgroup) {
	list<Node *>::iterator c;
	int pre = n->get_preorder_number();
	if (num_in.size() <= pre)
		num_in.resize(pre + 1, 0);
	if (num_out.size() <= pre)
		num_out.resize(pre + 1, 0);
	if (n->is_leaf()) {
		if (outgroup.find(n->get_name()) != outgroup.end()) {
				num_out[pre] = 1;
				num_in[pre] = 0;
		}
		else {
				num_out[pre] = 0;
				num_in[pre] = 1;
		}
	}
	else {
		for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
			count_in_out(*c, num_in, num_out, outgroup);
			num_in[pre] += num_in[(*c)->get_preorder_number()];
			num_out[pre] += num_out[(*c)->get_preorder_number()];
		}
	}
}

void modify_bipartition_support(Node *T1, Node *T2, enum RELAXATION relaxed) {
	Forest F1 = Forest(T1);
	Forest F2 = Forest(T2);
	if (!sync_twins(&F1, &F2)) {
		return;
	}
	Node *n = F1.get_component(0);
	vector<int> *F1_descendant_counts = n->find_leaf_counts();
	modify_bipartition_support(n, &F1, &F2, T1, T2, F1_descendant_counts,
			relaxed);
	delete F1_descendant_counts;
}

void modify_bipartition_support(Node *n, Forest *F1, Forest *F2,
		Node *T1, Node *T2, vector<int> *F1_descendant_counts, enum RELAXATION relaxed) {
	if (n->is_leaf())
		return;
	list<Node *>::iterator c;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		modify_bipartition_support(*c, F1, F2, T1, T2, F1_descendant_counts,
				relaxed);
	}
	Node *t = T1->find_by_prenum(n->get_preorder_number());
	if (t->parent() == NULL)
		return;
	int pre_start = n->get_edge_pre_start();
	int pre_end = n->get_edge_pre_end();
	int group_1_total = (*F1_descendant_counts)[n->get_preorder_number()];
	int group_2_total = (*F1_descendant_counts)[F1->get_component(0)->get_preorder_number()] - group_1_total;
	if (group_2_total >= 2) {
		if (contains_bipartition(F2->get_component(0), pre_start, pre_end,
				group_1_total, group_2_total, NULL, NULL)) {
			t->a_inc_support();
			t->a_inc_support_normalization();
			// relaxed
			if (false && relaxed == ALL_RELAXED) {
				int stop_pre = 0;
				if (n->parent() != NULL) {
					stop_pre = n->parent()->get_preorder_number();
				while ((t = t->parent()) != NULL
						&& t->get_preorder_number() != stop_pre)
					t->a_inc_support();
					t->a_inc_support_normalization();
				}
			}
		}
		else {
//			t->a_dec_support();
			t->a_inc_support_normalization();
			// relaxed
			if (relaxed == ALL_RELAXED || relaxed == NEGATIVE_RELAXED) {
				int stop_pre = 0;
				if (n->parent() != NULL) {
					stop_pre = n->parent()->get_preorder_number();
				while ((t = t->parent()) != NULL
						&& t->get_preorder_number() != stop_pre)
//					t->a_dec_support();
					t->a_inc_support_normalization();
				}
			}
		}
	}
}

int rf_distance(Node *T1, Node *T2) {
	Forest F1 = Forest(T1);
	Forest F2 = Forest(T2);
	if (!sync_twins(&F1, &F2))
		return 0;
	if (F1.get_component(0)->is_leaf())
		return 0;
	sync_interior_twins(&F1, &F2);
	int rf_d = 0;
	rf_d += count_differing_bipartitions(F1.get_component(0));
	rf_d += count_differing_bipartitions(F2.get_component(0));
	return rf_d;
}

int count_differing_bipartitions(Node *n) {
	//cout << "Start: " << n->str_subtree() << endl;
	int count = 0;
	list<Node *>::iterator c;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		count += count_differing_bipartitions(*c);
	}
	if (n->get_twin() == NULL ||
//			n->get_depth() > n->get_twin()->get_twin()->get_depth())
			n != n->get_twin()->get_twin()) {
		count++;
			}
	return count;
}

bool is_nonbranching(Forest *T1, Forest *T2, Node *T1_a, Node *T1_c, Node *T2_a, Node *T2_c) {
	if ((T2_a->get_depth() < T2_c->get_depth()
			&& T2_c->parent() != NULL)
			|| T2_a->parent() == NULL) {
		swap(&T1_a, &T1_c);
		swap(&T2_a, &T2_c);
	}
	else if (T2_a->get_depth() == T2_c->get_depth()) {
		if (T2_a->parent() && T2_c->parent() &&
				(T2_a->parent()->get_depth() <
				T2_c->parent()->get_depth()
				//|| (T2_c->parent()->parent()
				//&& T2_c->parent()->parent() == T2_a->parent())
				)) {
		swap(&T1_a, &T1_c);
		swap(&T2_a, &T2_c);
		}
	}
	int num_protected = T2_a->is_protected() + T2_c->is_protected();
	if (T2_a->parent()->get_children().size() == 2)
		num_protected += T2_a->get_sibling()->is_protected();
	if (num_protected >= 2)
		return true;
	if (CUT_ONE_B) {
		if (T2_a->parent()->parent() == T2_c->parent()
			&& T2_c->parent() != NULL
			&& T2_a->parent()->get_children().size() <= 2)
			return true;
	}
	if (CUT_TWO_B && T1_a->parent()->parent() != NULL) {
		Node *T1_s = T1_a->parent()->get_sibling();
		if (T1_s->is_leaf()) {
			Node *T2_l = T2_a->parent()->parent();
			if (T2_l != NULL && T2_l->get_children().size() <= 2) {
				if (T2_c->parent() != NULL && T2_c->parent()->parent() == T2_l
						&& ((T2_a->parent()->get_children().size() <= 2
						&& T2_c->parent()->get_children().size() <= 2)
						|| T1_s->get_twin()->is_protected())){
					if (T2_l->get_sibling() == T1_s->get_twin()) {
						return true;
					}
					else if (CUT_TWO_B_ROOT && T2_l->parent() == NULL &&
							(T2->contains_rho() ||
							 T2->get_component(0) != T2_l)) {
						return true;
					}
				}
				else if ((T2_l = T2_l->parent()) != NULL
						&& T2_c->parent() == T2_l
						&& ((T2_a->parent()->get_children().size() <= 2
						&& T2_a->parent()->parent()->get_children().size() <= 2
						&& T2_l->get_children().size() <= 2)
						|| T1_s->get_twin()->is_protected())){
					if (T2_l->get_sibling() == T1_s->get_twin()) {
						return true;
					}
					else if (CUT_TWO_B_ROOT && T2_l->parent() == NULL &&
							(T2->contains_rho() ||
							 T2->get_component(0) != T2_l)) {
						return true;
					}
				}
			}
		}
	}
	if (REVERSE_CUT_ONE_B && T1_a->parent()->parent() != NULL) {
		Node *T1_s = T1_a->parent()->get_sibling();
		Node *T2_s = T1_s->get_twin();
		if (T1_s->is_leaf()) {
			if (T2_s->parent() == T2_a->parent()) {
				return true;
			}
			else if (T2_s->parent() == T2_c->parent()
					&& T2_c->parent()->get_children().size() <= 2) {
				return true;
			}
//			else if (REVERSE_CUT_ONE_B_3
//							&& T2_s->is_protected()
//							&& T2_s->parent() != NULL
//							&& T2_s->parent()->parent() == T2_a->parent()
//							&& T2_s->parent()->get_children().size() <= 2) {
//				return true;
//			}
		}
		else if (REVERSE_CUT_ONE_B_2 && T2_c->parent() != NULL
				&& chain_match(T1_s, T2_c->get_sibling(), T2_a))
			return true;
	}
	return false;
}

void strip_whitespace(string &str) {
	std::string::iterator end_pos = std::remove_if(str.begin(), str.end(), ::isspace);
	str.erase(end_pos, str.end());
}

void strip_trailing_whitespace(string &str) {
	size_t start_pos = str.find_first_not_of(whitespaces);
	str.erase(0, start_pos);
	size_t end_pos = str.find_last_not_of(whitespaces)+1;
	str.erase(end_pos, str.size()-end_pos);
}

//randomizes T2 with count number of sprs. If T1 equals T2 at the beginning,
//then the spr distance between T1 and T2 is equal to (or possibly less than) count
//Assumes they are already synced, assumes there are valid sprs to be made
void randomize_tree_with_spr(Forest* T1, Forest* T2, int count) {
  for (int spr = 0; spr < count; spr++) {
    vector<Node*> all_nodes = T2->get_component(0)->find_nodes_in_subtree();
    Node* source = all_nodes[rand() % all_nodes.size()];
    Node* target = all_nodes[rand() % all_nodes.size()];
    bool target_in_subtree = false;
    Node* first_leaf = target->find_leaves()[0];
    vector<Node*> source_leaves = source->find_leaves();
    for (int i = 0; i < source_leaves.size(); i++) {
      if (source_leaves[i] == first_leaf) {
	target_in_subtree = true;
	break;
      }
    }

    //Get random node
    //if target is in source's subtree repick
    //spr
    while (source == target ||
	   source->is_sibling_of(target) ||
	   source->parent() == target || 
	   target_in_subtree) {	   
      source = all_nodes[rand() % all_nodes.size()];
      target = all_nodes[rand() % all_nodes.size()];
      //cout << "Trying: " << source->str_subtree() << " and "<<  target->str_subtree() << endl;
      target_in_subtree = false;
      Node* first_leaf = target->find_leaves()[0];
      
      int child_count = source->get_children().size();
      if (child_count > 2) {
	int rand_count = rand() % (child_count + 1);
	if (rand_count == child_count || rand_count == 0){
	  //cout << "moving whole tree" << endl;
	}
	else if (rand_count == 1) {
	  source = source->get_children().front();
	  //cout << "moving first child" << endl;
	}
	else {
	  list<Node*> to_expand = list<Node*>();
	  list<Node*>::iterator c = source->get_children().begin();
	  for (int i = 0; i < rand_count; i++) {
	    to_expand.push_back(*c);
	    c++;
	  }	
	  source = source->expand_children_out(to_expand);
	  //cout << "Moving part " << rand_count<< endl;
	}
      }
      
      vector<Node*> source_leaves = source->find_leaves();
      for (int i = 0; i < source_leaves.size(); i++) {
	if (source_leaves[i] == first_leaf) {
	  //cout << "target in subtree" << endl;
	  target_in_subtree = true;
	  break;
	}
      }
    } 

    //cout << "Moving : " << source->str_subtree() << " to " << target->str_subtree() << endl;
    
    Node* parent = source->parent();      
    if (parent != NULL) {
      source->cut_parent();
      if (parent->get_children().size() == 1) {
	parent->contract(true);
      }
    }
    else {
      continue;
    }
    //regraft
    if (target->parent() != NULL) {
      target->parent()->add_child(source);
    }
    else {
      Node* new_parent = target;//new Node();
      Node* replace = new Node(*target);
      new_parent->get_children().clear();
      new_parent->add_child(replace);
      new_parent->add_child(source);
    }

    //T1->print_components();
    //T2->print_components();
  }
}

// ===== END heuristic/clean_refactor/h36_whidden/rspr.h =====

using namespace std;
using Hash = uint64_t;

// Search policy constants.
static constexpr int SAFE_RECURSION_DEPTH = 1024;
#ifndef LOCAL_COMPONENT_WINDOW_VALUE
#define LOCAL_COMPONENT_WINDOW_VALUE 24
#endif
static constexpr int LOCAL_COMPONENT_WINDOW =
    LOCAL_COMPONENT_WINDOW_VALUE;
#ifndef LOOKAHEAD_NEIGHBOR_LIMIT_VALUE
#define LOOKAHEAD_NEIGHBOR_LIMIT_VALUE 16
#endif
static constexpr int LOOKAHEAD_NEIGHBOR_LIMIT =
    LOOKAHEAD_NEIGHBOR_LIMIT_VALUE;
#ifndef LOOKAHEAD_TOP_CANDIDATES_VALUE
#define LOOKAHEAD_TOP_CANDIDATES_VALUE 8
#endif
static constexpr int LOOKAHEAD_TOP_CANDIDATES =
    LOOKAHEAD_TOP_CANDIDATES_VALUE;
static constexpr int MAX_SPLIT_BLOCKS = 8;
static constexpr int CONFLICT_COMPONENT_WINDOW = 32;
static constexpr int MAX_CONFLICT_PAIRS_PER_TREE = 512;
static constexpr int MAX_CONFLICT_POOL_SIZE = 256;
static constexpr int LOCAL_REGION_COMPONENT_LIMIT = 48;
#ifndef NEIGHBOR_WALK_LIMIT_VALUE
#define NEIGHBOR_WALK_LIMIT_VALUE 20
#endif
static constexpr int NEIGHBOR_WALK_LIMIT =
    NEIGHBOR_WALK_LIMIT_VALUE;
static constexpr int SPLIT_APPROX_TARGET_COMPONENTS = 24;
static constexpr int SPLIT_APPROX_COMPONENT_LIMIT = 49;
static constexpr int SPLIT_APPROX_LEAF_LIMIT = 384;
#ifndef SPLIT_APPROX_START_ATTEMPT_VALUE
#define SPLIT_APPROX_START_ATTEMPT_VALUE 1000000000
#endif
static constexpr int SPLIT_APPROX_START_ATTEMPT =
    SPLIT_APPROX_START_ATTEMPT_VALUE;
#ifndef SPLIT_APPROX_INTERVAL_VALUE
#define SPLIT_APPROX_INTERVAL_VALUE 256
#endif
static constexpr int SPLIT_APPROX_INTERVAL =
    SPLIT_APPROX_INTERVAL_VALUE;
#ifndef CONFLICT_RESTART_INTERVAL_VALUE
#define CONFLICT_RESTART_INTERVAL_VALUE 2
#endif
static constexpr int CONFLICT_RESTART_INTERVAL =
    CONFLICT_RESTART_INTERVAL_VALUE;
#ifndef LOOKAHEAD_RESTART_INTERVAL_VALUE
#define LOOKAHEAD_RESTART_INTERVAL_VALUE 16
#endif
static constexpr int LOOKAHEAD_RESTART_INTERVAL =
    LOOKAHEAD_RESTART_INTERVAL_VALUE;
static constexpr int JITTER_CYCLE_LENGTH = 24;
#ifndef WHIDDEN_RESTART_BUDGET_VALUE
#define WHIDDEN_RESTART_BUDGET_VALUE 6144
#endif
static constexpr int WHIDDEN_RESTART_BUDGET =
    WHIDDEN_RESTART_BUDGET_VALUE;
#ifndef WHIDDEN_RESTART_SEED_OFFSET_VALUE
#define WHIDDEN_RESTART_SEED_OFFSET_VALUE 0
#endif
static constexpr int WHIDDEN_RESTART_SEED_OFFSET =
    WHIDDEN_RESTART_SEED_OFFSET_VALUE;
#ifndef WHIDDEN_EXTRA_SPLITS_VALUE
#define WHIDDEN_EXTRA_SPLITS_VALUE 0
#endif
static constexpr int WHIDDEN_EXTRA_SPLITS =
    WHIDDEN_EXTRA_SPLITS_VALUE;
#ifndef WHIDDEN_CLUSTER_TUNE_VALUE
// Preserve Whidden's original clustering setup for this experiment.
#define WHIDDEN_CLUSTER_TUNE_VALUE 30
#endif
static constexpr int WHIDDEN_CLUSTER_TUNE =
    WHIDDEN_CLUSTER_TUNE_VALUE;
#ifndef WHIDDEN_SPLIT_THRESHOLD_VALUE
#define WHIDDEN_SPLIT_THRESHOLD_VALUE 25
#endif
static constexpr int WHIDDEN_SPLIT_THRESHOLD =
    WHIDDEN_SPLIT_THRESHOLD_VALUE;
#ifndef WHIDDEN_RANDOM_SEED_VALUE
#define WHIDDEN_RANDOM_SEED_VALUE 7
#endif
static constexpr int WHIDDEN_RANDOM_SEED =
    WHIDDEN_RANDOM_SEED_VALUE;
#ifndef SECOND_WHIDDEN_CLUSTER_TUNE_VALUE
#define SECOND_WHIDDEN_CLUSTER_TUNE_VALUE -1
#endif
static constexpr int SECOND_WHIDDEN_CLUSTER_TUNE =
    SECOND_WHIDDEN_CLUSTER_TUNE_VALUE;
#ifndef SECOND_WHIDDEN_SPLIT_THRESHOLD_VALUE
#define SECOND_WHIDDEN_SPLIT_THRESHOLD_VALUE 45
#endif
static constexpr int SECOND_WHIDDEN_SPLIT_THRESHOLD =
    SECOND_WHIDDEN_SPLIT_THRESHOLD_VALUE;
#ifndef SECOND_WHIDDEN_RANDOM_SEED_VALUE
#define SECOND_WHIDDEN_RANDOM_SEED_VALUE 7
#endif
static constexpr int SECOND_WHIDDEN_RANDOM_SEED =
    SECOND_WHIDDEN_RANDOM_SEED_VALUE;
#ifndef DEDICATED_WHIDDEN_SCORE_THRESHOLD_VALUE
#define DEDICATED_WHIDDEN_SCORE_THRESHOLD_VALUE 3000
#endif
static constexpr int DEDICATED_WHIDDEN_SCORE_THRESHOLD =
    DEDICATED_WHIDDEN_SCORE_THRESHOLD_VALUE;
#ifndef DEDICATED_WHIDDEN_POLL_LIMIT_VALUE
#define DEDICATED_WHIDDEN_POLL_LIMIT_VALUE 4200
#endif
static constexpr int DEDICATED_WHIDDEN_POLL_LIMIT =
    DEDICATED_WHIDDEN_POLL_LIMIT_VALUE;
#ifndef DEDICATED_WHIDDEN_CLUSTER_TUNE_VALUE
#define DEDICATED_WHIDDEN_CLUSTER_TUNE_VALUE 30
#endif
static constexpr int DEDICATED_WHIDDEN_CLUSTER_TUNE =
    DEDICATED_WHIDDEN_CLUSTER_TUNE_VALUE;
#ifndef DEDICATED_WHIDDEN_SPLIT_THRESHOLD_VALUE
#define DEDICATED_WHIDDEN_SPLIT_THRESHOLD_VALUE 25
#endif
static constexpr int DEDICATED_WHIDDEN_SPLIT_THRESHOLD =
    DEDICATED_WHIDDEN_SPLIT_THRESHOLD_VALUE;
#ifndef KEEP_EARLY_WHIDDEN_VALUE
#define KEEP_EARLY_WHIDDEN_VALUE 0
#endif
static constexpr bool KEEP_EARLY_WHIDDEN =
    KEEP_EARLY_WHIDDEN_VALUE != 0;
#ifndef RANDOM_CANDIDATE_MULTIPLIER_VALUE
#define RANDOM_CANDIDATE_MULTIPLIER_VALUE 2
#endif
static constexpr int RANDOM_CANDIDATE_MULTIPLIER =
    RANDOM_CANDIDATE_MULTIPLIER_VALUE;
#ifndef ENABLE_COMMON_CONTRACTION_VALUE
#define ENABLE_COMMON_CONTRACTION_VALUE 1
#endif
static constexpr bool ENABLE_COMMON_CONTRACTION =
    ENABLE_COMMON_CONTRACTION_VALUE != 0;
static constexpr int MIN_MERGE_CHECKS = 10000;
static constexpr int PERIODIC_TASK_INTERVAL = 8;
#ifndef MERGE_CHECK_MULTIPLIER_VALUE
#define MERGE_CHECK_MULTIPLIER_VALUE 80
#endif
static constexpr int MERGE_CHECK_MULTIPLIER =
    MERGE_CHECK_MULTIPLIER_VALUE;
#ifndef MAX_LOOKAHEAD_CHOICES_VALUE
#define MAX_LOOKAHEAD_CHOICES_VALUE 4
#endif
static constexpr int MAX_LOOKAHEAD_CHOICES = MAX_LOOKAHEAD_CHOICES_VALUE;
#ifndef BIG_SCORE_MODE
#define BIG_SCORE_MODE 0
#endif
#ifndef PRIMARY_SINGLETON_SCORE_MODE_VALUE
#define PRIMARY_SINGLETON_SCORE_MODE_VALUE 0
#endif
static constexpr int PRIMARY_SINGLETON_SCORE_MODE =
    PRIMARY_SINGLETON_SCORE_MODE_VALUE;
#ifndef PRIMARY_RESTART_SCORE_MODE_VALUE
#define PRIMARY_RESTART_SCORE_MODE_VALUE 0
#endif
static constexpr int PRIMARY_RESTART_SCORE_MODE =
    PRIMARY_RESTART_SCORE_MODE_VALUE;
#ifndef PRIMARY_EXTRA_SPLITS_VALUE
#define PRIMARY_EXTRA_SPLITS_VALUE 0
#endif
static constexpr int PRIMARY_EXTRA_SPLITS =
    PRIMARY_EXTRA_SPLITS_VALUE;
#ifndef RESTART_ACCEPT_WORSE_MARGIN_VALUE
#define RESTART_ACCEPT_WORSE_MARGIN_VALUE 0
#endif
static constexpr int RESTART_ACCEPT_WORSE_MARGIN =
    RESTART_ACCEPT_WORSE_MARGIN_VALUE;
#ifndef STRUCTURED_RESTART_SPLIT_VALUE
#define STRUCTURED_RESTART_SPLIT_VALUE 0
#endif
static constexpr bool STRUCTURED_RESTART_SPLIT =
    STRUCTURED_RESTART_SPLIT_VALUE != 0;
#ifndef SIMPLE_SINGLETON_RESTART_BUDGET_VALUE
#define SIMPLE_SINGLETON_RESTART_BUDGET_VALUE 3072
#endif
static constexpr int SIMPLE_SINGLETON_RESTART_BUDGET =
    SIMPLE_SINGLETON_RESTART_BUDGET_VALUE;
#ifndef BEST_RESTART_SEED_OFFSET_VALUE
#define BEST_RESTART_SEED_OFFSET_VALUE 193
#endif
static constexpr int BEST_RESTART_SEED_OFFSET =
    BEST_RESTART_SEED_OFFSET_VALUE;
#ifndef PRIMARY_RESTART_SEED_OFFSET_VALUE
#define PRIMARY_RESTART_SEED_OFFSET_VALUE 0
#endif
static constexpr int PRIMARY_RESTART_SEED_OFFSET =
    PRIMARY_RESTART_SEED_OFFSET_VALUE;

// Signal handling.
static volatile sig_atomic_t g_stop_requested = 0;
static const function<void()> *g_periodic_task = nullptr;
static bool g_periodic_task_running = false;
static int g_periodic_task_tick = 0;
static constexpr size_t SIGNAL_BEST_OUTPUT_CAP = 8u * 1024u * 1024u;
static char g_signal_best_output[SIGNAL_BEST_OUTPUT_CAP];
static volatile sig_atomic_t g_signal_best_output_len = 0;

static void update_signal_best_output(const string &text) {
    size_t n = min(text.size(), SIGNAL_BEST_OUTPUT_CAP - 1);
    g_signal_best_output_len = 0;
    memcpy(g_signal_best_output, text.data(), n);
    if (n == 0 || g_signal_best_output[n - 1] != '\n') {
        g_signal_best_output[n++] = '\n';
    }
    g_signal_best_output_len = (sig_atomic_t)n;
}

static void handle_signal(int) {
    g_stop_requested = 1;
    sig_atomic_t n = g_signal_best_output_len;
    if (n > 0) {
        const char *data = g_signal_best_output;
        sig_atomic_t written = 0;
        while (written < n) {
            ssize_t rc = write(STDOUT_FILENO, data + written,
                               (size_t)(n - written));
            if (rc <= 0) break;
            written += (sig_atomic_t)rc;
        }
        _exit(0);
    }
}

static void maybe_run_periodic_task(bool force = false) {
    if (g_periodic_task == nullptr || g_periodic_task_running) return;
    if (!force && (++g_periodic_task_tick % PERIODIC_TASK_INTERVAL) != 0)
        return;
    g_periodic_task_running = true;
    (*g_periodic_task)();
    g_periodic_task_running = false;
}

#ifdef TRACE_PORTFOLIO
// Timestamp every portfolio event so STRIDE logs show when each route returns.
static const chrono::steady_clock::time_point g_trace_start =
    chrono::steady_clock::now();
static double trace_elapsed_seconds() {
    return chrono::duration<double>(
        chrono::steady_clock::now() - g_trace_start).count();
}
#define TRACE_PORTFOLIO_MSG(msg)                                             \
    do {                                                                     \
        cerr << "[t=" << trace_elapsed_seconds() << "s] " << msg << '\n';    \
    } while (0)
#else
#define TRACE_PORTFOLIO_MSG(msg) do {} while (0)
#endif

// Tree parsing and traversal.
struct TwoTreeInstance {
    int leaf_count = 0;
    array<string, 2> trees;
};

static bool read_two_tree_instance(istream &input, TwoTreeInstance &instance) {
    int declared_tree_count = 0;
    int tree_index = 0;
    string line;

    while (getline(input, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (line.size() >= 2 && line[1] == 'p') {
                string tag;
                stringstream header(line);
                header >> tag >> declared_tree_count >> instance.leaf_count;
            }
            continue;
        }
        if (tree_index < 2) instance.trees[tree_index++] = line;
    }

    return declared_tree_count == 2 && instance.leaf_count > 0 &&
           tree_index == 2;
}

struct TreeNode {
    int label = -1;
    int parent = -1;
    int depth = 0;
    vector<int> child;
};

struct Tree {
    vector<TreeNode> nodes;
    int root = -1;
    int n = 0;
    int LOG = 1;
    int max_depth = 0;
    vector<int> leaf_node;
    vector<vector<int>> up;
    vector<int> subtree_leaves;
    vector<int> congestion_prefix;

    int add_node(int label) {
        TreeNode x;
        x.label = label;
        nodes.push_back(x);
        return (int)nodes.size() - 1;
    }
};

static Tree parse_newick(const string &s, int n) {
    Tree T;
    T.n = n;
    T.leaf_node.assign(n + 1, -1);

    vector<int> parents;
    bool expect_subtree = true;
    for (int i = 0; i < (int)s.size();) {
        unsigned char ch = (unsigned char)s[i];
        if (isspace(ch)) {
            i++;
            continue;
        }
        if (s[i] == '(') {
            int v = T.add_node(-1);
            if (parents.empty()) {
                if (T.root < 0) T.root = v;
            } else {
                T.nodes[v].parent = parents.back();
                T.nodes[parents.back()].child.push_back(v);
            }
            parents.push_back(v);
            expect_subtree = true;
            i++;
            continue;
        }
        if (s[i] == ',') {
            expect_subtree = true;
            i++;
            continue;
        }
        if (s[i] == ')') {
            if (!parents.empty()) parents.pop_back();
            expect_subtree = false;
            i++;
            continue;
        }
        if (s[i] == ';') break;
        if (expect_subtree && isdigit(ch)) {
            int label = 0;
            while (i < (int)s.size() && isdigit((unsigned char)s[i])) {
                label = label * 10 + (s[i] - '0');
                i++;
            }
            int v = T.add_node(label);
            if (parents.empty()) {
                if (T.root < 0) T.root = v;
            } else {
                T.nodes[v].parent = parents.back();
                T.nodes[parents.back()].child.push_back(v);
            }
            expect_subtree = false;
            continue;
        }
        i++;
    }
    return T;
}

struct TraversalFrame {
    int node;
    size_t next_child;
};

static string restricted_newick_recursive(const Tree &T, int v,
                                          const vector<int> &mark, int token) {
    const TreeNode &node = T.nodes[v];
    if (node.label >= 1) {
        if (node.label < (int)mark.size() && mark[node.label] == token)
            return to_string(node.label);
        return "";
    }
    vector<string> children;
    for (int c : node.child) {
        string s = restricted_newick_recursive(T, c, mark, token);
        if (!s.empty()) children.push_back(std::move(s));
    }
    if (children.empty()) return "";
    if (children.size() == 1) return std::move(children[0]);
    string out = "(";
    for (size_t i = 0; i < children.size(); i++) {
        if (i) out.push_back(',');
        out += children[i];
    }
    out.push_back(')');
    return out;
}

static string restricted_newick_iterative(const Tree &T, int v,
                                          const vector<int> &mark, int token) {
    static vector<TraversalFrame> traversal;
    static vector<unsigned char> nonempty;
    static vector<int> kept_children;
    if (nonempty.size() < T.nodes.size()) nonempty.resize(T.nodes.size());
    if (kept_children.size() < T.nodes.size()) kept_children.resize(T.nodes.size());
    traversal.clear();
    if (v < 0 || v >= (int)T.nodes.size()) return "";
    traversal.push_back({v, 0});

    while (!traversal.empty()) {
        int u = traversal.back().node;
        const TreeNode &node = T.nodes[u];
        if (node.label >= 1) {
            nonempty[u] = node.label < (int)mark.size() &&
                          mark[node.label] == token;
            kept_children[u] = 0;
            traversal.pop_back();
            continue;
        }
        if (traversal.back().next_child < node.child.size()) {
            int c = node.child[traversal.back().next_child++];
            traversal.push_back({c, 0});
            continue;
        }

        int count = 0;
        for (int c : node.child)
            if (nonempty[c]) count++;
        kept_children[u] = count;
        nonempty[u] = count > 0;
        traversal.pop_back();
    }

    if (!nonempty[v]) return "";

    struct OutputTask {
        int node;
        char symbol;
    };
    static vector<OutputTask> tasks;
    tasks.clear();
    tasks.push_back({v, 0});
    string out;

    while (!tasks.empty()) {
        OutputTask task = tasks.back();
        tasks.pop_back();
        if (task.symbol != 0) {
            out.push_back(task.symbol);
            continue;
        }

        int u = task.node;
        while (T.nodes[u].label < 1 && kept_children[u] == 1) {
            for (int c : T.nodes[u].child) {
                if (nonempty[c]) {
                    u = c;
                    break;
                }
            }
        }

        const TreeNode &node = T.nodes[u];
        if (node.label >= 1) {
            out += to_string(node.label);
            continue;
        }

        out.push_back('(');
        tasks.push_back({-1, ')'});
        vector<int> children;
        children.reserve((size_t)kept_children[u]);
        for (int c : node.child)
            if (nonempty[c]) children.push_back(c);
        for (int i = (int)children.size() - 1; i >= 0; i--) {
            tasks.push_back({children[i], 0});
            if (i > 0) tasks.push_back({-1, ','});
        }
    }
    return out;
}

static string restricted_newick(const Tree &T, int v,
                                const vector<int> &mark, int token) {
    if (T.max_depth < SAFE_RECURSION_DEPTH)
        return restricted_newick_recursive(T, v, mark, token);
    return restricted_newick_iterative(T, v, mark, token);
}

static void build_lca(Tree &T) {
    T.LOG = 1;
    while ((1 << T.LOG) <= max(1, (int)T.nodes.size())) T.LOG++;
    T.up.assign(T.LOG, vector<int>(T.nodes.size(), -1));
    if (T.root < 0) return;

    T.nodes[T.root].parent = -1;
    T.nodes[T.root].depth = 0;
    T.max_depth = 0;
    vector<int> stack = {T.root};
    vector<int> order;
    order.reserve(T.nodes.size());
    while (!stack.empty()) {
        int v = stack.back();
        stack.pop_back();
        order.push_back(v);
        T.max_depth = max(T.max_depth, T.nodes[v].depth);
        int label = T.nodes[v].label;
        if (label >= 1 && label <= T.n) T.leaf_node[label] = v;
        for (int c : T.nodes[v].child) {
            T.nodes[c].parent = v;
            T.nodes[c].depth = T.nodes[v].depth + 1;
            stack.push_back(c);
        }
    }
    T.subtree_leaves.assign(T.nodes.size(), 0);
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        int v = *it;
        if (T.nodes[v].label >= 1) {
            T.subtree_leaves[v] = 1;
        } else {
            for (int c : T.nodes[v].child)
                T.subtree_leaves[v] += T.subtree_leaves[c];
        }
    }
    T.congestion_prefix.assign(T.nodes.size(), 0);
    for (int v : order) {
        for (int c : T.nodes[v].child) {
            long long s = T.subtree_leaves[c];
            long long demand = max(1LL, s * (T.n - s));
            int edge_price = 1;
            while ((1LL << min(edge_price, 60)) <= demand &&
                   edge_price < 60) {
                ++edge_price;
            }
            T.congestion_prefix[c] =
                T.congestion_prefix[v] + edge_price;
        }
    }
    for (int v = 0; v < (int)T.nodes.size(); v++) T.up[0][v] = T.nodes[v].parent;
    for (int j = 1; j < T.LOG; j++) {
        for (int v = 0; v < (int)T.nodes.size(); v++) {
            int a = T.up[j - 1][v];
            T.up[j][v] = (a < 0 ? -1 : T.up[j - 1][a]);
        }
    }
}

static int lca_node(const Tree &T, int a, int b) {
    if (a < 0 || b < 0) return -1;
    if (T.nodes[a].depth < T.nodes[b].depth) swap(a, b);
    int diff = T.nodes[a].depth - T.nodes[b].depth;
    for (int j = 0; j < T.LOG; j++)
        if (diff & (1 << j)) a = T.up[j][a];
    if (a == b) return a;
    for (int j = T.LOG - 1; j >= 0; j--) {
        if (T.up[j][a] != T.up[j][b]) {
            a = T.up[j][a];
            b = T.up[j][b];
        }
    }
    return T.nodes[a].parent;
}

static int leaf_dist(const Tree &T, int leaf_a, int leaf_b) {
    int a = T.leaf_node[leaf_a], b = T.leaf_node[leaf_b];
    int w = lca_node(T, a, b);
    if (w < 0) return 1 << 28;
    return T.nodes[a].depth + T.nodes[b].depth - 2 * T.nodes[w].depth;
}

static int congestion_dist(const Tree &T, int a, int b) {
    int w = lca_node(T, a, b);
    if (w < 0 || a < 0 || b < 0 ||
        a >= (int)T.congestion_prefix.size() ||
        b >= (int)T.congestion_prefix.size()) {
        return 1 << 28;
    }
    return T.congestion_prefix[a] + T.congestion_prefix[b] -
           2 * T.congestion_prefix[w];
}

static int lca_of_leaf_set(const Tree &T, const vector<int> &S) {
    if (S.empty()) return -1;
    int r = T.leaf_node[S[0]];
    if (r < 0) return -1;
    for (size_t i = 1; i < S.size(); i++) {
        int v = T.leaf_node[S[i]];
        if (v < 0) return -1;
        r = lca_node(T, r, v);
        if (r < 0) return -1;
    }
    return r;
}

// Canonical restricted-tree topology.
static Hash combine_hashes(vector<Hash> &v) {
    sort(v.begin(), v.end());
    Hash h = 1469598103934665603ULL;
    for (Hash x : v) h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
}

static Hash canonical_restriction_hash_recursive(const Tree &T, int v,
                                                 const vector<int> &mark,
                                                 int token) {
    const TreeNode &node = T.nodes[v];
    if (node.label >= 1) {
        if (node.label < (int)mark.size() && mark[node.label] == token)
            return (Hash)(node.label * 1315423911ULL);
        return 0;
    }
    vector<Hash> children;
    for (int c : node.child) {
        Hash h = canonical_restriction_hash_recursive(T, c, mark, token);
        if (h != 0) children.push_back(h);
    }
    if (children.empty()) return 0;
    if (children.size() == 1) return children[0];
    return combine_hashes(children);
}

static Hash canonical_restriction_hash_iterative(const Tree &T, int v,
                                                 const vector<int> &mark,
                                                 int token) {
    static vector<TraversalFrame> stack;
    static vector<Hash> value;
    if (value.size() < T.nodes.size()) value.resize(T.nodes.size());
    stack.clear();
    if (v < 0 || v >= (int)T.nodes.size()) return 0;
    stack.push_back({v, 0});

    while (!stack.empty()) {
        int u = stack.back().node;
        value[u] = 0;
        const TreeNode &node = T.nodes[u];
        if (node.label >= 1) {
            if (node.label < (int)mark.size() && mark[node.label] == token)
                value[u] = (Hash)(node.label * 1315423911ULL);
            stack.pop_back();
            continue;
        }
        if (stack.back().next_child < node.child.size()) {
            int c = node.child[stack.back().next_child++];
            stack.push_back({c, 0});
            continue;
        }

        vector<Hash> children;
        children.reserve(node.child.size());
        for (int c : node.child)
            if (value[c] != 0) children.push_back(value[c]);
        if (children.size() == 1)
            value[u] = children[0];
        else if (children.size() > 1)
            value[u] = combine_hashes(children);
        stack.pop_back();
    }
    return value[v];
}

static Hash canonical_restriction_hash(const Tree &T, int v,
                                       const vector<int> &mark, int token) {
    if (T.max_depth < SAFE_RECURSION_DEPTH)
        return canonical_restriction_hash_recursive(T, v, mark, token);
    return canonical_restriction_hash_iterative(T, v, mark, token);
}

struct VectorHash {
    size_t operator()(const vector<int> &v) const {
        uint64_t h = 1469598103934665603ULL;
        for (int x : v)
            h ^= (uint64_t)x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return (size_t)h;
    }
};

struct CanonicalInterner {
    unordered_map<vector<int>, int, VectorHash> ids;
    int next_id;

    CanonicalInterner(int n, size_t expected_internal_nodes)
        : next_id(n + 1) {
        ids.reserve(expected_internal_nodes * 2 + 1);
    }

    int get(vector<int> &key) {
        auto it = ids.find(key);
        if (it != ids.end()) return it->second;
        int id = next_id++;
        ids.emplace(key, id);
        return id;
    }
};

static int canonical_restriction_id_recursive(const Tree &T, int v,
                                              const vector<int> &mark,
                                              int token,
                                              CanonicalInterner &interner) {
    const TreeNode &node = T.nodes[v];
    if (node.label >= 1) {
        if (node.label < (int)mark.size() && mark[node.label] == token)
            return node.label;
        return 0;
    }
    vector<int> children;
    for (int c : node.child) {
        int id = canonical_restriction_id_recursive(
            T, c, mark, token, interner);
        if (id != 0) children.push_back(id);
    }
    if (children.empty()) return 0;
    if (children.size() == 1) return children[0];
    sort(children.begin(), children.end());
    return interner.get(children);
}

static int canonical_restriction_id_iterative(const Tree &T, int v,
                                              const vector<int> &mark,
                                              int token,
                                              CanonicalInterner &interner) {
    static vector<TraversalFrame> stack;
    static vector<int> value;
    if (value.size() < T.nodes.size()) value.resize(T.nodes.size());
    stack.clear();
    if (v < 0 || v >= (int)T.nodes.size()) return 0;
    stack.push_back({v, 0});

    while (!stack.empty()) {
        int u = stack.back().node;
        value[u] = 0;
        const TreeNode &node = T.nodes[u];
        if (node.label >= 1) {
            if (node.label < (int)mark.size() && mark[node.label] == token)
                value[u] = node.label;
            stack.pop_back();
            continue;
        }
        if (stack.back().next_child < node.child.size()) {
            int c = node.child[stack.back().next_child++];
            stack.push_back({c, 0});
            continue;
        }

        vector<int> children;
        children.reserve(node.child.size());
        for (int c : node.child)
            if (value[c] != 0) children.push_back(value[c]);
        if (children.size() == 1) {
            value[u] = children[0];
        } else if (children.size() > 1) {
            sort(children.begin(), children.end());
            value[u] = interner.get(children);
        }
        stack.pop_back();
    }
    return value[v];
}

static int canonical_restriction_id(const Tree &T, int v,
                                    const vector<int> &mark, int token,
                                    CanonicalInterner &interner) {
    if (T.max_depth < SAFE_RECURSION_DEPTH)
        return canonical_restriction_id_recursive(T, v, mark, token, interner);
    return canonical_restriction_id_iterative(T, v, mark, token, interner);
}

static bool component_topology_valid(const vector<int> &S,
                                     const Tree &T1, const Tree &T2,
                                     int lca1, int lca2) {
    if (S.empty() || lca1 < 0 || lca2 < 0) return false;
    vector<int> mark(T1.n + 1, 0);
    for (int x : S) {
        if (x <= 0 || x > T1.n) return false;
        mark[x] = 1;
    }
    Hash h1 = canonical_restriction_hash(T1, lca1, mark, 1);
    Hash h2 = canonical_restriction_hash(T2, lca2, mark, 1);
    return h1 != 0 && h1 == h2;
}

static bool component_topology_valid_exact_unhashed(
        const vector<int> &S, const Tree &T1, const Tree &T2,
        int lca1, int lca2) {
    if (S.empty() || lca1 < 0 || lca2 < 0) return false;
    vector<int> mark(T1.n + 1, 0);
    for (int x : S) {
        if (x <= 0 || x > T1.n) return false;
        mark[x] = 1;
    }
    CanonicalInterner interner(T1.n, S.size());
    int id1 = canonical_restriction_id(T1, lca1, mark, 1, interner);
    int id2 = canonical_restriction_id(T2, lca2, mark, 1, interner);
    return id1 != 0 && id1 == id2;
}

static bool component_topology_valid_exact(const vector<int> &S,
                                           const Tree &T1, const Tree &T2,
                                           int lca1, int lca2) {
    if (!component_topology_valid(S, T1, T2, lca1, lca2)) return false;
    return component_topology_valid_exact_unhashed(
        S, T1, T2, lca1, lca2);
}

// Agreement forest representation and validity.
struct Component {
    vector<int> leaves;
    bool active = true;
    int repr = -1;
    int root1 = -1;
    int root2 = -1;
};

struct AgreementForest {
    vector<Component> comp;
    vector<int> leaf_comp;
    vector<int> owner1, owner2;
};

// Passive component collection for the final HiGHS recombination stage.
// It does not participate in scoring, acceptance, or restart decisions.
static unordered_set<string> g_h45_component_pool;

static string h45_component_key(const vector<int> &raw_leaves) {
    vector<int> leaves = raw_leaves;
    sort(leaves.begin(), leaves.end());
    string key;
    key.reserve(leaves.size() * 6);
    for (size_t i = 0; i < leaves.size(); i++) {
        if (i) key.push_back(',');
        key += to_string(leaves[i]);
    }
    return key;
}

static void h45_collect_components(const AgreementForest &F) {
    for (const Component &component : F.comp) {
        if (!component.active || component.leaves.empty()) continue;
        g_h45_component_pool.insert(h45_component_key(component.leaves));
    }
}

static void h45_write_component_pool() {
    const char *path = getenv("H45_COMPONENT_POOL");
    if (path == nullptr || path[0] == '\0') return;
    ofstream output(path);
    if (!output) {
        cerr << "h45 component_pool_write_failed path=" << path << '\n';
        return;
    }
    vector<string> keys(g_h45_component_pool.begin(),
                        g_h45_component_pool.end());
    sort(keys.begin(), keys.end());
    for (const string &key : keys) output << key << '\n';
    cerr << "h45 component_pool_written candidates=" << keys.size()
         << " path=" << path << '\n';
}

static int h45_env_int(const char *name, int fallback) {
    const char *value = getenv(name);
    if (value == nullptr || value[0] == '\0') return fallback;
    char *end = nullptr;
    long parsed = strtol(value, &end, 10);
    if (end == value || *end != '\0' || parsed < INT_MIN || parsed > INT_MAX)
        return fallback;
    return (int)parsed;
}

static AgreementForest singleton_forest(int n) {
    AgreementForest F;
    F.leaf_comp.assign(n + 1, -1);
    F.comp.reserve(n);
    for (int x = 1; x <= n; x++) {
        Component c;
        c.leaves.push_back(x);
        c.repr = x;
        F.comp.push_back(c);
        F.leaf_comp[x] = (int)F.comp.size() - 1;
    }
    return F;
}

static int active_count(const AgreementForest &F) {
    int k = 0;
    for (const auto &c : F.comp) if (c.active) k++;
    return k;
}

static bool paths_allowed(const Tree &T, const vector<int> &owner,
                          const vector<int> &S, int root,
                          int allowed_a, int allowed_b) {
    for (int x : S) {
        int v = T.leaf_node[x];
        while (v >= 0 && v != root) {
            int o = owner[v];
            if (o != -1 && o != allowed_a && o != allowed_b) return false;
            v = T.nodes[v].parent;
        }
    }
    return true;
}

static void mark_paths(const Tree &T, vector<int> &owner,
                       const vector<int> &S, int root, int id) {
    for (int x : S) {
        int v = T.leaf_node[x];
        while (v >= 0 && v != root) {
            owner[v] = id;
            v = T.nodes[v].parent;
        }
    }
}

static bool rebuild_owners(AgreementForest &F, const Tree &T1, const Tree &T2) {
    F.owner1.assign(T1.nodes.size(), -1);
    F.owner2.assign(T2.nodes.size(), -1);
    for (int cid = 0; cid < (int)F.comp.size(); cid++) {
        if (!F.comp[cid].active) continue;
        const auto &lv = F.comp[cid].leaves;
        int l1 = lca_of_leaf_set(T1, lv);
        int l2 = lca_of_leaf_set(T2, lv);
        if (l1 < 0 || l2 < 0) return false;
        F.comp[cid].root1 = l1;
        F.comp[cid].root2 = l2;
        if (lv.size() <= 1) continue;
        if (!paths_allowed(T1, F.owner1, lv, l1, cid, cid)) return false;
        if (!paths_allowed(T2, F.owner2, lv, l2, cid, cid)) return false;
        mark_paths(T1, F.owner1, lv, l1, cid);
        mark_paths(T2, F.owner2, lv, l2, cid);
    }
    return true;
}

static int append_component(AgreementForest &R, const vector<int> &leaves, int repr) {
    Component c;
    c.leaves = leaves;
    c.active = true;
    c.repr = repr;
    int id = (int)R.comp.size();
    R.comp.push_back(c);
    for (int x : leaves) R.leaf_comp[x] = id;
    return id;
}

static int append_singleton(AgreementForest &F, int leaf) {
    return append_component(F, vector<int>{leaf}, leaf);
}

static bool try_append_exact_component(AgreementForest &R, const vector<int> &leaves,
                                       int repr,
                                       const Tree &T1, const Tree &T2) {
    if (leaves.empty()) return false;
    if (leaves.size() == 1) {
        append_component(R, leaves, repr);
        return true;
    }

    int l1 = lca_of_leaf_set(T1, leaves);
    int l2 = lca_of_leaf_set(T2, leaves);
    if (l1 < 0 || l2 < 0) return false;
    if (!component_topology_valid_exact(leaves, T1, T2, l1, l2))
        return false;
    if (!paths_allowed(T1, R.owner1, leaves, l1, -2, -2)) return false;
    if (!paths_allowed(T2, R.owner2, leaves, l2, -2, -2)) return false;

    int id = append_component(R, leaves, repr);
    mark_paths(T1, R.owner1, leaves, l1, id);
    mark_paths(T2, R.owner2, leaves, l2, id);
    return true;
}

static AgreementForest exact_output_repair(const AgreementForest &F,
                                  const Tree &T1, const Tree &T2) {
    AgreementForest R;
    R.leaf_comp.assign(T1.n + 1, -1);
    R.owner1.assign(T1.nodes.size(), -1);
    R.owner2.assign(T2.nodes.size(), -1);
    R.comp.reserve(F.comp.size());

    for (const Component &c : F.comp) {
        if (!c.active) continue;
        if (try_append_exact_component(R, c.leaves, c.repr, T1, T2))
            continue;
        for (int x : c.leaves) append_singleton(R, x);
    }
    return R;
}

struct ComponentRegion {
    vector<int> base_ids;
    vector<int> donor_ids;
};

static vector<ComponentRegion> component_overlap_regions(
        const AgreementForest &base, const AgreementForest &donor, int leaf_count) {
    int base_size = (int)base.comp.size();
    int total_size = base_size + (int)donor.comp.size();
    vector<int> parent(total_size);
    for (int i = 0; i < total_size; i++) parent[i] = i;

    auto find_root = [&](int x) {
        int root = x;
        while (parent[root] != root) root = parent[root];
        while (parent[x] != x) {
            int next = parent[x];
            parent[x] = root;
            x = next;
        }
        return root;
    };
    auto unite = [&](int a, int b) {
        int ra = find_root(a);
        int rb = find_root(b);
        if (ra != rb) parent[rb] = ra;
    };

    for (int leaf = 1; leaf <= leaf_count; leaf++) {
        int a = base.leaf_comp[leaf];
        int b = donor.leaf_comp[leaf];
        if (a >= 0 && b >= 0) unite(a, base_size + b);
    }

    unordered_map<int, int> region_index;
    vector<ComponentRegion> regions;
    auto region_for = [&](int vertex) -> ComponentRegion& {
        int root = find_root(vertex);
        auto [it, inserted] =
            region_index.emplace(root, (int)regions.size());
        if (inserted) regions.push_back({});
        return regions[it->second];
    };

    for (int cid = 0; cid < (int)base.comp.size(); cid++)
        if (base.comp[cid].active) region_for(cid).base_ids.push_back(cid);
    for (int cid = 0; cid < (int)donor.comp.size(); cid++)
        if (donor.comp[cid].active)
            region_for(base_size + cid).donor_ids.push_back(cid);
    return regions;
}

static AgreementForest validity_aware_crossover(AgreementForest base, const AgreementForest &donor,
                                       const Tree &T1, const Tree &T2) {
    vector<ComponentRegion> regions =
        component_overlap_regions(base, donor, T1.n);
    sort(regions.begin(), regions.end(),
         [](const ComponentRegion &a, const ComponentRegion &b) {
             int gain_a = (int)a.base_ids.size() - (int)a.donor_ids.size();
             int gain_b = (int)b.base_ids.size() - (int)b.donor_ids.size();
             if (gain_a != gain_b) return gain_a > gain_b;
             if (a.base_ids.size() != b.base_ids.size())
                 return a.base_ids.size() < b.base_ids.size();
             if (a.donor_ids.size() != b.donor_ids.size())
                 return a.donor_ids.size() < b.donor_ids.size();
             int a_base = a.base_ids.empty() ? INT_MAX : a.base_ids.front();
             int b_base = b.base_ids.empty() ? INT_MAX : b.base_ids.front();
             if (a_base != b_base) return a_base < b_base;
             int a_donor =
                 a.donor_ids.empty() ? INT_MAX : a.donor_ids.front();
             int b_donor =
                 b.donor_ids.empty() ? INT_MAX : b.donor_ids.front();
             return a_donor < b_donor;
         });

    vector<unsigned char> switched(regions.size(), 0);

    auto try_region = [&](int index) {
        const ComponentRegion &region = regions[index];
        bool available = true;
        for (int cid : region.base_ids) {
            if (cid < 0 || cid >= (int)base.comp.size() ||
                !base.comp[cid].active) {
                available = false;
                break;
            }
        }
        if (!available) return false;

        AgreementForest trial = base;
        for (int cid : region.base_ids) trial.comp[cid].active = false;
        for (int cid : region.donor_ids) {
            const Component &component = donor.comp[cid];
            append_component(trial, component.leaves, component.repr);
        }
        if (!rebuild_owners(trial, T1, T2)) return false;
        // region.base_ids always refers to the original base components.
        // Donor components appended during accepted swaps receive larger ids,
        // so later availability checks cannot accidentally target them.
        base = std::move(trial);
        switched[index] = 1;
        return true;
    };

    auto strict_pass = [&]() {
        bool changed = true;
        while (changed && !g_stop_requested) {
            changed = false;
            for (int index = 0; index < (int)regions.size(); index++) {
                if (switched[index]) continue;
                const ComponentRegion &region = regions[index];
                int gain =
                    (int)region.base_ids.size() - (int)region.donor_ids.size();
                if (gain <= 0 || !try_region(index)) continue;
                changed = true;
            }
        }
    };

    strict_pass();
    if (!g_stop_requested) {
        vector<int> equal_order(regions.size());
        for (int i = 0; i < (int)regions.size(); i++) equal_order[i] = i;
        for (int index : equal_order) {
            if (switched[index]) continue;
            const ComponentRegion &region = regions[index];
            if (region.base_ids.size() != region.donor_ids.size()) continue;
            if (region.base_ids.size() <= 1) continue;
            try_region(index);
        }
        strict_pass();
    }
    return base;
}

#ifdef H44_DIAG_MAIN
static bool g_audit_incremental_merges = false;
static long long g_audited_incremental_merges = 0;
static long long g_bad_incremental_merges = 0;
#endif

static int try_merge(AgreementForest &F, const Tree &T1, const Tree &T2, int A, int B) {
    if (A < 0 || B < 0 || A == B) return -1;
    if (A >= (int)F.comp.size() || B >= (int)F.comp.size()) return -1;
    if (!F.comp[A].active || !F.comp[B].active) return -1;

    vector<int> merged = F.comp[A].leaves;
    merged.insert(merged.end(), F.comp[B].leaves.begin(), F.comp[B].leaves.end());

    int l1 = lca_of_leaf_set(T1, merged);
    int l2 = lca_of_leaf_set(T2, merged);
    if (!component_topology_valid(merged, T1, T2, l1, l2)) return -1;
    if (!paths_allowed(T1, F.owner1, merged, l1, A, B)) return -1;
    if (!paths_allowed(T2, F.owner2, merged, l2, A, B)) return -1;

    Component c;
    c.leaves.swap(merged);
    c.repr = F.comp[A].repr;
    c.root1 = l1;
    c.root2 = l2;
    c.active = true;
    F.comp[A].active = false;
    F.comp[B].active = false;
    int C = (int)F.comp.size();
    F.comp.push_back(c);
    for (int x : F.comp[C].leaves) F.leaf_comp[x] = C;
    mark_paths(T1, F.owner1, F.comp[C].leaves, l1, C);
    mark_paths(T2, F.owner2, F.comp[C].leaves, l2, C);
#ifdef H44_DIAG_MAIN
    if (g_audit_incremental_merges) {
        g_audited_incremental_merges++;
        AgreementForest rebuilt = F;
        bool ok = rebuild_owners(rebuilt, T1, T2);
        if (!ok || rebuilt.owner1 != F.owner1 || rebuilt.owner2 != F.owner2) {
            g_bad_incremental_merges++;
            if (g_bad_incremental_merges <= 10) {
                cerr << "audit_incremental_owner_mismatch merge="
                     << A << "," << B << " result=" << C
                     << " rebuild_ok=" << ok << '\n';
            }
        }
    }
#endif
    return C;
}

// Greedy merge candidate search.
struct MergeCandidate {
    int score;
    int a, b;
};
struct MergeCandidateGreater {
    bool operator()(const MergeCandidate &x, const MergeCandidate &y) const {
        return x.score > y.score;
    }
};

// libstdc++ and libc++ choose different children while popping equal heap
// entries. That changed the complete greedy route on Linux. Use the libc++
// left-child-on-tie heap rule explicitly so the queue is portable and retains
// the stronger development route.
class MergeQueue {
    vector<MergeCandidate> heap;
    MergeCandidateGreater comp;

public:
    bool empty() const { return heap.empty(); }
    const MergeCandidate &top() const { return heap.front(); }

    void push(const MergeCandidate &value) {
        heap.push_back(value);
        size_t child = heap.size() - 1;
        while (child > 0) {
            size_t parent = (child - 1) / 2;
            if (!comp(heap[parent], heap[child])) break;
            swap(heap[parent], heap[child]);
            child = parent;
        }
    }

    void pop() {
        if (heap.empty()) return;
        MergeCandidate value = heap.back();
        heap.pop_back();
        if (heap.empty()) return;

        size_t hole = 0;
        while (true) {
            size_t left = 2 * hole + 1;
            if (left >= heap.size()) break;
            size_t right = left + 1;
            size_t child =
                right < heap.size() && comp(heap[left], heap[right])
                    ? right
                    : left;
            heap[hole] = heap[child];
            hole = child;
        }
        while (hole > 0) {
            size_t parent = (hole - 1) / 2;
            if (!comp(heap[parent], value)) break;
            heap[hole] = heap[parent];
            hole = parent;
        }
        heap[hole] = value;
    }
};

struct SearchState {
    int jitter = 0;
    int restart_seed_offset = 0;
    int extra_split_count = 0;
    int merge_score_mode = BIG_SCORE_MODE;
    mt19937 rng{1};
};

static int base_merge_score(int a, int b, const Tree &T1, const Tree &T2) {
    int d1 = leaf_dist(T1, a, b);
    int d2 = leaf_dist(T2, a, b);
    return max(d1, d2);
}

static bool component_roots_are_siblings(const Tree &T, int a, int b) {
    if (a < 0 || b < 0 || a == b) return false;
    int pa = T.nodes[a].parent;
    return pa >= 0 && pa == T.nodes[b].parent;
}

static int component_merge_score(const AgreementForest &F, int A, int B,
                                 const Tree &T1, const Tree &T2,
                                 SearchState &search) {
    int r1a = F.comp[A].root1, r1b = F.comp[B].root1;
    int r2a = F.comp[A].root2, r2b = F.comp[B].root2;
    int s;
    if (r1a >= 0 && r1b >= 0 && r2a >= 0 && r2b >= 0) {
        int l1 = lca_node(T1, r1a, r1b);
        int l2 = lca_node(T2, r2a, r2b);
        int d1 = T1.nodes[r1a].depth + T1.nodes[r1b].depth -
                 2 * T1.nodes[l1].depth;
        int d2 = T2.nodes[r2a].depth + T2.nodes[r2b].depth -
                 2 * T2.nodes[l2].depth;
        int max_distance = max(d1, d2);
        int total_distance = d1 + d2;
        int congestion_distance =
            congestion_dist(T1, r1a, r1b) +
            congestion_dist(T2, r2a, r2b);
        int combined_size =
            min(255, (int)F.comp[A].leaves.size() +
                     (int)F.comp[B].leaves.size());
        int mode = search.merge_score_mode;
        if (mode == 5)
            s = base_merge_score(F.comp[A].repr, F.comp[B].repr, T1, T2);
        else if (mode == 1)
            s = total_distance;
        else if (mode == 2)
            s = 32 * max_distance + total_distance;
        else if (mode == 3)
            s = 128 * max_distance + total_distance - combined_size;
        else if (mode == 4)
            s = 128 * total_distance - combined_size;
        else if (mode == 12 || mode == 13)
            s = congestion_distance;
        else if (mode == 6) {
            int common_sibling =
                (component_roots_are_siblings(T1, r1a, r1b) ? 1 : 0) +
                (component_roots_are_siblings(T2, r2a, r2b) ? 1 : 0);
            int imbalance = abs(d1 - d2);
            int local_depth = min(T1.nodes[l1].depth, T2.nodes[l2].depth);
            s = 512 * max_distance +
                32 * imbalance +
                total_distance -
                12 * local_depth -
                3 * combined_size;
            if (common_sibling == 2)
                s -= 1000000;
            else if (common_sibling == 1)
                s -= 128;
        }
        else
            s = max_distance;
    } else {
        s = base_merge_score(F.comp[A].repr, F.comp[B].repr, T1, T2);
    }
    if (search.jitter > 0)
        s += (int)(search.rng() % (unsigned)search.jitter);
    return s;
}

static void local_components_recursive_helper(const Tree &T, const AgreementForest &F,
                                              int v, int limit,
                                              unordered_set<int> &seen,
                                              vector<int> &result) {
    if (g_stop_requested || (int)result.size() >= limit) return;
    const TreeNode &node = T.nodes[v];
    if (node.label >= 1) {
        int cid = F.leaf_comp[node.label];
        if (cid >= 0 && F.comp[cid].active &&
            seen.insert(cid).second)
            result.push_back(cid);
        return;
    }
    for (int c : node.child) {
        if ((int)result.size() >= limit) return;
        local_components_recursive_helper(T, F, c, limit, seen, result);
    }
}

static vector<int> local_components_iterative(const Tree &T, const AgreementForest &F,
                                              int v, int limit) {
    vector<int> result;
    if (v < 0 || v >= (int)T.nodes.size() || limit <= 0) return result;
    unordered_set<int> seen;
    seen.reserve((size_t)min(limit, (int)F.comp.size()) * 2 + 1);

    vector<int> stack = {v};
    while (!stack.empty() && !g_stop_requested &&
           (int)result.size() < limit) {
        int u = stack.back();
        stack.pop_back();
        const TreeNode &nd = T.nodes[u];
        if (nd.label >= 1) {
            int cid = F.leaf_comp[nd.label];
            if (cid >= 0 && F.comp[cid].active &&
                seen.insert(cid).second)
                result.push_back(cid);
            continue;
        }
        for (auto it = nd.child.rbegin(); it != nd.child.rend(); ++it)
            stack.push_back(*it);
    }
    return result;
}

static vector<int> local_components(const Tree &T, const AgreementForest &F,
                                    int v, int limit) {
    vector<int> result;
    if (v < 0 || v >= (int)T.nodes.size() || limit <= 0) return result;
    if (T.max_depth < SAFE_RECURSION_DEPTH) {
        unordered_set<int> seen;
        seen.reserve((size_t)min(limit, (int)F.comp.size()) * 2 + 1);
        local_components_recursive_helper(T, F, v, limit, seen, result);
        return result;
    }
    return local_components_iterative(T, F, v, limit);
}

static uint64_t representative_pair_key(int a, int b) {
    int lo = min(a, b), hi = max(a, b);
    return ((uint64_t)(unsigned)lo << 32) | (unsigned)hi;
}

static void push_cand(MergeQueue &pq, unordered_set<uint64_t> &seen,
                      const AgreementForest &F, const Tree &T1, const Tree &T2,
                      SearchState &search, int A, int B) {
    if (A < 0 || B < 0 || A == B) return;
    if (!F.comp[A].active || !F.comp[B].active) return;
    int a = F.comp[A].repr, b = F.comp[B].repr;
    uint64_t key = representative_pair_key(a, b);
    if (!seen.insert(key).second) return;
    pq.push({component_merge_score(F, A, B, T1, T2, search), a, b});
}

static void seed_candidates(MergeQueue &pq, unordered_set<uint64_t> &seen,
                            const AgreementForest &F, const Tree &T1, const Tree &T2,
                            SearchState &search) {
    for (const Tree *T : {&T1, &T2}) {
        if (g_stop_requested) break;
        for (int v = 0; v < (int)T->nodes.size(); v++) {
            if (g_stop_requested) break;
            maybe_run_periodic_task();
            if (T->nodes[v].child.empty()) continue;
            vector<int> ids =
                local_components(*T, F, v, LOCAL_COMPONENT_WINDOW);
            if (ids.size() < 2) continue;
            for (size_t i = 0; i < ids.size(); i++)
                for (size_t j = i + 1; j < ids.size(); j++)
                    push_cand(pq, seen, F, T1, T2, search, ids[i], ids[j]);
        }
    }
    if (T1.n >= 2) {
        unsigned range = (unsigned)T1.n;
        unsigned mask = 1;
        while (mask < range - 1) mask = (mask << 1) | 1u;
        auto portable_leaf = [&]() {
            unsigned value;
            do {
                value = search.rng() & mask;
            } while (value >= range);
            return (int)value + 1;
        };
        int trials = T1.n * RANDOM_CANDIDATE_MULTIPLIER;
        for (int t = 0; t < trials && !g_stop_requested; t++) {
            maybe_run_periodic_task();
            int a = portable_leaf(), b = portable_leaf();
            if (a == b) continue;
            int A = F.leaf_comp[a], B = F.leaf_comp[b];
            if (A < 0 || B < 0 || A == B) continue;
            push_cand(pq, seen, F, T1, T2, search, A, B);
        }
    }
}

static void push_neighbors(MergeQueue &pq, unordered_set<uint64_t> &seen,
                           unordered_set<uint64_t> &failed,
                           const AgreementForest &F, const Tree &T1, const Tree &T2,
                           SearchState &search, int C) {
    if (C < 0 || !F.comp[C].active) return;
    for (const Tree *T : {&T1, &T2}) {
        if (g_stop_requested) break;
        unordered_set<int> touched;
        for (int leaf : F.comp[C].leaves) {
            if (g_stop_requested) break;
            int v = T->leaf_node[leaf];
            for (int walk = 0;
                 v >= 0 && !g_stop_requested &&
                 walk < NEIGHBOR_WALK_LIMIT;
                 walk++) {
                if (touched.insert(v).second) {
                    vector<int> ids =
                        local_components(*T, F, v, LOCAL_COMPONENT_WINDOW);
                    for (int D : ids) {
                        if (D == C || !F.comp[D].active) continue;
                        uint64_t key = representative_pair_key(F.comp[C].repr,
                                                               F.comp[D].repr);
                        if (failed.erase(key)) seen.erase(key);
                        push_cand(pq, seen, F, T1, T2, search, C, D);
                    }
                }
                v = T->nodes[v].parent;
            }
        }
    }
}

static bool active_candidate_components(const MergeCandidate &cand,
                                        const AgreementForest &F,
                                        int &A, int &B) {
    if (cand.a <= 0 || cand.b <= 0 ||
        cand.a >= (int)F.leaf_comp.size() ||
        cand.b >= (int)F.leaf_comp.size()) {
        return false;
    }
    A = F.leaf_comp[cand.a];
    B = F.leaf_comp[cand.b];
    if (A < 0 || B < 0 || A == B) return false;
    return F.comp[A].active && F.comp[B].active;
}

static bool pop_active_candidate(MergeQueue &pq, const AgreementForest &F,
                                 MergeCandidate &out) {
    while (!pq.empty() && !g_stop_requested) {
        out = pq.top();
        pq.pop();
        int A, B;
        if (active_candidate_components(out, F, A, B)) return true;
    }
    return false;
}

static vector<int> nearby_components(const AgreementForest &F, const Tree &T1,
                                     const Tree &T2, int focus) {
    vector<int> result;
    unordered_set<int> seen;
    if (focus < 0 || focus >= (int)F.comp.size() || !F.comp[focus].active)
        return result;

    for (const Tree *T : {&T1, &T2}) {
        if (g_stop_requested) break;
        unordered_set<int> touched;
        for (int leaf : F.comp[focus].leaves) {
            if (g_stop_requested) break;
            int v = T->leaf_node[leaf];
            while (v >= 0 && !g_stop_requested) {
                if (touched.insert(v).second) {
                    vector<int> ids =
                        local_components(*T, F, v, LOCAL_COMPONENT_WINDOW);
                    for (int D : ids) {
                        if (D == focus || !F.comp[D].active ||
                            !seen.insert(D).second) {
                            continue;
                        }
                        result.push_back(D);
                        if ((int)result.size() >= LOOKAHEAD_NEIGHBOR_LIMIT)
                            return result;
                    }
                }
                v = T->nodes[v].parent;
            }
        }
    }
    return result;
}

static int best_followup_merge(AgreementForest &F, const Tree &T1, const Tree &T2,
                               int focus) {
    vector<int> neighbors = nearby_components(F, T1, T2, focus);
    int best_neighbor = -1;
    int best_score = INT_MAX;
    int best_size = -1;

    for (int D : neighbors) {
        if (g_stop_requested) break;
        AgreementForest trial = F;
        int C = try_merge(trial, T1, T2, focus, D);
        if (C < 0) continue;

        int score = base_merge_score(F.comp[focus].repr, F.comp[D].repr, T1, T2);
        int size = (int)trial.comp[C].leaves.size();
        if (score < best_score ||
            (score == best_score && size > best_size) ||
            (score == best_score && size == best_size &&
             (best_neighbor < 0 ||
              F.comp[D].repr < F.comp[best_neighbor].repr))) {
            best_neighbor = D;
            best_score = score;
            best_size = size;
        }
    }

    if (best_neighbor < 0 || g_stop_requested) return -1;
    return try_merge(F, T1, T2, focus, best_neighbor);
}

static long long lookahead_score(const MergeCandidate &cand, const AgreementForest &F,
                                 const Tree &T1, const Tree &T2) {
    constexpr long long INVALID = (1LL << 60);
    int A, B;
    if (!active_candidate_components(cand, F, A, B)) return INVALID;

    AgreementForest trial = F;
    int focus = try_merge(trial, T1, T2, A, B);
    if (focus < 0) return INVALID;

    int followup = best_followup_merge(trial, T1, T2, focus);
    int followups = followup >= 0 ? 1 : 0;
    if (followup >= 0) focus = followup;
    int final_size = (int)trial.comp[focus].leaves.size();

    return cand.score - 1000000LL * followups - 10000LL * final_size;
}

static bool choose_lookahead_candidate(MergeQueue &pq, const AgreementForest &F,
                                       const Tree &T1, const Tree &T2,
                                       unordered_set<uint64_t> &failed,
                                       MergeCandidate &chosen) {
    constexpr long long INVALID = (1LL << 60);

    while (!pq.empty() && !g_stop_requested) {
        vector<MergeCandidate> batch;
        MergeCandidate cand;
        while ((int)batch.size() < LOOKAHEAD_TOP_CANDIDATES &&
               pop_active_candidate(pq, F, cand)) {
            batch.push_back(cand);
        }
        if (batch.empty()) return false;

        int best_index = -1;
        int best_base_score = batch.front().score;
        long long best_score = INVALID;
        vector<char> keep(batch.size(), 1);
        for (int i = 0; i < (int)batch.size() && !g_stop_requested; i++) {
            if (batch[i].score != best_base_score) continue;
            long long score = lookahead_score(batch[i], F, T1, T2);
            if (score >= INVALID) {
                keep[i] = 0;
                failed.insert(representative_pair_key(batch[i].a, batch[i].b));
                continue;
            }
            if (best_index < 0 || score < best_score ||
                (score == best_score &&
                 batch[i].score < batch[best_index].score) ||
                (score == best_score &&
                 batch[i].score == batch[best_index].score &&
                 batch[i].a < batch[best_index].a) ||
                (score == best_score &&
                 batch[i].score == batch[best_index].score &&
                 batch[i].a == batch[best_index].a &&
                 batch[i].b < batch[best_index].b)) {
                best_index = i;
                best_score = score;
            }
        }

        if (g_stop_requested) return false;
        if (best_index < 0) {
            for (int i = 0; i < (int)batch.size(); i++) {
                if (batch[i].score != best_base_score) pq.push(batch[i]);
            }
            continue;
        }
        for (int i = 0; i < (int)batch.size(); i++) {
            if (i != best_index && keep[i]) pq.push(batch[i]);
        }
        chosen = batch[best_index];
        return true;
    }
    return false;
}

// Common-cherry contraction.
struct SmallComponentSet {
    int ids[3] = {-1, -1, -1};
    int count = 0;
    bool overflow = false;
};

static void small_add(SmallComponentSet &s, int id) {
    if (id < 0 || s.overflow) return;
    for (int i = 0; i < s.count; i++)
        if (s.ids[i] == id) return;
    if (s.count == 3) {
        s.overflow = true;
        return;
    }
    s.ids[s.count++] = id;
}

static void small_merge(SmallComponentSet &a, const SmallComponentSet &b) {
    if (a.overflow) return;
    if (b.overflow) {
        a.overflow = true;
        return;
    }
    for (int i = 0; i < b.count; i++) small_add(a, b.ids[i]);
}

static uint64_t comp_pair_key(int A, int B) {
    int lo = min(A, B), hi = max(A, B);
    return ((uint64_t)(unsigned)lo << 32) | (unsigned)hi;
}

static SmallComponentSet collect_active_cherries_recursive(
        const Tree &T, const AgreementForest &F, int v,
        vector<pair<int, int>> &pairs, unordered_set<uint64_t> &seen) {
    if (g_stop_requested) return {};
    const TreeNode &node = T.nodes[v];
    SmallComponentSet here;

    if (node.label >= 1) {
        int cid = F.leaf_comp[node.label];
        if (cid >= 0 && F.comp[cid].active) small_add(here, cid);
        return here;
    }

    for (int c : node.child) {
        SmallComponentSet child =
            collect_active_cherries_recursive(T, F, c, pairs, seen);
        small_merge(here, child);
    }

    if (!here.overflow && here.count == 2) {
        int A = here.ids[0], B = here.ids[1];
        uint64_t key = comp_pair_key(A, B);
        if (seen.insert(key).second) pairs.push_back({A, B});
    }
    return here;
}

static vector<pair<int, int>> active_cherry_pairs_iterative(const Tree &T,
                                                            const AgreementForest &F) {
    vector<pair<int, int>> pairs;
    unordered_set<uint64_t> seen;
    seen.reserve((size_t)T.n * 2);

    static vector<SmallComponentSet> value;
    static vector<TraversalFrame> stack;
    if (value.size() < T.nodes.size()) value.resize(T.nodes.size());
    stack.clear();
    if (T.root < 0 || T.root >= (int)T.nodes.size()) return pairs;
    stack.push_back({T.root, 0});

    while (!stack.empty() && !g_stop_requested) {
        int v = stack.back().node;
        if (g_stop_requested) break;
        const TreeNode &node = T.nodes[v];
        if (node.label < 1 && stack.back().next_child < node.child.size()) {
            int c = node.child[stack.back().next_child++];
            stack.push_back({c, 0});
            continue;
        }

        SmallComponentSet here;
        if (node.label >= 1) {
            int cid = F.leaf_comp[node.label];
            if (cid >= 0 && F.comp[cid].active) small_add(here, cid);
        } else {
            for (int c : node.child) small_merge(here, value[c]);
        }
        value[v] = here;

        if (!here.overflow && here.count == 2) {
            int A = here.ids[0], B = here.ids[1];
            uint64_t key = comp_pair_key(A, B);
            if (seen.insert(key).second) pairs.push_back({A, B});
        }
        stack.pop_back();
    }
    return pairs;
}

static vector<pair<int, int>> active_cherry_pairs(const Tree &T,
                                                  const AgreementForest &F) {
    if (T.max_depth >= SAFE_RECURSION_DEPTH)
        return active_cherry_pairs_iterative(T, F);

    vector<pair<int, int>> pairs;
    unordered_set<uint64_t> seen;
    seen.reserve((size_t)T.n * 2);
    collect_active_cherries_recursive(T, F, T.root, pairs, seen);
    return pairs;
}

static bool contract_common_cherries(AgreementForest &F, const Tree &T1,
                                     const Tree &T2) {
    bool changed = true;
    bool any_changed = false;
    while (changed && !g_stop_requested) {
        changed = false;

        vector<pair<int, int>> p1 = active_cherry_pairs(T1, F);
        vector<pair<int, int>> p2 = active_cherry_pairs(T2, F);
        unordered_set<uint64_t> in_t2;
        in_t2.reserve(p2.size() * 2 + 1);
        for (auto [A, B] : p2) in_t2.insert(comp_pair_key(A, B));

        for (auto [A, B] : p1) {
            if (g_stop_requested) break;
            if (!in_t2.count(comp_pair_key(A, B))) continue;
            if (try_merge(F, T1, T2, A, B) >= 0) {
                changed = true;
                any_changed = true;
            }
        }
    }
    return any_changed;
}

static bool contract_common_groups(AgreementForest &F, const Tree &T1,
                                   const Tree &T2) {
    return contract_common_cherries(F, T1, T2);
}

static AgreementForest greedy_merge(AgreementForest F, const Tree &T1, const Tree &T2,
                           SearchState &search, bool use_lookahead,
                           bool reconsider_failed) {
    if (!rebuild_owners(F, T1, T2)) return F;
    if (ENABLE_COMMON_CONTRACTION)
        contract_common_groups(F, T1, T2);
    if (g_stop_requested) return F;

    MergeQueue pq;
    unordered_set<uint64_t> seen;
    unordered_set<uint64_t> failed;
    seen.reserve((size_t)T1.n * 8);
    failed.reserve((size_t)T1.n * 2);
    seed_candidates(pq, seen, F, T1, T2, search);

    const int max_checks =
        max(MIN_MERGE_CHECKS, T1.n * MERGE_CHECK_MULTIPLIER);
    const int max_lookahead_choices = use_lookahead ? MAX_LOOKAHEAD_CHOICES : 0;
    int checks = 0;
    int lookahead_choices = 0;

    while (!pq.empty() && !g_stop_requested && checks < max_checks) {
        maybe_run_periodic_task();
        MergeCandidate e;
        if (lookahead_choices < max_lookahead_choices) {
            if (!choose_lookahead_candidate(pq, F, T1, T2, failed, e)) break;
        } else if (!pop_active_candidate(pq, F, e)) {
            break;
        }

        int A, B;
        if (!active_candidate_components(e, F, A, B)) continue;
        checks++;
        int C = try_merge(F, T1, T2, A, B);
        if (C >= 0) {
            failed.erase(representative_pair_key(e.a, e.b));
            if (lookahead_choices < max_lookahead_choices) lookahead_choices++;
            if (ENABLE_COMMON_CONTRACTION &&
                contract_common_groups(F, T1, T2)) {
                pq = MergeQueue();
                seen.clear();
                failed.clear();
                seed_candidates(pq, seen, F, T1, T2, search);
            } else {
                push_neighbors(pq, seen, failed, F, T1, T2, search, C);
            }
        } else {
            if (reconsider_failed) {
                failed.insert(representative_pair_key(e.a, e.b));
            }
        }
    }
    return F;
}

static bool h25_collect_atomic_leaves(
    Node *root, vector<Node *> &leaves) {
    if (root == nullptr) return false;
    vector<Node *> stack = {root};
    while (!stack.empty()) {
        Node *node = stack.back();
        stack.pop_back();
        if (node->get_children().empty()) {
            leaves.push_back(node);
            continue;
        }
        for (Node *child : node->get_children()) stack.push_back(child);
    }
    return true;
}

static string h25_local_newick(
    Node *root, const unordered_map<Node *, int> &local_label) {
    if (root == nullptr) return "";
    if (root->get_children().empty()) {
        auto found = local_label.find(root);
        return found == local_label.end() ? "" : to_string(found->second);
    }

    string result = "(";
    bool first = true;
    for (Node *child : root->get_children()) {
        string child_newick = h25_local_newick(child, local_label);
        if (child_newick.empty()) return "";
        if (!first) result.push_back(',');
        result += child_newick;
        first = false;
    }
    result.push_back(')');
    return result;
}

// Whidden clusters are passed after common-subtree contractions. Contracted
// nodes are intentionally treated as atomic leaves because both trees already
// agree on their hidden topology.
static int h25_local_greedy_distance(
    Forest *first, Forest *second,
    Forest **out_first, Forest **out_second) {
    if (out_first != nullptr) *out_first = nullptr;
    if (out_second != nullptr) *out_second = nullptr;
    if (first == nullptr || second == nullptr ||
        first->num_components() != 1 || second->num_components() != 1) {
        return -1;
    }

    vector<Node *> first_leaves;
    vector<Node *> second_leaves;
    if (!h25_collect_atomic_leaves(first->get_component(0), first_leaves) ||
        !h25_collect_atomic_leaves(second->get_component(0), second_leaves)) {
        return -1;
    }
    if (first_leaves.size() != second_leaves.size() || first_leaves.empty())
        return -1;

    unordered_map<string, int> first_name_count;
    unordered_map<string, vector<Node *>> second_by_name;
    first_name_count.reserve(first_leaves.size() * 2 + 1);
    second_by_name.reserve(second_leaves.size() * 2 + 1);
    for (Node *first_leaf : first_leaves)
        first_name_count[first_leaf->get_name()]++;
    for (Node *second_leaf : second_leaves)
        second_by_name[second_leaf->get_name()].push_back(second_leaf);

    unordered_map<Node *, Node *> second_by_original_first;
    second_by_original_first.reserve(second_leaves.size() * 2 + 1);
    for (Node *second_leaf : second_leaves) {
        Node *original_first = second_leaf->get_twin();
        if (original_first == nullptr ||
            second_by_original_first.count(original_first)) {
            return -1;
        }
        second_by_original_first[original_first] = second_leaf;
    }
    unordered_map<Node *, int> first_label;
    unordered_map<Node *, int> second_label;
    vector<string> first_original(first_leaves.size() + 1);
    vector<string> second_original(first_leaves.size() + 1);
    first_label.reserve(first_leaves.size() * 2 + 1);
    second_label.reserve(second_leaves.size() * 2 + 1);
    for (int i = 0; i < (int)first_leaves.size(); i++) {
        Node *matched_second = nullptr;
        const string name = first_leaves[i]->get_name();
        auto same_name = second_by_name.find(name);
        if (first_name_count[name] == 1 &&
            same_name != second_by_name.end() &&
            same_name->second.size() == 1) {
            matched_second = same_name->second.front();
        } else {
            Node *original_second = first_leaves[i]->get_twin();
            Node *original_first = original_second == nullptr
                ? nullptr
                : original_second->get_twin();
            auto found = second_by_original_first.find(original_first);
            if (original_first != nullptr &&
                found != second_by_original_first.end()) {
                matched_second = found->second;
            }
        }
        if (matched_second == nullptr || second_label.count(matched_second)) {
            return -1;
        }
        first_label[first_leaves[i]] = i + 1;
        second_label[matched_second] = i + 1;
        first_original[i + 1] = first_leaves[i]->get_name();
        second_original[i + 1] = matched_second->get_name();
    }

    string first_newick =
        h25_local_newick(first->get_component(0), first_label);
    string second_newick =
        h25_local_newick(second->get_component(0), second_label);
    if (first_newick.empty() || second_newick.empty()) return -1;

    const int leaf_count = (int)first_leaves.size();
    Tree T1 = parse_newick(first_newick, leaf_count);
    Tree T2 = parse_newick(second_newick, leaf_count);
    if (T1.root < 0 || T2.root < 0) return -1;
    build_lca(T1);
    build_lca(T2);

    SearchState local_search;
    local_search.merge_score_mode = PRIMARY_SINGLETON_SCORE_MODE;
    AgreementForest greedy = greedy_merge(
        singleton_forest(leaf_count), T1, T2, local_search, false, false);
    const int greedy_distance = max(0, active_count(greedy) - 1);

    if (out_first == nullptr || out_second == nullptr)
        return greedy_distance;

    vector<int> component_order;
    component_order.reserve((size_t)active_count(greedy));
    int root_component = -1;
    for (int cid = 0; cid < (int)greedy.comp.size(); cid++) {
        const Component &component = greedy.comp[cid];
        if (!component.active) continue;
        if (root_component < 0 && component.root1 == T1.root &&
            component.root2 == T2.root) {
            root_component = cid;
        }
    }
    if (root_component < 0) {
        for (int cid = 0; cid < (int)greedy.comp.size(); cid++) {
            if (greedy.comp[cid].active &&
                greedy.comp[cid].root1 == T1.root) {
                root_component = cid;
                break;
            }
        }
    }
    if (root_component < 0) {
        for (int cid = 0; cid < (int)greedy.comp.size(); cid++) {
            if (greedy.comp[cid].active) {
                root_component = cid;
                break;
            }
        }
    }
    component_order.push_back(root_component);
    for (int cid = 0; cid < (int)greedy.comp.size(); cid++) {
        if (greedy.comp[cid].active && cid != root_component)
            component_order.push_back(cid);
    }

    vector<int> mark(leaf_count + 1, 0);
    vector<Node *> first_components;
    vector<Node *> second_components;
    first_components.reserve(component_order.size());
    second_components.reserve(component_order.size());
    int token = 0;
    for (int cid : component_order) {
        token++;
        for (int leaf : greedy.comp[cid].leaves) mark[leaf] = token;
        string first_component =
            restricted_newick(T1, T1.root, mark, token);
        string second_component =
            restricted_newick(T2, T2.root, mark, token);
        if (first_component.empty() || second_component.empty()) {
            for (Node *node : first_components) node->delete_tree();
            for (Node *node : second_components) node->delete_tree();
            return greedy_distance;
        }
        first_components.push_back(build_tree(first_component));
        second_components.push_back(build_tree(second_component));
    }

    Forest *native_first = new Forest(first_components);
    Forest *native_second = new Forest(second_components);
    if (!sync_twins(native_first, native_second)) {
        delete native_first;
        delete native_second;
        return greedy_distance;
    }
    for (int i = 0; i < native_first->num_components(); i++) {
        for (Node *leaf : native_first->get_component(i)->find_leaves()) {
            int local = atoi(leaf->get_name().c_str());
            if (local >= 1 && local <= leaf_count)
                leaf->set_name(first_original[local]);
        }
        for (Node *leaf : native_second->get_component(i)->find_leaves()) {
            int local = atoi(leaf->get_name().c_str());
            if (local >= 1 && local <= leaf_count)
                leaf->set_name(second_original[local]);
        }
    }

    *out_first = native_first;
    *out_second = native_second;
    return greedy_distance;
}

static string h25_restricted_local_newick(
    Node *root, const unordered_map<Node *, int> &local_label) {
    if (root == nullptr) return "";
    if (root->get_children().empty()) {
        auto found = local_label.find(root);
        return found == local_label.end() ? "" : to_string(found->second);
    }

    vector<string> children;
    children.reserve(root->get_children().size());
    for (Node *child : root->get_children()) {
        string restricted = h25_restricted_local_newick(child, local_label);
        if (!restricted.empty()) children.push_back(std::move(restricted));
    }
    if (children.empty()) return "";
    if (children.size() == 1) return children.front();

    string result = "(";
    for (int i = 0; i < (int)children.size(); i++) {
        if (i > 0) result.push_back(',');
        result += children[i];
    }
    result.push_back(')');
    return result;
}

// Score a Whidden split candidate using our singleton-greedy algorithm on the
// induced local problem. This only guides subtree selection; Whidden still
// performs every cut and therefore retains its native forest representation.
static int h25_greedy_subtree_score(
    Node *candidate, Forest *first, Forest *second) {
    if (candidate == nullptr || first == nullptr || second == nullptr)
        return -1;

    vector<Node *> first_leaves = candidate->find_leaves();
    if (first_leaves.size() <= 1) return 0;

    vector<Node *> all_first_leaves;
    vector<Node *> all_second_leaves;
    for (int cid = 0; cid < first->num_components(); cid++) {
        vector<Node *> leaves = first->get_component(cid)->find_leaves();
        all_first_leaves.insert(
            all_first_leaves.end(), leaves.begin(), leaves.end());
    }
    for (int cid = 0; cid < second->num_components(); cid++) {
        vector<Node *> leaves = second->get_component(cid)->find_leaves();
        all_second_leaves.insert(
            all_second_leaves.end(), leaves.begin(), leaves.end());
    }

    unordered_map<string, int> first_name_count;
    unordered_map<string, vector<Node *>> second_by_name;
    unordered_map<Node *, Node *> second_by_original_first;
    first_name_count.reserve(all_first_leaves.size() * 2 + 1);
    second_by_name.reserve(all_second_leaves.size() * 2 + 1);
    second_by_original_first.reserve(all_second_leaves.size() * 2 + 1);
    for (Node *leaf : all_first_leaves) first_name_count[leaf->get_name()]++;
    for (Node *leaf : all_second_leaves) {
        second_by_name[leaf->get_name()].push_back(leaf);
        Node *original_first = leaf->get_twin();
        if (original_first != nullptr &&
            !second_by_original_first.count(original_first)) {
            second_by_original_first[original_first] = leaf;
        }
    }

    unordered_map<Node *, int> first_label;
    unordered_map<Node *, int> second_label;
    first_label.reserve(first_leaves.size() * 2 + 1);
    second_label.reserve(first_leaves.size() * 2 + 1);
    for (int i = 0; i < (int)first_leaves.size(); i++) {
        Node *matched_second = nullptr;
        const string name = first_leaves[i]->get_name();
        auto same_name = second_by_name.find(name);
        if (first_name_count[name] == 1 &&
            same_name != second_by_name.end() &&
            same_name->second.size() == 1) {
            matched_second = same_name->second.front();
        } else {
            Node *original_second = first_leaves[i]->get_twin();
            Node *original_first = original_second == nullptr
                ? nullptr
                : original_second->get_twin();
            auto found = second_by_original_first.find(original_first);
            if (original_first != nullptr &&
                found != second_by_original_first.end()) {
                matched_second = found->second;
            }
        }
        if (matched_second == nullptr || second_label.count(matched_second))
            return -1;
        first_label[first_leaves[i]] = i + 1;
        second_label[matched_second] = i + 1;
    }

    string first_newick =
        h25_restricted_local_newick(candidate, first_label);
    vector<string> second_parts;
    second_parts.reserve(second->num_components());
    for (int cid = 0; cid < second->num_components(); cid++) {
        string part = h25_restricted_local_newick(
            second->get_component(cid), second_label);
        if (!part.empty()) second_parts.push_back(std::move(part));
    }
    if (first_newick.empty() || second_parts.empty()) return -1;

    string second_newick;
    if (second_parts.size() == 1) {
        second_newick = std::move(second_parts.front());
    } else {
        second_newick.push_back('(');
        for (int i = 0; i < (int)second_parts.size(); i++) {
            if (i > 0) second_newick.push_back(',');
            second_newick += second_parts[i];
        }
        second_newick.push_back(')');
    }

    const int leaf_count = (int)first_leaves.size();
    Tree T1 = parse_newick(first_newick, leaf_count);
    Tree T2 = parse_newick(second_newick, leaf_count);
    if (T1.root < 0 || T2.root < 0) return -1;
    build_lca(T1);
    build_lca(T2);

    SearchState local_search;
    local_search.merge_score_mode = PRIMARY_SINGLETON_SCORE_MODE;
    AgreementForest greedy = greedy_merge(
        singleton_forest(leaf_count), T1, T2, local_search, false, false);
    return max(0, active_count(greedy) - 1);
}

// Restart destruction and local split approximation.
static void split_to_singletons(AgreementForest &F, const vector<int> &ids) {
    for (int ci : ids) {
        if (ci < 0 || ci >= (int)F.comp.size()) continue;
        if (!F.comp[ci].active || F.comp[ci].leaves.size() <= 1) continue;
        vector<int> leaves = F.comp[ci].leaves;
        F.comp[ci].active = false;
        for (int x : leaves) append_singleton(F, x);
    }
}

static bool exact_component_valid(const vector<int> &leaves,
                                  const Tree &T1, const Tree &T2) {
    if (leaves.empty()) return false;
    if (leaves.size() == 1) return true;
    int l1 = lca_of_leaf_set(T1, leaves);
    int l2 = lca_of_leaf_set(T2, leaves);
    return component_topology_valid_exact(leaves, T1, T2, l1, l2);
}

static void collect_marked_leaves_under(const Tree &T, int root,
                                        const vector<unsigned char> &in_comp,
                                        vector<int> &leaves) {
    leaves.clear();
    if (root < 0 || root >= (int)T.nodes.size()) return;
    vector<int> nodes = {root};
    while (!nodes.empty()) {
        int v = nodes.back();
        nodes.pop_back();
        const TreeNode &node = T.nodes[v];
        if (node.label >= 1) {
            if (node.label <= T.n && in_comp[node.label])
                leaves.push_back(node.label);
        } else {
            for (int c : node.child) nodes.push_back(c);
        }
    }
}

struct SplitCandidate {
    int node;
    int size;
    int depth;
};

static bool component_split_blocks(const Component &comp, const Tree &guide,
                                   const Tree &T1, const Tree &T2,
                                   vector<vector<int>> &blocks,
                                   vector<int> &leftover) {
    blocks.clear();
    leftover.clear();

    int total = (int)comp.leaves.size();
    if (total <= 1 || guide.root < 0) return false;

    vector<unsigned char> in_comp(guide.n + 1, 0);
    for (int x : comp.leaves)
        if (x >= 1 && x <= guide.n) in_comp[x] = 1;

    static vector<TraversalFrame> stack;
    static vector<int> count;
    if (count.size() < guide.nodes.size()) count.resize(guide.nodes.size());
    stack.clear();
    stack.push_back({guide.root, 0});

    vector<SplitCandidate> candidates;
    int max_block = max(2, (total * 3) / 4);
    while (!stack.empty() && !g_stop_requested) {
        int v = stack.back().node;
        const TreeNode &node = guide.nodes[v];
        count[v] = 0;
        if (node.label >= 1) {
            count[v] = (node.label <= guide.n && in_comp[node.label]) ? 1 : 0;
            stack.pop_back();
        } else if (stack.back().next_child < node.child.size()) {
            int c = node.child[stack.back().next_child++];
            stack.push_back({c, 0});
            continue;
        } else {
            for (int c : node.child) count[v] += count[c];
            stack.pop_back();
        }

        int c = count[v];
        if (c >= 2 && c < total && c <= max_block)
            candidates.push_back({v, c, guide.nodes[v].depth});
    }

    sort(candidates.begin(), candidates.end(),
         [](const SplitCandidate &a, const SplitCandidate &b) {
             if (a.size != b.size) return a.size > b.size;
             if (a.depth != b.depth) return a.depth > b.depth;
             return a.node < b.node;
         });

    vector<unsigned char> assigned(T1.n + 1, 0);
    vector<int> leaves;
    for (const SplitCandidate &cand : candidates) {
        if (g_stop_requested || (int)blocks.size() >= MAX_SPLIT_BLOCKS)
            break;
        collect_marked_leaves_under(guide, cand.node, in_comp, leaves);
        if (leaves.size() < 2 || (int)leaves.size() >= total) continue;
        bool overlaps = false;
        for (int x : leaves) {
            if (x < 1 || x > T1.n || assigned[x]) {
                overlaps = true;
                break;
            }
        }
        if (overlaps || !exact_component_valid(leaves, T1, T2)) continue;
        for (int x : leaves) assigned[x] = 1;
        blocks.push_back(leaves);
    }

    for (int x : comp.leaves)
        if (x >= 1 && x <= T1.n && !assigned[x]) leftover.push_back(x);

    return !blocks.empty() && !leftover.empty();
}

static bool split_component_structured(AgreementForest &F, int ci,
                                       const Tree &guide,
                                       const Tree &T1, const Tree &T2) {
    if (ci < 0 || ci >= (int)F.comp.size()) return false;
    if (!F.comp[ci].active || F.comp[ci].leaves.size() <= 1) return false;

    vector<vector<int>> blocks;
    vector<int> leftover;
    Component component = F.comp[ci];
    if (!component_split_blocks(component, guide, T1, T2, blocks, leftover))
        return false;

    F.comp[ci].active = false;
    for (const vector<int> &block : blocks)
        append_component(F, block, block.front());
    for (int x : leftover) append_singleton(F, x);
    return true;
}

static void split_components_for_restart(AgreementForest &F,
                                         const vector<int> &ids,
                                         const Tree &T1, const Tree &T2,
                                         int attempt) {
    if (!STRUCTURED_RESTART_SPLIT) {
        split_to_singletons(F, ids);
        return;
    }

    for (int ci : ids) {
        const Tree &first_guide = ((attempt + ci) & 1) ? T2 : T1;
        const Tree &second_guide = ((attempt + ci) & 1) ? T1 : T2;
        if (!split_component_structured(F, ci, first_guide, T1, T2) &&
            !split_component_structured(F, ci, second_guide, T1, T2)) {
            split_to_singletons(F, vector<int>{ci});
        }
    }
}

static int component_leaf_total(const AgreementForest &F, const vector<int> &ids) {
    int total = 0;
    for (int cid : ids) {
        if (cid >= 0 && cid < (int)F.comp.size() && F.comp[cid].active)
            total += (int)F.comp[cid].leaves.size();
    }
    return total;
}

static int cached_component_root(const AgreementForest &F, const Tree &T, int cid,
                                 vector<int> &cache) {
    if (cid < 0 || cid >= (int)F.comp.size() || !F.comp[cid].active)
        return -1;
    int &root = cache[cid];
    if (root == -2) root = lca_of_leaf_set(T, F.comp[cid].leaves);
    return root;
}

static void score_conflict_split_pair(const AgreementForest &F, const Tree &other,
                                      int A, int B, vector<int> &root_cache,
                                      vector<int> &score) {
    if (A < 0 || B < 0 || A >= (int)F.comp.size() ||
        B >= (int)F.comp.size() || !F.comp[A].active ||
        !F.comp[B].active) {
        return;
    }

    int ra = cached_component_root(F, other, A, root_cache);
    int rb = cached_component_root(F, other, B, root_cache);
    int l = lca_node(other, ra, rb);
    if (l < 0) return;

    vector<int> ids =
        local_components(other, F, l, CONFLICT_COMPONENT_WINDOW);
    bool scored = false;
    for (int cid : ids) {
        if (cid < 0 || cid >= (int)F.comp.size() || !F.comp[cid].active ||
            F.comp[cid].leaves.size() <= 1) {
            continue;
        }
        int bonus = (cid == A || cid == B) ? 4 : 16;
        score[cid] += bonus + min(64, (int)F.comp[cid].leaves.size());
        scored = true;
    }

    if (!scored) {
        if (F.comp[A].leaves.size() > 1)
            score[A] += 4 + min(64, (int)F.comp[A].leaves.size());
        if (F.comp[B].leaves.size() > 1)
            score[B] += 4 + min(64, (int)F.comp[B].leaves.size());
    }
}

static vector<int> collect_conflict_split_candidates(const AgreementForest &F,
                                                     const Tree &T1,
                                                     const Tree &T2) {
    vector<pair<int, int>> p1 = active_cherry_pairs(T1, F);
    vector<pair<int, int>> p2 = active_cherry_pairs(T2, F);
    unordered_set<uint64_t> in_t1, in_t2;
    in_t1.reserve(p1.size() * 2 + 1);
    in_t2.reserve(p2.size() * 2 + 1);
    for (auto [A, B] : p1) in_t1.insert(comp_pair_key(A, B));
    for (auto [A, B] : p2) in_t2.insert(comp_pair_key(A, B));

    vector<int> score(F.comp.size(), 0);
    vector<int> root1(F.comp.size(), -2), root2(F.comp.size(), -2);
    int examined = 0;
    for (auto [A, B] : p1) {
        if (g_stop_requested ||
            examined >= MAX_CONFLICT_PAIRS_PER_TREE) {
            break;
        }
        if (in_t2.count(comp_pair_key(A, B))) continue;
        score_conflict_split_pair(F, T2, A, B, root2, score);
        examined++;
    }
    for (auto [A, B] : p2) {
        if (g_stop_requested ||
            examined >= MAX_CONFLICT_PAIRS_PER_TREE * 2) {
            break;
        }
        if (in_t1.count(comp_pair_key(A, B))) continue;
        score_conflict_split_pair(F, T1, A, B, root1, score);
        examined++;
    }

    vector<int> ids;
    for (int cid = 0; cid < (int)F.comp.size(); cid++) {
        if (score[cid] > 0 && F.comp[cid].active &&
            F.comp[cid].leaves.size() > 1) {
            ids.push_back(cid);
        }
    }
    sort(ids.begin(), ids.end(), [&](int a, int b) {
        if (score[a] != score[b]) return score[a] > score[b];
        if (F.comp[a].leaves.size() != F.comp[b].leaves.size())
            return F.comp[a].leaves.size() > F.comp[b].leaves.size();
        return a < b;
    });
    if (ids.size() > MAX_CONFLICT_POOL_SIZE)
        ids.resize(MAX_CONFLICT_POOL_SIZE);
    return ids;
}

static int find_split_approx_region(const Tree &guide, const AgreementForest &F,
                                    int target_components,
                                    int component_limit,
                                    int leaf_limit) {
    if (guide.root < 0) return -1;
    int node = guide.root;
    while (!g_stop_requested) {
        const TreeNode &nd = guide.nodes[node];
        if (nd.child.empty()) return node;

        int best_child = -1;
        int best_score = 0;
        int best_depth = -1;
        int best_leaves = 0;
        for (int c : nd.child) {
            vector<int> ids = local_components(guide, F, c, component_limit);
            int score = (int)ids.size();
            int leaves = component_leaf_total(F, ids);
            if (score > best_score ||
                (score == best_score && guide.nodes[c].depth > best_depth)) {
                best_child = c;
                best_score = score;
                best_depth = guide.nodes[c].depth;
                best_leaves = leaves;
            }
        }

        if (best_child < 0 || best_score <= 1) return node;
        if (best_score < target_components && best_leaves <= leaf_limit)
            return best_child;
        node = best_child;
    }
    return -1;
}

static bool local_region_components(const Tree &guide, const AgreementForest &F,
                                    int node, vector<int> &ids) {
    ids = local_components(guide, F, node, LOCAL_REGION_COMPONENT_LIMIT + 1);
    if ((int)ids.size() > LOCAL_REGION_COMPONENT_LIMIT) return false;
    return ids.size() >= 2;
}

static bool build_split_approx_trial(const AgreementForest &base,
                                     const vector<int> &ids,
                                     const Tree &guide,
                                     const Tree &T1, const Tree &T2,
                                     AgreementForest &trial_out) {
    Component region;
    region.active = true;
    unordered_set<int> used_leaves;
    used_leaves.reserve(ids.size() * 16 + 1);

    for (int cid : ids) {
        if (cid < 0 || cid >= (int)base.comp.size() || !base.comp[cid].active)
            return false;
        for (int x : base.comp[cid].leaves) {
            if (used_leaves.insert(x).second) region.leaves.push_back(x);
        }
    }
    if (region.leaves.size() <= 1) return false;
    if (region.leaves.size() > SPLIT_APPROX_LEAF_LIMIT) return false;
    region.repr = region.leaves[0];

    vector<vector<int>> blocks;
    vector<int> remainder;
    bool have_blocks =
        component_split_blocks(region, guide, T1, T2, blocks, remainder);

    auto build_trial = [&](bool use_blocks) -> bool {
        AgreementForest trial = base;
        for (int cid : ids) trial.comp[cid].active = false;

        if (use_blocks) {
            for (const vector<int> &block : blocks)
                append_component(trial, block, block[0]);
            for (int x : remainder) append_singleton(trial, x);
        } else {
            for (int x : region.leaves) append_singleton(trial, x);
        }

        if (!rebuild_owners(trial, T1, T2)) return false;
        trial_out = std::move(trial);
        return true;
    };

    if (have_blocks && build_trial(true)) return true;
    return build_trial(false);
}

static bool split_approx_local_improve(AgreementForest &best,
                                       const Tree &T1, const Tree &T2,
                                       SearchState &search,
                                       int attempt) {
    int incumbent_k = active_count(best);

    const Tree *guides[2] = {&T1, &T2};
    for (int g = 0; g < 2 && !g_stop_requested; g++) {
        const Tree &guide = *guides[(attempt + g) & 1];
        int node = find_split_approx_region(guide, best,
                                            SPLIT_APPROX_TARGET_COMPONENTS,
                                            SPLIT_APPROX_COMPONENT_LIMIT,
                                            SPLIT_APPROX_LEAF_LIMIT);
        if (node < 0) continue;

        vector<int> ids;
        if (!local_region_components(guide, best, node, ids)) continue;

        AgreementForest trial;
        if (!build_split_approx_trial(best, ids, guide, T1, T2, trial))
            continue;

        int old_jitter = search.jitter;
        search.jitter = 0;
        AgreementForest repaired = greedy_merge(trial, T1, T2, search, true, true);
        search.jitter = old_jitter;

        int repaired_k = active_count(repaired);
        if (repaired_k < incumbent_k) {
            best = std::move(repaired);
            return true;
        }
    }
    return false;
}

static AgreementForest compact_forest(const AgreementForest &F,
                             const Tree &T1, const Tree &T2) {
    AgreementForest compact;
    compact.leaf_comp.assign(T1.n + 1, -1);
    compact.comp.reserve(active_count(F));
    for (const Component &component : F.comp) {
        if (component.active)
            append_component(compact, component.leaves, component.repr);
    }
    if (!rebuild_owners(compact, T1, T2)) return F;
    return compact;
}

static AgreementForest split_restart(AgreementForest initial, const Tree &T1, const Tree &T2,
                            SearchState &search,
                            const function<void(const AgreementForest &)> &on_improvement,
                            const function<void()> &on_iteration,
                            int max_attempts = 0) {
    AgreementForest best = initial;
    AgreementForest current = std::move(initial);
    int best_k = active_count(best);
    int current_k = best_k;
    int attempt = 0;
    const int old_jitter = search.jitter;
    const int old_merge_score_mode = search.merge_score_mode;
    int restart_seed_offset = search.restart_seed_offset;
    mt19937 pick_rng(987654321u +
                     104729u * (unsigned)restart_seed_offset);
    vector<int> conflict_pool;
    int conflict_pool_current_k = -1;

    while (!g_stop_requested) {
        if (max_attempts > 0 && attempt >= max_attempts) break;
        if (on_iteration) on_iteration();

        vector<int> splittable;
        for (int i = 0; i < (int)current.comp.size(); i++)
            if (current.comp[i].active &&
                current.comp[i].leaves.size() > 1) {
                splittable.push_back(i);
            }
        if (splittable.empty()) break;

        AgreementForest trial = current;
        if (g_stop_requested) break;
        vector<int> split_ids;

        if ((attempt % CONFLICT_RESTART_INTERVAL) == 0) {
            if (conflict_pool_current_k != current_k) {
                conflict_pool =
                    collect_conflict_split_candidates(current, T1, T2);
                conflict_pool_current_k = current_k;
            }
            for (size_t s = 0; s < conflict_pool.size(); s++) {
                int cid =
                    conflict_pool[((attempt / CONFLICT_RESTART_INTERVAL) + s +
                                   (size_t)restart_seed_offset) %
                                  conflict_pool.size()];
                if (cid >= 0 && cid < (int)current.comp.size() &&
                    current.comp[cid].active &&
                    current.comp[cid].leaves.size() > 1) {
                    split_ids.push_back(cid);
                    break;
                }
            }
        }

        if (split_ids.empty()) {
            sort(splittable.begin(), splittable.end(), [&](int a, int b) {
                return current.comp[a].leaves.size() >
                       current.comp[b].leaves.size();
            });
            split_ids.push_back(
                splittable[(attempt + restart_seed_offset) %
                           splittable.size()]);
        }
        if (splittable.size() > 1 &&
            (attempt % CONFLICT_RESTART_INTERVAL) != 0) {
            int j = (int)(pick_rng() % splittable.size());
            if (splittable[j] != split_ids[0])
                split_ids.push_back(splittable[j]);
        }
        for (int extra = 0;
             extra < search.extra_split_count &&
             split_ids.size() < splittable.size();
             extra++) {
            int j = (int)(pick_rng() % splittable.size());
            if (find(split_ids.begin(), split_ids.end(), splittable[j]) ==
                split_ids.end()) {
                split_ids.push_back(splittable[j]);
            }
        }
        split_components_for_restart(trial, split_ids, T1, T2, attempt);

        int mode = attempt % JITTER_CYCLE_LENGTH;
        if      (mode < 8)  search.jitter = 3 + mode;
        else if (mode < 16) search.jitter = 12 + 3 * (mode - 8);
        else                search.jitter = 40 + 8 * (mode - 16);
        search.rng.seed(7654321u +
                        1000003u * (unsigned)restart_seed_offset +
                        65537u * (unsigned)attempt);

        bool use_lookahead =
            (attempt % LOOKAHEAD_RESTART_INTERVAL) ==
            LOOKAHEAD_RESTART_INTERVAL - 1;
        trial = greedy_merge(trial, T1, T2, search,
                             use_lookahead, use_lookahead);
        attempt++;

        int trial_k = active_count(trial);
        h45_collect_components(trial);
        // Equal-size moves let the search cross plateaus. Keep best separate
        // so timeout output can never be worse than the best forest found.
        bool accept_current =
            trial_k <= current_k ||
            (RESTART_ACCEPT_WORSE_MARGIN > 0 &&
             trial_k <= best_k + RESTART_ACCEPT_WORSE_MARGIN);
        if (accept_current) {
            current = std::move(trial);
            current_k = trial_k;
            if (current.comp.size() > (size_t)current_k * 4)
                current = compact_forest(current, T1, T2);
            conflict_pool_current_k = -1;
            if (current_k < best_k) {
                best = current;
                best_k = current_k;
                on_improvement(best);
            }
        }

        if (attempt >= SPLIT_APPROX_START_ATTEMPT &&
            (attempt % SPLIT_APPROX_INTERVAL) == 0 &&
            !g_stop_requested) {
            cerr << "restart_local_repair attempt=" << attempt
                 << " components=" << current_k << '\n';
            if (split_approx_local_improve(
                    current, T1, T2, search, attempt)) {
                current_k = active_count(current);
                conflict_pool_current_k = -1;
                if (current_k < best_k) {
                    best = current;
                    best_k = current_k;
                    on_improvement(best);
                }
            }
        }

        if (on_iteration) on_iteration();
    }

    search.jitter = old_jitter;
    search.merge_score_mode = old_merge_score_mode;
    return best;
}

// Output and entry point.
static string forest_to_output(const AgreementForest &F, const Tree &T1) {
    vector<int> mark(T1.n + 1, 0);
    int token = 0;
    string out;
    for (const Component &c : F.comp) {
        if (!c.active) continue;
        token++;
        for (int x : c.leaves) mark[x] = token;
        string s = restricted_newick(T1, T1.root, mark, token);
        if (!s.empty()) out += s + ";\n";
    }
    return out;
}

static void configure_whidden(int cluster_tune, int split_threshold,
                              int random_seed) {
    CUT_ALL_B = true;
    CUT_ONE_B = true;
    REVERSE_CUT_ONE_B = true;
    REVERSE_CUT_ONE_B_3 = true;
    CUT_TWO_B = true;
    CUT_AC_SEPARATE_COMPONENTS = true;
    EDGE_PROTECTION = true;
    EDGE_PROTECTION_TWO_B = true;
    NEAR_PREORDER_SIBLING_PAIRS = true;
    LEAF_REDUCTION = true;
    LEAF_REDUCTION2 = true;
    PREFER_NONBRANCHING = true;
    APPROX_CUT_ONE_B = true;
    APPROX_CUT_TWO_B = true;
    APPROX_REVERSE_CUT_ONE_B = true;
    DEEPEST_PROTECTED_ORDER = true;
    DEEPEST_ORDER = true;
    PREORDER_SIBLING_PAIRS = true;
    BB = true;
    PREFER_RHO = true;
    CLUSTER_TUNE = cluster_tune;
    SPLIT_APPROX = true;
    SPLIT_APPROX_THRESHOLD = split_threshold;
    srand(random_seed);
}

static bool parse_leaf_label_exact(const string &label, int leaf_count,
                                   int &value) {
    if (label.empty()) return false;
    long long parsed = 0;
    for (unsigned char ch : label) {
        if (!isdigit(ch)) return false;
        parsed = parsed * 10 + (ch - '0');
        if (parsed > leaf_count) return false;
    }
    if (parsed < 1) return false;
    value = (int)parsed;
    return true;
}

static bool component_leaf_labels(Node *component, int leaf_count,
                                  vector<int> &leaves) {
    leaves.clear();
    if (component == nullptr) return false;

    vector<Node *> whidden_leaves = component->find_leaves();
    leaves.reserve(whidden_leaves.size());
    for (Node *leaf : whidden_leaves) {
        int value = 0;
        if (leaf == nullptr ||
            !parse_leaf_label_exact(leaf->get_name(), leaf_count, value)) {
            return false;
        }
        leaves.push_back(value);
    }

    sort(leaves.begin(), leaves.end());
    return adjacent_find(leaves.begin(), leaves.end()) == leaves.end();
}

static bool forest_to_agreement_seed(Forest *source, int leaf_count,
                                     map<int, string> *reverse_label_map,
                                     const Tree &T1, const Tree &T2,
                                     AgreementForest &seed,
                                     bool repair_seed) {
    if (source == nullptr) return false;
    source->numbers_to_labels(reverse_label_map);
    expand_contracted_nodes(source);

    AgreementForest converted;
    converted.leaf_comp.assign(leaf_count + 1, -1);
    vector<unsigned char> seen(leaf_count + 1, 0);
    bool complete = true;
    for (int i = 0; i < source->num_components(); i++) {
        vector<int> leaves;
        if (!component_leaf_labels(source->get_component(i), leaf_count,
                                   leaves)) {
            return false;
        }
        if (leaves.empty()) continue;
        for (int leaf : leaves) {
            if (seen[leaf]) complete = false;
            seen[leaf] = 1;
        }
        append_component(converted, leaves, leaves.front());
    }
    for (int leaf = 1; leaf <= leaf_count; leaf++)
        if (!seen[leaf]) complete = false;
    if (!complete) return false;

    if (repair_seed) {
        converted = exact_output_repair(converted, T1, T2);
        if (!rebuild_owners(converted, T1, T2)) return false;
    }
    seed = std::move(converted);
    return true;
}

static bool whidden_seed_in_process(const TwoTreeInstance &instance,
                                    const Tree &T1, const Tree &T2,
                                    AgreementForest &seed,
                                    bool repair_seed = true,
                                    int cluster_tune = WHIDDEN_CLUSTER_TUNE,
                                    int split_threshold =
                                        WHIDDEN_SPLIT_THRESHOLD,
                                    int random_seed =
                                        WHIDDEN_RANDOM_SEED,
                                    int global_min_k = -1,
                                    int global_max_k = -1) {
    configure_whidden(cluster_tune, split_threshold, random_seed);

    Node *first = build_tree(instance.trees[0]);
    Node *second = build_tree(instance.trees[1]);
    if (first == nullptr || second == nullptr) return false;

    map<string, int> label_map;
    map<int, string> reverse_label_map;
    first->labels_to_numbers(&label_map, &reverse_label_map);
    second->labels_to_numbers(&label_map, &reverse_label_map);

    Forest *out_first = nullptr;
    Forest *out_second = nullptr;
    const int previous_min_spr = MIN_SPR;
    const int previous_max_spr = MAX_SPR;
    const int previous_cluster_max_spr = CLUSTER_MAX_SPR;
    if (global_min_k >= 0) MIN_SPR = global_min_k;
    if (global_max_k >= 0) {
        MAX_SPR = max(MAX_SPR, global_max_k);
        CLUSTER_MAX_SPR = max(CLUSTER_MAX_SPR, global_max_k);
    }
    int whidden_k = rSPR_branch_and_bound_simple_clustering(
        first, second, false, &label_map, &reverse_label_map,
        -1, global_max_k, &out_first, &out_second);
    MIN_SPR = previous_min_spr;
    MAX_SPR = previous_max_spr;
    CLUSTER_MAX_SPR = previous_cluster_max_spr;
    TRACE_PORTFOLIO_MSG("whidden child k " << whidden_k
                         << " out " << (out_first != nullptr));

    first->delete_tree();
    second->delete_tree();
    if (out_first == nullptr) {
        delete out_second;
        return false;
    }

    bool ok = forest_to_agreement_seed(
        out_first, instance.leaf_count, &reverse_label_map,
        T1, T2, seed, repair_seed);
    delete out_first;
    delete out_second;
    return ok;
}

static string serialize_component_forest(const AgreementForest &F) {
    string out;
    for (const Component &component : F.comp) {
        if (!component.active) continue;
        for (size_t i = 0; i < component.leaves.size(); i++) {
            if (i) out.push_back(' ');
            out += to_string(component.leaves[i]);
        }
        out.push_back('\n');
    }
    return out;
}

static bool parse_component_forest_text(const string &text, int n,
                                        const Tree &T1, const Tree &T2,
                                        AgreementForest &F) {
    AgreementForest parsed;
    parsed.leaf_comp.assign(n + 1, -1);
    vector<unsigned char> seen(n + 1, 0);
    stringstream input(text);
    string line;
    while (getline(input, line)) {
        vector<int> leaves;
        for (size_t i = 0; i < line.size();) {
            if (!isdigit((unsigned char)line[i])) {
                i++;
                continue;
            }
            int value = 0;
            while (i < line.size() && isdigit((unsigned char)line[i])) {
                value = value * 10 + (line[i] - '0');
                i++;
            }
            if (value < 1 || value > n || seen[value]) return false;
            seen[value] = 1;
            leaves.push_back(value);
        }
        if (!leaves.empty()) append_component(parsed, leaves, leaves.front());
    }
    for (int leaf = 1; leaf <= n; leaf++)
        if (!seen[leaf]) return false;
    parsed = exact_output_repair(parsed, T1, T2);
    if (!rebuild_owners(parsed, T1, T2)) return false;
    F = std::move(parsed);
    return true;
}

static bool write_all_fd(int fd, const string &data) {
    const char *ptr = data.data();
    size_t remaining = data.size();
    while (remaining > 0) {
        ssize_t written = write(fd, ptr, remaining);
        if (written < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        if (written == 0) return false;
        ptr += written;
        remaining -= (size_t)written;
    }
    return true;
}

static void drain_pipe(int fd, string &buffer) {
    char chunk[8192];
    while (true) {
        ssize_t got = read(fd, chunk, sizeof(chunk));
        if (got > 0) {
            buffer.append(chunk, chunk + got);
            continue;
        }
        if (got == 0) return;
        if (errno == EINTR) continue;
        if (errno == EAGAIN || errno == EWOULDBLOCK) return;
        return;
    }
}

struct WhiddenChild {
    pid_t pid = -1;
    int read_fd = -1;
    string data;
    int status = 0;
    bool running = false;
};

enum class WhiddenPollResult {
    Running,
    Success,
    Failed
};

static bool start_whidden_child(const TwoTreeInstance &instance,
                                const Tree &T1, const Tree &T2,
                                WhiddenChild &child,
                                int cluster_tune,
                                int split_threshold,
                                int random_seed) {
    int pipe_fd[2];
    if (pipe(pipe_fd) != 0) return false;

    pid_t pid = fork();
    if (pid < 0) {
        close(pipe_fd[0]);
        close(pipe_fd[1]);
        return false;
    }

    if (pid == 0) {
        close(pipe_fd[0]);
        signal(SIGTERM, SIG_DFL);
        signal(SIGINT, SIG_DFL);
        g_stop_requested = 0;

        AgreementForest child_seed;
        bool ok = whidden_seed_in_process(
            instance, T1, T2, child_seed, true,
            cluster_tune, split_threshold, random_seed);
        if (ok) {
            string data = serialize_component_forest(child_seed);
            ok = write_all_fd(pipe_fd[1], data);
        }
        close(pipe_fd[1]);
        _exit(ok ? 0 : 2);
    }

    close(pipe_fd[1]);
    int flags = fcntl(pipe_fd[0], F_GETFL, 0);
    if (flags >= 0) fcntl(pipe_fd[0], F_SETFL, flags | O_NONBLOCK);

    child.pid = pid;
    child.read_fd = pipe_fd[0];
    child.data.clear();
    child.status = 0;
    child.running = true;
    return true;
}

static bool start_exact_whidden_seed(const TwoTreeInstance &instance,
                                     const Tree &T1, const Tree &T2,
                                     WhiddenChild &child,
                                     int cluster_tune,
                                     int split_threshold,
                                     int random_seed) {
    return start_whidden_child(instance, T1, T2, child,
                               cluster_tune, split_threshold,
                               random_seed);
}

static void stop_whidden_seed(WhiddenChild &child) {
    if (child.read_fd >= 0) {
        close(child.read_fd);
        child.read_fd = -1;
    }
    if (child.running && child.pid > 0) {
        kill(child.pid, SIGKILL);
        int status = 0;
        while (waitpid(child.pid, &status, 0) < 0 && errno == EINTR) {}
    }
    child.pid = -1;
    child.running = false;
}

static WhiddenPollResult poll_whidden_seed(WhiddenChild &child,
                                           int leaf_count,
                                           const Tree &T1, const Tree &T2,
                                           AgreementForest &seed) {
    if (!child.running) return WhiddenPollResult::Failed;
    if (child.read_fd >= 0) drain_pipe(child.read_fd, child.data);

    int status = 0;
    pid_t result = waitpid(child.pid, &status, WNOHANG);
    if (result == 0) {
        if (g_stop_requested) {
            TRACE_PORTFOLIO_MSG("whidden killed on stop");
            stop_whidden_seed(child);
            return WhiddenPollResult::Failed;
        }
        return WhiddenPollResult::Running;
    }
    if (result < 0) {
        if (errno == EINTR) return WhiddenPollResult::Running;
        TRACE_PORTFOLIO_MSG("whidden waitpid errno " << errno);
        stop_whidden_seed(child);
        return WhiddenPollResult::Failed;
    }

    child.status = status;
    if (child.read_fd >= 0) {
        drain_pipe(child.read_fd, child.data);
        close(child.read_fd);
        child.read_fd = -1;
    }
    child.running = false;
    child.pid = -1;

    if (!WIFEXITED(child.status) || WEXITSTATUS(child.status) != 0) {
        TRACE_PORTFOLIO_MSG("whidden exit status " << child.status
                             << " data " << child.data.size());
        return WhiddenPollResult::Failed;
    }
    if (!parse_component_forest_text(child.data, leaf_count, T1, T2, seed)) {
        TRACE_PORTFOLIO_MSG("whidden parse failed data " << child.data.size());
        return WhiddenPollResult::Failed;
    }
    return WhiddenPollResult::Success;
}

static int legacy_portfolio_main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    TwoTreeInstance instance;
    if (!read_two_tree_instance(cin, instance)) return 1;

    Tree T1 = parse_newick(instance.trees[0], instance.leaf_count);
    Tree T2 = parse_newick(instance.trees[1], instance.leaf_count);
    if (T1.root < 0 || T2.root < 0) return 1;
    build_lca(T1);
    build_lca(T2);

    signal(SIGTERM, handle_signal);
    signal(SIGINT,  handle_signal);

#ifdef RAW_WHIDDEN_ONLY
    AgreementForest raw_whidden;
    if (!whidden_seed_in_process(instance, T1, T2, raw_whidden, false))
        return 1;
    cout << forest_to_output(raw_whidden, T1);
    return 0;
#endif

    WhiddenChild whidden_child;
    WhiddenChild second_whidden_child;
    bool use_second_whidden = SECOND_WHIDDEN_CLUSTER_TUNE > 0;
    bool started_exact_whidden = false;
    bool started_extra_whidden = false;
    if (!g_stop_requested) {
        started_exact_whidden = start_exact_whidden_seed(
            instance, T1, T2, whidden_child,
            DEDICATED_WHIDDEN_CLUSTER_TUNE,
            DEDICATED_WHIDDEN_SPLIT_THRESHOLD,
            WHIDDEN_RANDOM_SEED);
        TRACE_PORTFOLIO_MSG("early whidden start "
                            << started_exact_whidden
                            << " cluster "
                            << DEDICATED_WHIDDEN_CLUSTER_TUNE
                            << " split "
                            << DEDICATED_WHIDDEN_SPLIT_THRESHOLD);
    }

    SearchState singleton_search;
    singleton_search.merge_score_mode = PRIMARY_SINGLETON_SCORE_MODE;
    AgreementForest singleton_raw = greedy_merge(
        singleton_forest(instance.leaf_count), T1, T2,
        singleton_search, false, false);
    AgreementForest singleton =
        exact_output_repair(singleton_raw, T1, T2);
    TRACE_PORTFOLIO_MSG("singleton greedy returned "
                        << active_count(singleton));
    AgreementForest best = singleton;
    AgreementForest best_restart = singleton;
    AgreementForest best_whidden;
    AgreementForest best_crossover;
    bool have_restart = true;
    bool have_whidden = false;
    bool have_crossover = false;

    auto accept_candidate =
        [&](AgreementForest candidate, AgreementForest *bucket,
            bool *have_bucket) {
            candidate = exact_output_repair(candidate, T1, T2);
            int candidate_k = active_count(candidate);
            bool accepted = false;
            if (bucket != nullptr && have_bucket != nullptr &&
                (!*have_bucket || candidate_k < active_count(*bucket))) {
                *bucket = candidate;
                *have_bucket = true;
                accepted = true;
            }
            if (candidate_k < active_count(best)) {
                int previous_k = active_count(best);
                best = std::move(candidate);
                accepted = true;
                TRACE_PORTFOLIO_MSG("global best improved "
                                    << previous_k << " -> "
                                    << candidate_k);
            }
            return accepted;
    };

    int crossover_attempt = 0;
    auto try_crossover_pair =
        [&](const AgreementForest &first, const AgreementForest &second) {
            if (g_stop_requested) return;
            int first_k = active_count(first);
            int second_k = active_count(second);
            const AgreementForest *base = &first;
            const AgreementForest *donor = &second;
            if (second_k < first_k ||
                (second_k == first_k && (crossover_attempt & 1))) {
                base = &second;
                donor = &first;
            }

            AgreementForest crossed =
                validity_aware_crossover(*base, *donor, T1, T2);
            if (g_stop_requested) return;
            SearchState crossover_search;
            crossover_search.merge_score_mode = 1;
            crossover_search.rng.seed(
                24681357u + 131071u * (unsigned)crossover_attempt++);
            crossed = greedy_merge(
                crossed, T1, T2, crossover_search, true, true);
            accept_candidate(
                std::move(crossed), &best_crossover, &have_crossover);
        };

    auto try_pool_crossovers = [&](const AgreementForest &donor) {
        if (g_stop_requested) return;
        if (have_whidden) try_crossover_pair(best_whidden, donor);
        bool had_crossover = have_crossover;
        if (!g_stop_requested && had_crossover)
            try_crossover_pair(best_crossover, donor);
        if (!g_stop_requested && have_whidden && had_crossover)
            try_crossover_pair(best_whidden, best_crossover);
    };

    auto any_whidden_running = [&]() {
        return whidden_child.running ||
               (use_second_whidden && second_whidden_child.running);
    };

    auto run_whidden_restart = [&](const AgreementForest &repaired_donor) {
        if (g_stop_requested || !have_whidden) return;
        SearchState whidden_restart_search;
        whidden_restart_search.merge_score_mode = 1;
        whidden_restart_search.restart_seed_offset =
            WHIDDEN_RESTART_SEED_OFFSET;
        whidden_restart_search.extra_split_count =
            WHIDDEN_EXTRA_SPLITS;
        auto whidden_donor_improved =
            [&](const AgreementForest &raw_donor) {
                accept_candidate(raw_donor, &best_whidden, &have_whidden);
                TRACE_PORTFOLIO_MSG("whidden restart improved "
                                    << active_count(best_whidden));
                try_pool_crossovers(best_restart);
            };
        AgreementForest whidden_restart_result =
            split_restart(repaired_donor, T1, T2,
                          whidden_restart_search,
                          whidden_donor_improved,
                          function<void()>(),
                          WHIDDEN_RESTART_BUDGET);
        accept_candidate(std::move(whidden_restart_result),
                         &best_whidden, &have_whidden);
        try_pool_crossovers(best_restart);
    };

    AgreementForest deferred_whidden_restart;
    bool have_deferred_whidden_restart = false;
    auto defer_whidden_restart =
        [&](const AgreementForest &repaired_donor) {
            if (!have_deferred_whidden_restart ||
                active_count(repaired_donor) <
                    active_count(deferred_whidden_restart)) {
                deferred_whidden_restart = repaired_donor;
                have_deferred_whidden_restart = true;
                TRACE_PORTFOLIO_MSG("whidden restart deferred "
                                    << active_count(deferred_whidden_restart));
            }
    };

    auto flush_deferred_whidden_restart = [&]() {
        if (g_stop_requested || any_whidden_running() ||
            !have_deferred_whidden_restart) {
            return;
        }
        AgreementForest donor = std::move(deferred_whidden_restart);
        have_deferred_whidden_restart = false;
        TRACE_PORTFOLIO_MSG("whidden restart deferred run "
                            << active_count(donor));
        run_whidden_restart(donor);
    };

    auto process_whidden_seed = [&](AgreementForest whidden) {
        int seed_k = active_count(whidden);
        TRACE_PORTFOLIO_MSG(
            "Whidden structured/split-approx returned " << seed_k);
        accept_candidate(whidden, &best_whidden, &have_whidden);
        if (g_stop_requested) return;

        SearchState repair_search;
        repair_search.merge_score_mode = 1;
        AgreementForest repaired =
            greedy_merge(whidden, T1, T2, repair_search, true, true);
        int repaired_k = active_count(repaired);
        TRACE_PORTFOLIO_MSG("greedy on Whidden returned " << repaired_k);
        AgreementForest repaired_donor =
            exact_output_repair(repaired, T1, T2);
        accept_candidate(
            repaired_donor, &best_whidden, &have_whidden);

        if (!g_stop_requested && have_whidden)
            try_crossover_pair(best_whidden, repaired_donor);
        if (!g_stop_requested && have_whidden)
            try_pool_crossovers(best_restart);

        if (!g_stop_requested && have_whidden) {
            if (any_whidden_running())
                defer_whidden_restart(repaired_donor);
            else
                run_whidden_restart(repaired_donor);
        }
    };

    bool dedicated_whidden =
        active_count(singleton) >= DEDICATED_WHIDDEN_SCORE_THRESHOLD;
    if (!dedicated_whidden && whidden_child.running &&
        !KEEP_EARLY_WHIDDEN) {
        TRACE_PORTFOLIO_MSG("early whidden stopped for baseline "
                            << active_count(singleton));
        stop_whidden_seed(whidden_child);
    }
    if (!g_stop_requested && !whidden_child.running) {
        int exact_cluster = dedicated_whidden
            ? DEDICATED_WHIDDEN_CLUSTER_TUNE
            : WHIDDEN_CLUSTER_TUNE;
        int exact_split = dedicated_whidden
            ? DEDICATED_WHIDDEN_SPLIT_THRESHOLD
            : WHIDDEN_SPLIT_THRESHOLD;
        started_exact_whidden = start_exact_whidden_seed(
            instance, T1, T2, whidden_child,
            exact_cluster, exact_split, WHIDDEN_RANDOM_SEED);
        TRACE_PORTFOLIO_MSG("whidden start " << started_exact_whidden
                            << " cluster " << exact_cluster
                            << " split " << exact_split);
    }

    auto poll_one_whidden = [&](WhiddenChild &child) {
        if (!child.running) return;
        AgreementForest whidden;
        WhiddenPollResult result =
            poll_whidden_seed(child, instance.leaf_count,
                              T1, T2, whidden);
        if (result == WhiddenPollResult::Success)
            process_whidden_seed(std::move(whidden));
        else if (result == WhiddenPollResult::Failed)
            TRACE_PORTFOLIO_MSG("whidden failed");
    };

    function<void()> poll_whidden = [&]() {
        poll_one_whidden(whidden_child);
        if (use_second_whidden)
            poll_one_whidden(second_whidden_child);
        flush_deferred_whidden_restart();
    };

    if (!g_stop_requested && dedicated_whidden && started_exact_whidden) {
        TRACE_PORTFOLIO_MSG("dedicated whidden wait start best "
                            << active_count(best));
        for (int tick = 0;
             !g_stop_requested && whidden_child.running &&
             tick < DEDICATED_WHIDDEN_POLL_LIMIT;
             tick++) {
            poll_one_whidden(whidden_child);
            if (whidden_child.running) usleep(50000);
        }
        TRACE_PORTFOLIO_MSG("dedicated whidden wait end best "
                            << active_count(best)
                            << " running " << whidden_child.running);
    }

    auto start_extra_whidden_children = [&]() {
        if (started_extra_whidden || g_stop_requested) return;
        started_extra_whidden = true;
        if (use_second_whidden) {
            bool second_started = start_exact_whidden_seed(
                instance, T1, T2, second_whidden_child,
                SECOND_WHIDDEN_CLUSTER_TUNE,
                SECOND_WHIDDEN_SPLIT_THRESHOLD,
                SECOND_WHIDDEN_RANDOM_SEED);
            TRACE_PORTFOLIO_MSG("whidden2 start " << second_started);
        }
    };
    start_extra_whidden_children();

    if (!g_stop_requested) {
        SearchState restart_search;
        restart_search.merge_score_mode = PRIMARY_RESTART_SCORE_MODE;
        restart_search.extra_split_count = PRIMARY_EXTRA_SPLITS;
        g_periodic_task = &poll_whidden;
        g_periodic_task_tick = 0;
        auto donor_improved = [&](const AgreementForest &raw_donor) {
            accept_candidate(raw_donor, &best_restart, &have_restart);
            TRACE_PORTFOLIO_MSG("restart improved " << active_count(best_restart));
            poll_whidden();
            try_pool_crossovers(best_restart);
        };
        AgreementForest restart_result =
            split_restart(singleton_raw, T1, T2, restart_search,
                          donor_improved, poll_whidden);
        accept_candidate(
            std::move(restart_result), &best_restart, &have_restart);
        poll_whidden();
        try_pool_crossovers(best_restart);
        g_periodic_task = nullptr;
    }
    g_periodic_task = nullptr;

    while (!g_stop_requested && any_whidden_running()) {
        poll_whidden();
        if (any_whidden_running()) usleep(50000);
    }
    if (whidden_child.running) stop_whidden_seed(whidden_child);
    if (use_second_whidden && second_whidden_child.running)
        stop_whidden_seed(second_whidden_child);

    AgreementForest output_forest = exact_output_repair(best, T1, T2);
    cout << forest_to_output(output_forest, T1);
    return 0;
}

#ifndef H44_DIAG_MAIN
static void h45_write_best_snapshot(const AgreementForest &F, const Tree &T1) {
    string serialized = forest_to_output(F, T1);
    update_signal_best_output(serialized);
    const char *path = getenv("H45_BEST_SNAPSHOT");
    if (path == nullptr || path[0] == '\0') return;
    string tmp_path = string(path) + ".tmp";
    ofstream output(tmp_path, ios::binary);
    if (!output) return;
    output << serialized;
    output.close();
    rename(tmp_path.c_str(), path);
}

static void h45_write_previous_best_snapshot(const AgreementForest &F,
                                             const Tree &T1) {
    const char *path = getenv("H45_PREVIOUS_BEST_SNAPSHOT");
    if (path == nullptr || path[0] == '\0') return;
    string tmp_path = string(path) + ".tmp";
    ofstream output(tmp_path, ios::binary);
    if (!output) return;
    output << forest_to_output(F, T1);
    output.close();
    rename(tmp_path.c_str(), path);
}

static vector<int> h45_milestone_thresholds() {
    const char *value = getenv("H45_MILESTONE_THRESHOLDS");
    vector<int> thresholds;
    if (value != nullptr && value[0] != '\0') {
        string text(value);
        size_t start = 0;
        while (start <= text.size()) {
            size_t comma = text.find(',', start);
            string token = text.substr(start, comma == string::npos
                                               ? string::npos
                                               : comma - start);
            if (!token.empty()) {
                char *end = nullptr;
                long parsed = strtol(token.c_str(), &end, 10);
                if (end != token.c_str() && parsed > 0 && parsed < INT_MAX)
                    thresholds.push_back((int)parsed);
            }
            if (comma == string::npos) break;
            start = comma + 1;
        }
    }
    if (thresholds.empty()) thresholds = {430, 428, 427, 426};
    sort(thresholds.begin(), thresholds.end(), greater<int>());
    thresholds.erase(unique(thresholds.begin(), thresholds.end()),
                     thresholds.end());
    return thresholds;
}

static void h45_write_milestone_snapshots(const AgreementForest &F,
                                          const Tree &T1) {
    const char *prefix = getenv("H45_MILESTONE_PREFIX");
    if (prefix == nullptr || prefix[0] == '\0') return;
    int k = active_count(F);
    for (int threshold : h45_milestone_thresholds()) {
        if (k > threshold) continue;
        string path = string(prefix) + "_" + to_string(threshold) + ".out";
        if (access(path.c_str(), F_OK) == 0) continue;
        string tmp_path = path + ".tmp";
        ofstream output(tmp_path, ios::binary);
        if (!output) continue;
        output << forest_to_output(F, T1);
        output.close();
        rename(tmp_path.c_str(), path.c_str());
        cerr << "h45 milestone_written threshold=" << threshold
             << " components=" << k
             << " path=" << path << '\n';
    }
}

static int h45_solve_streams(istream &lane_input, ostream &lane_output) {
    TwoTreeInstance instance;
    if (!read_two_tree_instance(lane_input, instance)) return 1;

    Tree T1 = parse_newick(instance.trees[0], instance.leaf_count);
    Tree T2 = parse_newick(instance.trees[1], instance.leaf_count);
    if (T1.root < 0 || T2.root < 0) return 1;
    build_lca(T1);
    build_lca(T2);

    signal(SIGTERM, handle_signal);
    signal(SIGINT, handle_signal);

    const auto total_start = chrono::steady_clock::now();
    const int common_score_mode =
        h45_env_int("H45_SCORE_MODE", PRIMARY_RESTART_SCORE_MODE);

    SearchState singleton_search;
    singleton_search.merge_score_mode = common_score_mode;
    AgreementForest singleton_raw = greedy_merge(
        singleton_forest(instance.leaf_count), T1, T2,
        singleton_search, false, false);
    AgreementForest singleton = exact_output_repair(singleton_raw, T1, T2);
    AgreementForest best = singleton;
    h45_write_best_snapshot(best, T1);
    h45_write_milestone_snapshots(best, T1);
    AgreementForest previous_best_context;
    bool have_previous_best_context = false;
    auto remember_previous_best = [&]() {
        previous_best_context = best;
        have_previous_best_context = true;
    };
    AgreementForest whidden_repaired_seed;
    bool have_whidden_repaired_seed = false;

    const int singleton_components = active_count(singleton);
    h45_collect_components(singleton);
    const int greedy_ub = max(0, singleton_components - 1);
    cerr << "simple singleton components=" << singleton_components
         << " score_mode=" << common_score_mode << '\n';

    if (!g_stop_requested) {
        H25_GLOBAL_GREEDY_DISTANCE = greedy_ub;
        const auto whidden_start = chrono::steady_clock::now();
        AgreementForest whidden;
        bool whidden_ok = whidden_seed_in_process(
            instance, T1, T2, whidden, false,
            WHIDDEN_CLUSTER_TUNE, WHIDDEN_SPLIT_THRESHOLD,
            WHIDDEN_RANDOM_SEED, -1, -1);
        H25_GLOBAL_GREEDY_DISTANCE = -1;

        const double whidden_seconds = chrono::duration<double>(
            chrono::steady_clock::now() - whidden_start).count();
        if (whidden_ok) {
            whidden = exact_output_repair(whidden, T1, T2);
            h45_collect_components(whidden);
            cerr << "simple whidden components=" << active_count(whidden)
                 << " seconds=" << whidden_seconds << '\n';
            if (active_count(whidden) < active_count(best)) {
                remember_previous_best();
                best = whidden;
                h45_write_best_snapshot(best, T1);
                h45_write_milestone_snapshots(best, T1);
            }

            SearchState repair_search;
            repair_search.merge_score_mode = common_score_mode;
            AgreementForest repaired = greedy_merge(
                whidden, T1, T2, repair_search, true, true);
            repaired = exact_output_repair(repaired, T1, T2);
            h45_collect_components(repaired);
            cerr << "simple whidden_repaired components="
                 << active_count(repaired) << '\n';
            whidden_repaired_seed = repaired;
            have_whidden_repaired_seed = true;
            if (active_count(repaired) < active_count(best)) {
                remember_previous_best();
                best = std::move(repaired);
                h45_write_best_snapshot(best, T1);
                h45_write_milestone_snapshots(best, T1);
            }
        } else {
            cerr << "simple whidden failed seconds="
                 << whidden_seconds << '\n';
        }
    }

    if (!g_stop_requested) {
        SearchState restart_search;
        restart_search.merge_score_mode = common_score_mode;
        restart_search.restart_seed_offset =
            h45_env_int("H45_PRIMARY_SEED_OFFSET",
                        PRIMARY_RESTART_SEED_OFFSET);
        auto on_restart_improvement = [&](const AgreementForest &candidate) {
            AgreementForest repaired = exact_output_repair(candidate, T1, T2);
            h45_collect_components(repaired);
            if (active_count(repaired) < active_count(best)) {
                remember_previous_best();
                best = std::move(repaired);
                const double seconds = chrono::duration<double>(
                    chrono::steady_clock::now() - total_start).count();
                cerr << "simple restart_improvement components="
                     << active_count(best)
                     << " seconds=" << seconds << '\n';
                h45_write_best_snapshot(best, T1);
                h45_write_milestone_snapshots(best, T1);
            }
        };
        AgreementForest restarted = split_restart(
            singleton_raw, T1, T2, restart_search,
            on_restart_improvement, function<void()>(),
            SIMPLE_SINGLETON_RESTART_BUDGET);
        restarted = exact_output_repair(restarted, T1, T2);
        h45_collect_components(restarted);
        if (active_count(restarted) < active_count(best)) {
            remember_previous_best();
            best = std::move(restarted);
            h45_write_best_snapshot(best, T1);
            h45_write_milestone_snapshots(best, T1);
        }
    }

    if (!g_stop_requested) {
        AgreementForest best_seed = best;
        SearchState best_restart_search;
        best_restart_search.merge_score_mode = common_score_mode;
        best_restart_search.restart_seed_offset =
            h45_env_int("H45_BEST_SEED_OFFSET", BEST_RESTART_SEED_OFFSET);
        auto on_best_restart_improvement =
            [&](const AgreementForest &candidate) {
                AgreementForest repaired =
                    exact_output_repair(candidate, T1, T2);
                h45_collect_components(repaired);
                if (active_count(repaired) < active_count(best)) {
                    remember_previous_best();
                    best = std::move(repaired);
                    const double seconds = chrono::duration<double>(
                        chrono::steady_clock::now() - total_start).count();
                    cerr << "simple best_restart_improvement components="
                         << active_count(best)
                         << " seconds=" << seconds << '\n';
                    h45_write_best_snapshot(best, T1);
                    h45_write_milestone_snapshots(best, T1);
                }
        };
        AgreementForest restarted = split_restart(
            best_seed, T1, T2, best_restart_search,
            on_best_restart_improvement, function<void()>(), 0);
        restarted = exact_output_repair(restarted, T1, T2);
        h45_collect_components(restarted);
        if (active_count(restarted) < active_count(best)) {
            remember_previous_best();
            best = std::move(restarted);
            h45_write_best_snapshot(best, T1);
            h45_write_milestone_snapshots(best, T1);
        }
    }

    AgreementForest output_forest = exact_output_repair(best, T1, T2);
    h45_write_best_snapshot(output_forest, T1);
    h45_write_milestone_snapshots(output_forest, T1);
    if (have_previous_best_context)
        h45_write_previous_best_snapshot(previous_best_context, T1);
    h45_collect_components(output_forest);
    h45_write_component_pool();
    cerr << "simple final components=" << active_count(output_forest)
         << '\n';
    lane_output << forest_to_output(output_forest, T1);
    return 0;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    return h45_solve_streams(cin, cout);
}
#endif

#ifdef H44_DIAG_MAIN
static string read_whole_file(const string &path) {
    ifstream in(path);
    stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

static string strip_comment_lines(const string &text) {
    stringstream input(text);
    stringstream output;
    string line;
    while (getline(input, line)) {
        size_t p = line.find_first_not_of(" \t\r\n");
        if (p != string::npos && line[p] == '#') continue;
        output << line << '\n';
    }
    return output.str();
}

static bool read_forest_from_path(const string &path, int n,
                                  const Tree &T1, const Tree &T2,
                                  AgreementForest &F) {
    string forest_text = strip_comment_lines(read_whole_file(path));
    if (!parse_component_forest_text(forest_text, n, T1, T2, F))
        return false;
    F = exact_output_repair(F, T1, T2);
    return rebuild_owners(F, T1, T2);
}

static bool write_candidate_output(const AgreementForest &F,
                                   const Tree &T1,
                                   const string &path) {
    ofstream out(path);
    if (!out) return false;
    out << forest_to_output(F, T1);
    return true;
}

static bool try_write_if_valid(AgreementForest trial,
                               const Tree &T1, const Tree &T2,
                               int before, const string &path,
                               const string &why,
                               int max_after = -1) {
    if (!rebuild_owners(trial, T1, T2)) return false;
    trial = exact_output_repair(trial, T1, T2);
    if (!rebuild_owners(trial, T1, T2)) return false;
    int after = active_count(trial);
    if (after >= before) return false;
    if (max_after >= 0 && after > max_after) return false;
    cerr << "diag_improvement " << why << " " << before
         << " -> " << after << "\n";
    write_candidate_output(trial, T1, path);
    return true;
}

static bool try_local_exact_exchange(const AgreementForest &F,
                                     const Tree &T1, const Tree &T2,
                                     const vector<int> &ids,
                                     int max_leaves,
                                     int before,
                                     const string &out_path,
                                     const string &why,
                                     int max_after = -1) {
    vector<int> leaves;
    for (int cid : ids) {
        if (cid < 0 || cid >= (int)F.comp.size() || !F.comp[cid].active)
            return false;
        leaves.insert(leaves.end(), F.comp[cid].leaves.begin(),
                      F.comp[cid].leaves.end());
    }
    sort(leaves.begin(), leaves.end());
    leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());

    int m = (int)leaves.size();
    int target = (int)ids.size() - 1;
    if (target <= 0 || m < target || m > max_leaves || m >= 62)
        return false;

    AgreementForest base = F;
    for (int cid : ids) base.comp[cid].active = false;
    if (!rebuild_owners(base, T1, T2)) return false;

    struct CandidateSubset {
        uint64_t mask;
        vector<int> leaves;
    };

    vector<CandidateSubset> candidates;
    vector<vector<int>> by_first(m);
    uint64_t full = (uint64_t(1) << m) - 1;

    for (uint64_t mask = 1; mask <= full; mask++) {
        vector<int> subset;
        subset.reserve(__builtin_popcountll(mask));
        int first = -1;
        for (int bit = 0; bit < m; bit++) {
            if ((mask & (uint64_t(1) << bit)) == 0) continue;
            if (first < 0) first = bit;
            subset.push_back(leaves[bit]);
        }
        if (!exact_component_valid(subset, T1, T2)) continue;

        AgreementForest test = base;
        if (!try_append_exact_component(test, subset, subset.front(), T1, T2))
            continue;
        by_first[first].push_back((int)candidates.size());
        candidates.push_back({mask, std::move(subset)});
    }

    if ((int)candidates.size() < target) return false;

    function<bool(uint64_t, AgreementForest&, int)> dfs =
        [&](uint64_t covered, AgreementForest &trial, int used) -> bool {
            if (covered == full) {
                if (used == target) {
                    return try_write_if_valid(trial, T1, T2, before, out_path,
                                              why, max_after);
                }
                return false;
            }
            if (used >= target) return false;

            int first = 0;
            while (first < m && (covered & (uint64_t(1) << first))) first++;
            if (first >= m) return false;

            for (int idx : by_first[first]) {
                const CandidateSubset &cand = candidates[idx];
                if (cand.mask & covered) continue;
                AgreementForest next = trial;
                if (!try_append_exact_component(next, cand.leaves,
                                                cand.leaves.front(),
                                                T1, T2)) {
                    continue;
                }
                if (dfs(covered | cand.mask, next, used + 1)) {
                    cerr << "diag_region_components " << ids.size()
                         << " diag_region_leaves " << m << "\n";
                    return true;
                }
            }
            return false;
        };

    return dfs(0, base, 0);
}

static int audit_merge_validator(const AgreementForest &input,
                                 const Tree &T1, const Tree &T2,
                                 const string &out_path) {
    AgreementForest F = input;
    AgreementForest rebuilt = F;
    bool initial_rebuild_ok = rebuild_owners(rebuilt, T1, T2);
    bool owners_match = initial_rebuild_ok &&
        F.owner1 == rebuilt.owner1 && F.owner2 == rebuilt.owner2;

    long long pairs = 0;
    long long fast_topology_ok = 0;
    long long exact_topology_ok = 0;
    long long hash_false_negative = 0;
    long long state_paths_ok = 0;
    long long fresh_forest_ok = 0;
    long long try_merge_ok = 0;
    long long false_rejections = 0;
    long long false_acceptances = 0;
    long long incremental_owner_failures = 0;
    bool wrote_topology_only_candidate = false;

    vector<int> ids;
    for (int cid = 0; cid < (int)F.comp.size(); cid++)
        if (F.comp[cid].active) ids.push_back(cid);

    for (size_t i = 0; i < ids.size(); i++) {
        for (size_t j = i + 1; j < ids.size(); j++) {
            int A = ids[i], B = ids[j];
            pairs++;
            vector<int> merged = F.comp[A].leaves;
            merged.insert(merged.end(), F.comp[B].leaves.begin(),
                          F.comp[B].leaves.end());
            int l1 = lca_of_leaf_set(T1, merged);
            int l2 = lca_of_leaf_set(T2, merged);
            bool fast = component_topology_valid(
                merged, T1, T2, l1, l2);
            bool exact = component_topology_valid_exact_unhashed(
                merged, T1, T2, l1, l2);
            fast_topology_ok += fast;
            exact_topology_ok += exact;
            if (!fast && exact) hash_false_negative++;

            bool state_paths = exact &&
                paths_allowed(T1, F.owner1, merged, l1, A, B) &&
                paths_allowed(T2, F.owner2, merged, l2, A, B);
            state_paths_ok += state_paths;

            AgreementForest fresh_trial = F;
            fresh_trial.comp[A].active = false;
            fresh_trial.comp[B].active = false;
            append_component(fresh_trial, merged, F.comp[A].repr);
            bool fresh_ok = exact && rebuild_owners(fresh_trial, T1, T2);
            fresh_forest_ok += fresh_ok;
            if (exact && !fresh_ok && !wrote_topology_only_candidate) {
                write_candidate_output(fresh_trial, T1, out_path);
                wrote_topology_only_candidate = true;
                cerr << "audit_topology_only_candidate A=" << A
                     << " B=" << B
                     << " leaves=" << merged.size() << '\n';
            }

            AgreementForest incremental_trial = F;
            bool merged_ok =
                try_merge(incremental_trial, T1, T2, A, B) >= 0;
            try_merge_ok += merged_ok;

            bool incremental_rebuild_ok = false;
            if (merged_ok) {
                AgreementForest check = incremental_trial;
                incremental_rebuild_ok = rebuild_owners(check, T1, T2);
                if (!incremental_rebuild_ok ||
                    check.owner1 != incremental_trial.owner1 ||
                    check.owner2 != incremental_trial.owner2) {
                    incremental_owner_failures++;
                }
            }

            if (fresh_ok && !merged_ok) {
                false_rejections++;
                cerr << "audit_false_rejection A=" << A
                     << " B=" << B
                     << " leaves=" << merged.size()
                     << " fast=" << fast
                     << " state_paths=" << state_paths << '\n';
                if (false_rejections == 1)
                    write_candidate_output(fresh_trial, T1, out_path);
            }
            if (merged_ok && !fresh_ok) {
                false_acceptances++;
                cerr << "audit_false_acceptance A=" << A
                     << " B=" << B
                     << " leaves=" << merged.size()
                     << " exact=" << exact
                     << " incremental_rebuild=" << incremental_rebuild_ok
                     << '\n';
            }
        }
    }

    cerr << "audit_initial_rebuild_ok " << initial_rebuild_ok << '\n';
    cerr << "audit_initial_owners_match " << owners_match << '\n';
    cerr << "audit_pairs " << pairs << '\n';
    cerr << "audit_fast_topology_ok " << fast_topology_ok << '\n';
    cerr << "audit_exact_topology_ok " << exact_topology_ok << '\n';
    cerr << "audit_hash_false_negative " << hash_false_negative << '\n';
    cerr << "audit_state_paths_ok " << state_paths_ok << '\n';
    cerr << "audit_fresh_forest_ok " << fresh_forest_ok << '\n';
    cerr << "audit_try_merge_ok " << try_merge_ok << '\n';
    cerr << "audit_false_rejections " << false_rejections << '\n';
    cerr << "audit_false_acceptances " << false_acceptances << '\n';
    cerr << "audit_incremental_owner_failures "
         << incremental_owner_failures << '\n';
    return (false_rejections || false_acceptances ||
            incremental_owner_failures || !owners_match) ? 1 : 0;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        cerr << "usage: h44_diag INSTANCE FOREST_OUT CANDIDATE_OUT\n"
             << "   or: h44_diag INSTANCE --audit FOREST_OUT CANDIDATE_OUT\n"
             << "   or: h44_diag INSTANCE --audit-greedy CANDIDATE_OUT\n";
        return 2;
    }

    ifstream instance_in(argv[1]);
    TwoTreeInstance instance;
    if (!read_two_tree_instance(instance_in, instance)) return 1;
    Tree T1 = parse_newick(instance.trees[0], instance.leaf_count);
    Tree T2 = parse_newick(instance.trees[1], instance.leaf_count);
    build_lca(T1);
    build_lca(T2);

    if (string(argv[2]) == "--audit-greedy") {
        g_audit_incremental_merges = true;
        SearchState search;
        search.merge_score_mode = 0;
        AgreementForest F = greedy_merge(
            singleton_forest(instance.leaf_count), T1, T2,
            search, false, false);
        g_audit_incremental_merges = false;
        AgreementForest rebuilt = F;
        bool rebuild_ok = rebuild_owners(rebuilt, T1, T2);
        bool owners_match = rebuild_ok &&
            rebuilt.owner1 == F.owner1 && rebuilt.owner2 == F.owner2;
        cerr << "audit_greedy_components " << active_count(F) << '\n';
        cerr << "audit_greedy_merges " << g_audited_incremental_merges
             << '\n';
        cerr << "audit_greedy_bad_incremental_merges "
             << g_bad_incremental_merges << '\n';
        cerr << "audit_greedy_final_rebuild_ok " << rebuild_ok << '\n';
        cerr << "audit_greedy_final_owners_match " << owners_match << '\n';
        write_candidate_output(F, T1, argv[3]);
        return (g_bad_incremental_merges || !owners_match) ? 1 : 0;
    }

    if (argc >= 5 && string(argv[2]) == "--audit") {
        AgreementForest F;
        if (!read_forest_from_path(argv[3], instance.leaf_count,
                                   T1, T2, F)) {
            cerr << "audit_parse_forest_failed\n";
            return 1;
        }
        cerr << "audit_components " << active_count(F) << '\n';
        return audit_merge_validator(F, T1, T2, argv[4]);
    }

    if (argc >= 5 && string(argv[2]) == "--pool") {
        string out_path = argv[3];
        vector<pair<string, AgreementForest>> pool;
        for (int arg = 4; arg < argc; arg++) {
            AgreementForest F;
            if (!read_forest_from_path(argv[arg], instance.leaf_count,
                                       T1, T2, F)) {
                cerr << "pool_parse_failed " << argv[arg] << "\n";
                continue;
            }
            cerr << "pool_loaded " << active_count(F) << " "
                 << argv[arg] << "\n";
            pool.push_back({argv[arg], std::move(F)});
        }
        if (pool.empty()) {
            cerr << "pool_empty\n";
            return 1;
        }

        int best_k = INT_MAX;
        AgreementForest best;
        string best_name;
        for (auto &entry : pool) {
            int k = active_count(entry.second);
            if (k < best_k) {
                best_k = k;
                best = entry.second;
                best_name = entry.first;
            }
        }
        cerr << "pool_initial_best " << best_k << " " << best_name << "\n";
        if (best_k <= 424) {
            write_candidate_output(best, T1, out_path);
            return 0;
        }

        int attempt = 0;
        for (int i = 0; i < (int)pool.size(); i++) {
            for (int j = i + 1; j < (int)pool.size(); j++) {
                for (int dir = 0; dir < 2; dir++) {
                    if (g_stop_requested) break;
                    const AgreementForest &base =
                        dir == 0 ? pool[i].second : pool[j].second;
                    const AgreementForest &donor =
                        dir == 0 ? pool[j].second : pool[i].second;
                    AgreementForest crossed =
                        validity_aware_crossover(base, donor, T1, T2);
                    if (!rebuild_owners(crossed, T1, T2)) continue;

                    for (int mode : {0, 1, 5}) {
                        SearchState search;
                        search.merge_score_mode = mode;
                        search.rng.seed(911382323u +
                                        131071u * (unsigned)attempt++);
                        AgreementForest repaired =
                            greedy_merge(crossed, T1, T2, search, true, true);
                        repaired = exact_output_repair(repaired, T1, T2);
                        int k = active_count(repaired);
                        if (k < best_k) {
                            best_k = k;
                            best = repaired;
                            cerr << "pool_improvement " << best_k
                                 << " pair=" << i << "," << j
                                 << " dir=" << dir
                                 << " mode=" << mode << "\n";
                            write_candidate_output(best, T1, out_path);
                            if (best_k <= 424) return 0;
                        }
                    }
                }
            }
        }

        struct PoolRegionTask {
            int base_index;
            vector<int> ids;
            int leaves;
            int cap;
            string why;
        };
        vector<PoolRegionTask> region_tasks;
        unordered_set<string> pool_region_seen;
        for (int i = 0; i < (int)pool.size(); i++) {
            for (int j = 0; j < (int)pool.size(); j++) {
                if (i == j) continue;
                vector<ComponentRegion> regions =
                    component_overlap_regions(pool[i].second,
                                              pool[j].second,
                                              instance.leaf_count);
                for (const ComponentRegion &region : regions) {
                    int k = (int)region.base_ids.size();
                    if (k < 3 || k > 8) continue;
                    int cap = 0;
                    if (k == 3) cap = 28;
                    else if (k == 4) cap = 26;
                    else if (k == 5) cap = 24;
                    else if (k == 6) cap = 22;
                    else if (k == 7) cap = 18;
                    else cap = 16;

                    int leaves = component_leaf_total(pool[i].second,
                                                      region.base_ids);
                    if (leaves > cap) continue;
                    vector<int> ids = region.base_ids;
                    sort(ids.begin(), ids.end());
                    string key = to_string(i) + ":";
                    for (int id : ids) key += to_string(id) + ",";
                    if (!pool_region_seen.insert(key).second) continue;
                    region_tasks.push_back(
                        {i, ids, leaves, cap,
                         "pool_overlap_" + to_string(k) + "_to_" +
                         to_string(k - 1)});
                }
            }
        }
        sort(region_tasks.begin(), region_tasks.end(),
             [](const PoolRegionTask &a, const PoolRegionTask &b) {
                 if (a.leaves != b.leaves) return a.leaves < b.leaves;
                 if (a.ids.size() != b.ids.size())
                     return a.ids.size() > b.ids.size();
                 if (a.base_index != b.base_index)
                     return a.base_index < b.base_index;
                 return a.ids < b.ids;
             });
        cerr << "pool_overlap_tasks " << region_tasks.size() << "\n";
        int overlap_limit = min<int>((int)region_tasks.size(), 3000);
        for (int idx = 0; idx < overlap_limit; idx++) {
            if (idx > 0 && idx % 500 == 0)
                cerr << "pool_overlap_checked " << idx << "\n";
            const PoolRegionTask &task = region_tasks[idx];
            if (try_local_exact_exchange(pool[task.base_index].second,
                                         T1, T2, task.ids, task.cap,
                                         active_count(pool[task.base_index].second),
                                         out_path, task.why,
                                         best_k - 1)) {
                return 0;
            }
        }

        cerr << "pool_final_best " << best_k << "\n";
        write_candidate_output(best, T1, out_path);
        return 0;
    }

    AgreementForest F;
    string forest_text = strip_comment_lines(read_whole_file(argv[2]));
    if (!parse_component_forest_text(
            forest_text, instance.leaf_count, T1, T2, F)) {
        cerr << "diag_parse_forest_failed\n";
        return 1;
    }
    int before = active_count(F);
    cerr << "diag_components " << before << "\n";

    vector<int> active_ids;
    for (int cid = 0; cid < (int)F.comp.size(); cid++)
        if (F.comp[cid].active) active_ids.push_back(cid);

    for (size_t i = 0; i < active_ids.size(); i++) {
        for (size_t j = i + 1; j < active_ids.size(); j++) {
            AgreementForest trial = F;
            int C = try_merge(trial, T1, T2, active_ids[i], active_ids[j]);
            if (C >= 0 &&
                try_write_if_valid(trial, T1, T2, before, argv[3],
                                   "pair")) {
                return 0;
            }
        }
    }
    cerr << "diag_no_pair_merge\n";

    unordered_set<string> seen_regions;
    auto region_key = [](vector<int> ids) {
        sort(ids.begin(), ids.end());
        string key;
        for (int id : ids) key += to_string(id) + ",";
        return key;
    };

    for (const Tree *guide : {&T1, &T2}) {
        for (int node = 0; node < (int)guide->nodes.size(); node++) {
            vector<int> ids = local_components(*guide, F, node, 4);
            if (ids.size() != 3) continue;
            string key = region_key(ids);
            if (!seen_regions.insert(key).second) continue;

            vector<int> leaves;
            for (int cid : ids) {
                leaves.insert(leaves.end(),
                              F.comp[cid].leaves.begin(),
                              F.comp[cid].leaves.end());
            }
            sort(leaves.begin(), leaves.end());
            leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());
            int m = (int)leaves.size();
            if (m < 3 || m > 22) continue;

            int full = (1 << m) - 1;
            for (int mask = 1; mask < full; mask++) {
                if ((mask & 1) == 0) continue;
                vector<int> A, B;
                for (int bit = 0; bit < m; bit++) {
                    if (mask & (1 << bit)) A.push_back(leaves[bit]);
                    else B.push_back(leaves[bit]);
                }
                if (A.empty() || B.empty()) continue;
                if (!exact_component_valid(A, T1, T2) ||
                    !exact_component_valid(B, T1, T2)) {
                    continue;
                }

                AgreementForest trial = F;
                for (int cid : ids) trial.comp[cid].active = false;
                append_component(trial, A, A.front());
                append_component(trial, B, B.front());
                if (try_write_if_valid(trial, T1, T2, before, argv[3],
                                       "three_to_two")) {
                    cerr << "diag_region_leaves " << m << "\n";
                    return 0;
                }
            }
        }
    }
    cerr << "diag_no_three_to_two\n";

    seen_regions.clear();
    for (const Tree *guide : {&T1, &T2}) {
        for (int node = 0; node < (int)guide->nodes.size(); node++) {
            vector<int> ids = local_components(*guide, F, node, 5);
            if (ids.size() != 4) continue;
            string key = region_key(ids);
            if (!seen_regions.insert(key).second) continue;
            if (try_local_exact_exchange(F, T1, T2, ids, 22, before,
                                         argv[3], "four_to_three")) {
                return 0;
            }
        }
    }
    cerr << "diag_no_four_to_three\n";

    seen_regions.clear();
    for (const Tree *guide : {&T1, &T2}) {
        for (int node = 0; node < (int)guide->nodes.size(); node++) {
            vector<int> ids = local_components(*guide, F, node, 6);
            if (ids.size() != 5) continue;
            string key = region_key(ids);
            if (!seen_regions.insert(key).second) continue;
            if (try_local_exact_exchange(F, T1, T2, ids, 18, before,
                                         argv[3], "five_to_four")) {
                return 0;
            }
        }
    }
    cerr << "diag_no_five_to_four\n";

    struct RegionTask {
        vector<int> ids;
        int leaves;
        int cap;
        string why;
    };
    vector<RegionTask> tasks;
    seen_regions.clear();
    for (const Tree *guide : {&T1, &T2}) {
        for (int node = 0; node < (int)guide->nodes.size(); node++) {
            vector<int> around = local_components(*guide, F, node, 8);
            int s = (int)around.size();
            if (s < 3) continue;
            int masks = 1 << s;
            for (int mask = 1; mask < masks; mask++) {
                int k = __builtin_popcount((unsigned)mask);
                if (k < 3 || k > 6) continue;
                vector<int> ids;
                ids.reserve(k);
                for (int bit = 0; bit < s; bit++)
                    if (mask & (1 << bit)) ids.push_back(around[bit]);
                string key = region_key(ids);
                if (!seen_regions.insert(key).second) continue;

                int cap = 0;
                if (k == 3) cap = 24;
                else if (k == 4) cap = 22;
                else if (k == 5) cap = 18;
                else cap = 16;

                int leaf_count = 0;
                for (int cid : ids) leaf_count += (int)F.comp[cid].leaves.size();
                if (leaf_count > cap) continue;

                tasks.push_back({ids, leaf_count, cap,
                                 "subset_" + to_string(k) + "_to_" +
                                 to_string(k - 1)});
            }
        }
    }
    sort(tasks.begin(), tasks.end(), [](const RegionTask &a, const RegionTask &b) {
        if (a.leaves != b.leaves) return a.leaves < b.leaves;
        if (a.ids.size() != b.ids.size()) return a.ids.size() > b.ids.size();
        return a.ids < b.ids;
    });
    cerr << "diag_subset_tasks " << tasks.size() << "\n";
    int task_limit = min<int>((int)tasks.size(), 2500);
    for (int idx = 0; idx < task_limit; idx++) {
        if (idx > 0 && idx % 500 == 0)
            cerr << "diag_subset_checked " << idx << "\n";
        const RegionTask &task = tasks[idx];
        if (try_local_exact_exchange(F, T1, T2, task.ids, task.cap, before,
                                     argv[3], task.why)) {
            return 0;
        }
    }
    cerr << "diag_no_subset_exchange\n";
    return 0;
}
#endif

// ===== END heuristic/clean_refactor/h45.cpp =====
#undef main

#include "Highs.h"
#include <numeric>

struct H46LaneFiles {
    pid_t pid = -1;
    string output_path;
    string pool_path;
    string snapshot_path;
    bool running = false;
};

static string h46_read_file(const string &path) {
    ifstream input(path, ios::binary);
    stringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

static bool h46_write_file(const string &path, const string &text) {
    ofstream output(path, ios::binary);
    if (!output) return false;
    output << text;
    return (bool)output;
}

#ifndef H46_SEARCH_SECONDS_DEFAULT
#define H46_SEARCH_SECONDS_DEFAULT 250.0
#endif
#ifndef H46_ENABLE_MODE1_VALUE
#define H46_ENABLE_MODE1_VALUE 0
#endif
#ifndef H46_ENABLE_MODE2_VALUE
#define H46_ENABLE_MODE2_VALUE 1
#endif
#ifndef H46_HIGHS_SECONDS_VALUE
#define H46_HIGHS_SECONDS_VALUE 45.0
#endif

static double h46_search_seconds() {
    const char *value = getenv("H46_SEARCH_SECONDS");
    if (value == nullptr || value[0] == '\0') return H46_SEARCH_SECONDS_DEFAULT;
    char *end = nullptr;
    double seconds = strtod(value, &end);
    return end != value && seconds > 0.0 ? seconds : H46_SEARCH_SECONDS_DEFAULT;
}

static double h46_highs_seconds() {
    const char *value = getenv("H46_HIGHS_SECONDS");
    if (value == nullptr || value[0] == '\0') return H46_HIGHS_SECONDS_VALUE;
    char *end = nullptr;
    double seconds = strtod(value, &end);
    return end != value && seconds > 0.0 ? seconds : H46_HIGHS_SECONDS_VALUE;
}

static H46LaneFiles h46_start_lane(const string &input_path,
                                   const string &prefix,
                                   int lane, int score_mode) {
    H46LaneFiles files;
    files.output_path = prefix + "_lane" + to_string(lane) + ".out";
    files.pool_path = prefix + "_lane" + to_string(lane) + ".pool";
    files.snapshot_path = prefix + "_lane" + to_string(lane) + ".best";
    pid_t pid = fork();
    if (pid != 0) {
        files.pid = pid;
        files.running = pid > 0;
        return files;
    }

    ifstream lane_input(input_path);
    ofstream lane_output(files.output_path);
    if (!lane_input || !lane_output) _exit(2);

    string mode = to_string(score_mode);
    setenv("H45_SCORE_MODE", mode.c_str(), 1);
    setenv("H45_PRIMARY_SEED_OFFSET", "0", 1);
    setenv("H45_BEST_SEED_OFFSET", "193", 1);
    setenv("H45_COMPONENT_POOL", files.pool_path.c_str(), 1);
    setenv("H45_BEST_SNAPSHOT", files.snapshot_path.c_str(), 1);
    g_stop_requested = 0;
    int result = h45_solve_streams(lane_input, lane_output);
    lane_output.flush();
    cerr.flush();
    _exit(result);
}

static void h46_stop_lanes(vector<H46LaneFiles> &lanes) {
    for (H46LaneFiles &lane : lanes)
        if (lane.running) kill(lane.pid, SIGTERM);

    auto deadline = chrono::steady_clock::now() + chrono::seconds(8);
    while (chrono::steady_clock::now() < deadline) {
        bool any_running = false;
        for (H46LaneFiles &lane : lanes) {
            if (!lane.running) continue;
            int status = 0;
            pid_t result = waitpid(lane.pid, &status, WNOHANG);
            if (result == lane.pid) lane.running = false;
            else any_running = true;
        }
        if (!any_running) return;
        usleep(20000);
    }
    for (H46LaneFiles &lane : lanes) {
        if (!lane.running) continue;
        kill(lane.pid, SIGKILL);
        waitpid(lane.pid, nullptr, 0);
        lane.running = false;
    }
}

static bool h46_parse_pool(const string &path, int leaf_count,
                           unordered_map<string, vector<int>> &pool) {
    ifstream input(path);
    if (!input) return false;
    string line;
    while (getline(input, line)) {
        if (line.empty()) continue;
        vector<int> leaves;
        string token;
        stringstream parser(line);
        while (getline(parser, token, ',')) {
            if (token.empty()) return false;
            char *end = nullptr;
            long value = strtol(token.c_str(), &end, 10);
            if (end == token.c_str() || *end != '\0' ||
                value <= 0 || value > leaf_count) {
                return false;
            }
            leaves.push_back((int)value);
        }
        sort(leaves.begin(), leaves.end());
        leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());
        if (leaves.empty()) return false;
        pool.emplace(h45_component_key(leaves), std::move(leaves));
    }
    return true;
}

static bool h46_parse_forest(const string &text, int leaf_count,
                             const Tree &T1, const Tree &T2,
                             AgreementForest &F) {
    if (!parse_component_forest_text(text, leaf_count, T1, T2, F))
        return false;
    F = exact_output_repair(F, T1, T2);
    return rebuild_owners(F, T1, T2);
}

static vector<int> h46_component_edges(const Tree &T,
                                       const vector<int> &leaves) {
    vector<int> edges;
    if (leaves.size() <= 1) return edges;
    int root = lca_of_leaf_set(T, leaves);
    for (int leaf : leaves) {
        int node = T.leaf_node[leaf];
        while (node >= 0 && node != root) {
            edges.push_back(node);
            node = T.nodes[node].parent;
        }
    }
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

static bool h46_solve_pool(const unordered_map<string, vector<int>> &raw_pool,
                           int leaf_count, const Tree &T1, const Tree &T2,
                           int incumbent_components,
                           AgreementForest &solution_forest) {
    vector<vector<int>> components;
    components.reserve(raw_pool.size() + leaf_count);
    unordered_set<string> seen;
    seen.reserve(raw_pool.size() * 2 + leaf_count * 2);

    auto add_component = [&](vector<int> leaves) {
        sort(leaves.begin(), leaves.end());
        string key = h45_component_key(leaves);
        if (!seen.insert(key).second) return;
        int l1 = lca_of_leaf_set(T1, leaves);
        int l2 = lca_of_leaf_set(T2, leaves);
        if (!component_topology_valid_exact(leaves, T1, T2, l1, l2)) return;
        components.push_back(std::move(leaves));
    };
    for (const auto &entry : raw_pool) add_component(entry.second);
    for (int leaf = 1; leaf <= leaf_count; leaf++)
        add_component(vector<int>{leaf});
    cerr << "h46 pool_raw=" << raw_pool.size()
         << " pool_valid=" << components.size() << '\n';

    vector<vector<int>> by_leaf(leaf_count + 1);
    vector<vector<int>> by_first_edge(T1.nodes.size());
    vector<vector<int>> by_second_edge(T2.nodes.size());
    for (int column = 0; column < (int)components.size(); column++) {
        for (int leaf : components[column]) by_leaf[leaf].push_back(column);
        for (int edge : h46_component_edges(T1, components[column]))
            by_first_edge[edge].push_back(column);
        for (int edge : h46_component_edges(T2, components[column]))
            by_second_edge[edge].push_back(column);
    }

    string model_path = "/tmp/pace_h46_model_" +
                        to_string((long long)getpid()) + ".lp";
    ofstream model_file(model_path);
    if (!model_file) return false;
    auto write_terms = [&](const vector<int> &columns) {
        if (columns.empty()) {
            model_file << " 0\n";
            return;
        }
        int width = 1;
        model_file << ' ';
        for (int column : columns) {
            string term = "+ x" + to_string(column) + " ";
            if (width + (int)term.size() > 220) {
                model_file << "\n  ";
                width = 2;
            }
            model_file << term;
            width += (int)term.size();
        }
        model_file << '\n';
    };

    model_file << "Minimize\n obj:\n";
    vector<int> all_columns(components.size());
    iota(all_columns.begin(), all_columns.end(), 0);
    write_terms(all_columns);
    model_file << "Subject To\n";
    for (int leaf = 1; leaf <= leaf_count; leaf++) {
        model_file << " leaf_" << leaf << ":\n";
        write_terms(by_leaf[leaf]);
        model_file << " = 1\n";
    }
    for (int edge = 0; edge < (int)by_first_edge.size(); edge++) {
        if (by_first_edge[edge].size() <= 1) continue;
        model_file << " tree1_edge_" << edge << ":\n";
        write_terms(by_first_edge[edge]);
        model_file << " <= 1\n";
    }
    for (int edge = 0; edge < (int)by_second_edge.size(); edge++) {
        if (by_second_edge[edge].size() <= 1) continue;
        model_file << " tree2_edge_" << edge << ":\n";
        write_terms(by_second_edge[edge]);
        model_file << " <= 1\n";
    }
    if (incumbent_components > 1) {
        model_file << " better_than_incumbent:\n";
        write_terms(all_columns);
        model_file << " <= " << (incumbent_components - 1) << "\n";
    }
    model_file << "Binary\n";
    for (int column = 0; column < (int)components.size(); column++)
        model_file << " x" << column << '\n';
    model_file << "End\n";
    model_file.close();

    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("time_limit", h46_highs_seconds());
    highs.setOptionValue("threads", 1);
    highs.setOptionValue("mip_rel_gap", 0.0);
    highs.setOptionValue("random_seed", 1);
    HighsStatus pass_status = highs.readModel(model_path);
    unlink(model_path.c_str());
    if (pass_status == HighsStatus::kError) {
        cerr << "h46 highs_pass_error\n";
        return false;
    }
    const char *write_model = getenv("H46_WRITE_MODEL");
    if (write_model && write_model[0]) highs.writeModel(write_model);
    HighsStatus run_status = highs.run();
    if (run_status == HighsStatus::kError) {
        cerr << "h46 highs_run_error\n";
        return false;
    }
    cerr << "h46 highs_status="
         << highs.modelStatusToString(highs.getModelStatus())
         << " objective=" << highs.getInfo().objective_function_value
         << '\n';
    const HighsSolution &solution = highs.getSolution();
    if (!solution.value_valid ||
        solution.col_value.size() != components.size()) {
        cerr << "h46 highs_no_values value_valid=" << solution.value_valid
             << " columns=" << solution.col_value.size() << '\n';
        highs.resetGlobalScheduler(true);
        return false;
    }

    AgreementForest candidate;
    candidate.leaf_comp.assign(leaf_count + 1, -1);
    int selected = 0;
    vector<int> selected_leaf_count(leaf_count + 1, 0);
    for (size_t column = 0; column < components.size(); column++) {
        if (solution.col_value[column] <= 0.5) continue;
        string column_name;
        if (highs.getColName((HighsInt)column, column_name) !=
                HighsStatus::kOk ||
            column_name.size() < 2 || column_name[0] != 'x') {
            return false;
        }
        size_t original_column = stoull(column_name.substr(1));
        if (original_column >= components.size()) return false;
        append_component(candidate, components[original_column],
                         components[original_column].front());
        for (int leaf : components[original_column])
            selected_leaf_count[leaf]++;
        selected++;
    }
    int missing_leaves = 0;
    int duplicate_leaves = 0;
    for (int leaf = 1; leaf <= leaf_count; leaf++) {
        missing_leaves += selected_leaf_count[leaf] == 0;
        duplicate_leaves += selected_leaf_count[leaf] > 1;
    }
    if (!solution.row_value.empty()) {
        for (int leaf = 1; leaf <= leaf_count; leaf++) {
            if (fabs(solution.row_value[leaf - 1] -
                     selected_leaf_count[leaf]) > 1e-6) {
                cerr << "h46 first_row_mismatch leaf=" << leaf
                     << " solver=" << solution.row_value[leaf - 1]
                     << " recomputed=" << selected_leaf_count[leaf] << '\n';
                break;
            }
        }
    }
    cerr << "h46 highs_selected=" << selected
         << " incumbent=" << incumbent_components
         << " missing_leaves=" << missing_leaves
         << " duplicate_leaves=" << duplicate_leaves << '\n';
    if (selected <= 0 || selected >= incumbent_components) {
        highs.resetGlobalScheduler(true);
        return false;
    }
    const char *debug_selected = getenv("H46_DEBUG_SELECTED");
    if (debug_selected && debug_selected[0])
        h46_write_file(debug_selected, forest_to_output(candidate, T1));
    if (!rebuild_owners(candidate, T1, T2)) {
        cerr << "h46 highs_selected_rebuild_failed\n";
        highs.resetGlobalScheduler(true);
        return false;
    }
    candidate = exact_output_repair(candidate, T1, T2);
    if (active_count(candidate) != selected ||
        !rebuild_owners(candidate, T1, T2)) {
        cerr << "h46 highs_selected_exact_repair_failed after="
             << active_count(candidate) << '\n';
        highs.resetGlobalScheduler(true);
        return false;
    }
    solution_forest = std::move(candidate);
    cerr << "h46 highs_improvement " << incumbent_components
         << " -> " << selected
         << " candidates=" << components.size() << '\n';
    highs.resetGlobalScheduler(true);
    return true;
}

#ifndef H46_NO_MAIN
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    stringstream input_buffer;
    input_buffer << cin.rdbuf();
    string input_text = input_buffer.str();
    if (input_text.empty()) return 1;

    TwoTreeInstance instance;
    stringstream instance_input(input_text);
    if (!read_two_tree_instance(instance_input, instance)) return 1;
    if (instance.leaf_count > 500) {
        while (true) pause();
    }
    Tree T1 = parse_newick(instance.trees[0], instance.leaf_count);
    Tree T2 = parse_newick(instance.trees[1], instance.leaf_count);
    if (T1.root < 0 || T2.root < 0) return 1;
    build_lca(T1);
    build_lca(T2);

    string prefix = "/tmp/pace_h46_" + to_string((long long)getpid());
    string input_path = prefix + ".in";
    if (!h46_write_file(input_path, input_text)) return 1;

    vector<H46LaneFiles> lanes;
    const char *existing_out0 = getenv("H46_EXISTING_OUT0");
    const char *existing_pool0 = getenv("H46_EXISTING_POOL0");
    const char *existing_out2 = getenv("H46_EXISTING_OUT2");
    const char *existing_pool2 = getenv("H46_EXISTING_POOL2");
    bool use_existing = existing_out0 && existing_pool0 &&
                        existing_out2 && existing_pool2;
    if (use_existing) {
        H46LaneFiles lane0;
        lane0.output_path = existing_out0;
        lane0.pool_path = existing_pool0;
        H46LaneFiles lane2;
        lane2.output_path = existing_out2;
        lane2.pool_path = existing_pool2;
        lanes.push_back(std::move(lane0));
        lanes.push_back(std::move(lane2));
    } else {
        lanes.push_back(h46_start_lane(input_path, prefix, 0, 0));
        int next_lane = 1;
        if (H46_ENABLE_MODE1_VALUE) {
            lanes.push_back(h46_start_lane(input_path, prefix, next_lane++, 1));
        }
        if (H46_ENABLE_MODE2_VALUE) {
            lanes.push_back(h46_start_lane(input_path, prefix, next_lane++, 2));
        }
        bool any_lane_running = false;
        for (const H46LaneFiles &lane : lanes)
            any_lane_running = any_lane_running || lane.running;
        if (!any_lane_running) return 1;

        signal(SIGTERM, handle_signal);
        signal(SIGINT, handle_signal);
        auto deadline = chrono::steady_clock::now() + chrono::milliseconds(
            (long long)(1000.0 * h46_search_seconds()));
        while (!g_stop_requested && chrono::steady_clock::now() < deadline)
            usleep(50000);
        h46_stop_lanes(lanes);
    }

    AgreementForest best;
    bool have_best = false;
    unordered_map<string, vector<int>> component_pool;
    for (const H46LaneFiles &lane : lanes) {
        vector<string> candidate_paths;
        if (!lane.output_path.empty()) candidate_paths.push_back(lane.output_path);
        if (!lane.snapshot_path.empty()) candidate_paths.push_back(lane.snapshot_path);
        for (const string &candidate_path : candidate_paths) {
            AgreementForest candidate;
            string output = h46_read_file(candidate_path);
            if (h46_parse_forest(output, instance.leaf_count, T1, T2,
                                 candidate)) {
                int candidate_components = active_count(candidate);
                if (!have_best || candidate_components < active_count(best)) {
                    best = std::move(candidate);
                    have_best = true;
                }
            }
        }
        h46_parse_pool(lane.pool_path, instance.leaf_count, component_pool);
    }
    if (!have_best) return 1;
    for (const Component &component : best.comp)
        if (component.active)
            component_pool.emplace(h45_component_key(component.leaves),
                                   component.leaves);

    AgreementForest optimized;
    if (!g_stop_requested &&
        h46_solve_pool(component_pool, instance.leaf_count, T1, T2,
                       active_count(best), optimized)) {
        best = std::move(optimized);
    }

    cout << forest_to_output(best, T1);
    cout.flush();
    if (!use_existing) {
        for (const H46LaneFiles &lane : lanes) {
            unlink(lane.output_path.c_str());
            unlink(lane.pool_path.c_str());
            unlink(lane.snapshot_path.c_str());
        }
    }
    unlink(input_path.c_str());
    return 0;
}
#endif

// ===== END heuristic/clean_refactor/h46.cpp =====

#include <atomic>
#include <thread>

#ifndef H56_ROUTE_SECONDS_VALUE
#define H56_ROUTE_SECONDS_VALUE 165.0
#endif

#ifndef H56_TOTAL_SECONDS_VALUE
#define H56_TOTAL_SECONDS_VALUE 285.0
#endif

#ifndef H56_OUTPUT_SAFETY_SECONDS_VALUE
#define H56_OUTPUT_SAFETY_SECONDS_VALUE 8.0
#endif

#ifndef H56_HIGHS_SETUP_SECONDS_VALUE
#define H56_HIGHS_SETUP_SECONDS_VALUE 10.0
#endif

#ifndef H56_MIN_HIGHS_SECONDS_VALUE
#define H56_MIN_HIGHS_SECONDS_VALUE 8.0
#endif

#ifndef H56_MAX_POOL_VALUE
#define H56_MAX_POOL_VALUE 7900
#endif

#ifndef H56_MAX_COMPONENT_SIZE_VALUE
#define H56_MAX_COMPONENT_SIZE_VALUE 16
#endif

#ifndef H56_LOCAL_POOL_LIMIT_VALUE
#define H56_LOCAL_POOL_LIMIT_VALUE 4000
#endif

#ifndef H56_ENABLE_EXTERNAL_POOL_VALUE
#define H56_ENABLE_EXTERNAL_POOL_VALUE 0
#endif

#ifndef H56_EXTERNAL_POOL_LIMIT_VALUE
#define H56_EXTERNAL_POOL_LIMIT_VALUE H56_MAX_POOL_VALUE
#endif

#ifndef H56_ILP_PASSES_VALUE
#define H56_ILP_PASSES_VALUE 1
#endif

#ifndef H56_REENRICH_SECONDS_VALUE
#define H56_REENRICH_SECONDS_VALUE 12.0
#endif

#ifndef H56_SECOND_MAX_POOL_VALUE
#define H56_SECOND_MAX_POOL_VALUE 0
#endif

#ifndef H56_SECOND_EXTERNAL_POOL_LIMIT_VALUE
#define H56_SECOND_EXTERNAL_POOL_LIMIT_VALUE 0
#endif

#ifndef H56_ENABLE_COMPLETION_VALUE
#define H56_ENABLE_COMPLETION_VALUE 0
#endif

#ifndef H56_ENABLE_BATCH_COMPLETION_VALUE
#define H56_ENABLE_BATCH_COMPLETION_VALUE 0
#endif

struct H56RouteFiles {
    string output_path;
    string pool_path;
    string snapshot_path;
    string previous_snapshot_path;
    string milestone_prefix;
    vector<string> milestone_paths;
};

struct H56Region {
    vector<int> leaves;
    int score = 0;
};

static int h56_env_int(const char *name, int fallback) {
    const char *value = getenv(name);
    if (value == nullptr || value[0] == '\0') return fallback;
    char *end = nullptr;
    long parsed = strtol(value, &end, 10);
    if (end == value || *end != '\0' || parsed < INT_MIN ||
        parsed > INT_MAX) {
        return fallback;
    }
    return (int)parsed;
}

static double h56_env_double(const char *name, double fallback) {
    const char *value = getenv(name);
    if (value == nullptr || value[0] == '\0') return fallback;
    char *end = nullptr;
    double parsed = strtod(value, &end);
    return end != value && parsed > 0.0 ? parsed : fallback;
}

static vector<int> h56_env_modes(const char *name,
                                 const vector<int> &fallback) {
    const char *value = getenv(name);
    if (value == nullptr || value[0] == '\0') return fallback;
    vector<int> modes;
    string text(value);
    size_t start = 0;
    while (start <= text.size()) {
        size_t comma = text.find(',', start);
        string token = text.substr(start, comma == string::npos
                                           ? string::npos
                                           : comma - start);
        if (!token.empty()) {
            char *end = nullptr;
            long parsed = strtol(token.c_str(), &end, 10);
            if (end != token.c_str()) modes.push_back((int)parsed);
        }
        if (comma == string::npos) break;
        start = comma + 1;
    }
    if (modes.empty()) return fallback;
    sort(modes.begin(), modes.end());
    modes.erase(unique(modes.begin(), modes.end()), modes.end());
    return modes;
}

static double h56_elapsed_seconds(
    const chrono::steady_clock::time_point &start) {
    return chrono::duration<double>(chrono::steady_clock::now() - start)
        .count();
}

static string h56_key(vector<int> leaves) {
    sort(leaves.begin(), leaves.end());
    leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());
    return h45_component_key(leaves);
}

static vector<int> h56_union_component_leaves(const AgreementForest &F,
                                              const vector<int> &ids) {
    vector<int> leaves;
    for (int cid : ids) {
        if (cid < 0 || cid >= (int)F.comp.size() || !F.comp[cid].active)
            continue;
        leaves.insert(leaves.end(), F.comp[cid].leaves.begin(),
                      F.comp[cid].leaves.end());
    }
    sort(leaves.begin(), leaves.end());
    leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());
    return leaves;
}

static int h56_leaf_distance(const Tree &T, int a, int b) {
    int na = T.leaf_node[a];
    int nb = T.leaf_node[b];
    int root = lca_node(T, na, nb);
    if (root < 0) return 1 << 28;
    return T.nodes[na].depth + T.nodes[nb].depth -
           2 * T.nodes[root].depth;
}

static int h56_component_edge_count(const Tree &T,
                                    const vector<int> &leaves) {
    if (leaves.size() <= 1) return 0;
    int root = lca_of_leaf_set(T, leaves);
    int count = 0;
    for (int leaf : leaves) {
        int node = T.leaf_node[leaf];
        while (node >= 0 && node != root) {
            count++;
            node = T.nodes[node].parent;
        }
    }
    return count;
}

static bool h56_component_valid(const vector<int> &leaves,
                                const Tree &T1, const Tree &T2) {
    if (leaves.empty()) return false;
    if (leaves.size() == 1) return true;
    int root1 = lca_of_leaf_set(T1, leaves);
    int root2 = lca_of_leaf_set(T2, leaves);
    return root1 >= 0 && root2 >= 0 &&
           component_topology_valid_exact(leaves, T1, T2, root1, root2);
}

static int h56_component_score(const vector<int> &leaves,
                               const Tree &T1, const Tree &T2) {
    int root1 = lca_of_leaf_set(T1, leaves);
    int root2 = lca_of_leaf_set(T2, leaves);
    int depth = (root1 >= 0 ? T1.nodes[root1].depth : 0) +
                (root2 >= 0 ? T2.nodes[root2].depth : 0);
    int edges = h56_component_edge_count(T1, leaves) +
                h56_component_edge_count(T2, leaves);
    int size = (int)leaves.size();
    return 100000 * (size - 1) + 1200 * depth - 80 * edges;
}

static bool h56_add_candidate(unordered_map<string, vector<int>> &pool,
                              vector<int> leaves,
                              const Tree &T1, const Tree &T2,
                              int max_pool) {
    if ((int)pool.size() >= max_pool) return false;
    sort(leaves.begin(), leaves.end());
    leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());
    if (leaves.empty()) return false;
    string key = h45_component_key(leaves);
    if (pool.count(key)) return false;
    if (!h56_component_valid(leaves, T1, T2)) return false;
    pool.emplace(key, std::move(leaves));
    return true;
}

static int h56_node_distance(const Tree &T, int a, int b) {
    int root = lca_node(T, a, b);
    if (root < 0) return 1 << 28;
    return T.nodes[a].depth + T.nodes[b].depth -
           2 * T.nodes[root].depth;
}

static bool h56_components_disjoint(const vector<int> &a,
                                    const vector<int> &b) {
    int i = 0;
    int j = 0;
    while (i < (int)a.size() && j < (int)b.size()) {
        if (a[i] == b[j]) return false;
        if (a[i] < b[j]) i++;
        else j++;
    }
    return true;
}

static vector<vector<int>>
h56_valid_pool_components(const unordered_map<string, vector<int>> &pool,
                          const Tree &T1, const Tree &T2) {
    vector<vector<int>> components;
    components.reserve(pool.size());
    for (const auto &entry : pool) {
        vector<int> leaves = entry.second;
        sort(leaves.begin(), leaves.end());
        leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());
        if (!h56_component_valid(leaves, T1, T2)) continue;
        components.push_back(std::move(leaves));
    }
    sort(components.begin(), components.end(),
         [](const vector<int> &a, const vector<int> &b) {
             if (a.size() != b.size()) return a.size() < b.size();
             return a < b;
         });
    components.erase(unique(components.begin(), components.end()),
                     components.end());
    return components;
}

static int h56_add_forest_components(
    unordered_map<string, vector<int>> &pool,
    const AgreementForest &F) {
    int added = 0;
    for (const Component &component : F.comp) {
        if (!component.active || component.leaves.empty()) continue;
        string key = h45_component_key(component.leaves);
        if (pool.emplace(key, component.leaves).second) added++;
    }
    return added;
}

static void h56_add_fast_greedy_pool(
    unordered_map<string, vector<int>> &pool,
    AgreementForest &best,
    const Tree &T1, const Tree &T2) {
    if (!h56_env_int("H56_ENABLE_FAST_GREEDY_POOL", 1)) return;
    vector<int> modes =
        h56_env_modes("H56_FAST_GREEDY_MODES", vector<int>{0, 1, 2, 12});
    auto start = chrono::steady_clock::now();
    for (int mode : modes) {
        if (g_stop_requested) break;
        SearchState search;
        search.merge_score_mode = mode;
        AgreementForest raw = greedy_merge(
            singleton_forest(T1.n), T1, T2, search, false, false);
        AgreementForest candidate = exact_output_repair(raw, T1, T2);
        int added = h56_add_forest_components(pool, candidate);
        int k = active_count(candidate);
        if (k < active_count(best)) best = candidate;
        double elapsed = h56_elapsed_seconds(start);
        cerr << "h56 fast_greedy mode=" << mode
             << " components=" << k
             << " added=" << added
             << " elapsed=" << elapsed << '\n';
    }
}

static bool h56_run_h45_route(const string &input_text,
                              const string &prefix,
                              int score_mode,
                              int primary_seed,
                              int best_seed,
                              double route_seconds,
                              H56RouteFiles &files) {
    files.output_path = prefix + ".out";
    files.pool_path = prefix + ".pool";
    files.snapshot_path = prefix + ".best";
    files.previous_snapshot_path = prefix + ".prev";
    files.milestone_prefix = prefix + ".milestone";
    files.milestone_paths.clear();
    for (int threshold : {430, 428, 427, 426}) {
        files.milestone_paths.push_back(
            files.milestone_prefix + "_" + to_string(threshold) + ".out");
    }

    setenv("H45_SCORE_MODE", to_string(score_mode).c_str(), 1);
    setenv("H45_PRIMARY_SEED_OFFSET", to_string(primary_seed).c_str(), 1);
    setenv("H45_BEST_SEED_OFFSET", to_string(best_seed).c_str(), 1);
    setenv("H45_COMPONENT_POOL", files.pool_path.c_str(), 1);
    setenv("H45_BEST_SNAPSHOT", files.snapshot_path.c_str(), 1);
    setenv("H45_PREVIOUS_BEST_SNAPSHOT",
           files.previous_snapshot_path.c_str(), 1);
    setenv("H45_MILESTONE_PREFIX", files.milestone_prefix.c_str(), 1);
    setenv("H45_MILESTONE_THRESHOLDS", "430,428,427,426", 1);

    g_stop_requested = 0;
    g_h45_component_pool.clear();

    stringstream input(input_text);
    ofstream output(files.output_path, ios::binary);
    if (!output) return false;

    atomic<bool> route_done(false);
    thread timer([&]() {
        auto deadline = chrono::steady_clock::now() +
                        chrono::milliseconds((long long)(route_seconds * 1000.0));
        while (!route_done.load(memory_order_relaxed) &&
               chrono::steady_clock::now() < deadline) {
            this_thread::sleep_for(chrono::milliseconds(20));
        }
        if (!route_done.load(memory_order_relaxed)) {
            g_stop_requested = 1;
            cerr << "h56 route_timer_stop seconds=" << route_seconds
                 << " score_mode=" << score_mode << '\n';
        }
    });

    int rc = h45_solve_streams(input, output);
    route_done.store(true, memory_order_relaxed);
    timer.join();
    g_stop_requested = 0;
    output.flush();
    cerr << "h56 route_done rc=" << rc
         << " score_mode=" << score_mode
         << " primary_seed=" << primary_seed
         << " best_seed=" << best_seed << '\n';
    return rc == 0;
}

static vector<H56Region> h56_collect_regions(const AgreementForest &best,
                                             const Tree &T1,
                                             const Tree &T2) {
    const int region_limit = h56_env_int("H56_REGION_LIMIT", 90);
    const int component_limit = h56_env_int("H56_REGION_COMPONENT_LIMIT", 18);
    const int min_leaves = h56_env_int("H56_REGION_MIN_LEAVES", 6);
    const int max_leaves = h56_env_int("H56_REGION_MAX_LEAVES", 42);
    unordered_map<string, H56Region> by_key;

    auto consider = [&](vector<int> ids, const Tree &guide, int node,
                        int conflict_bonus = 0) {
        if ((int)ids.size() < 2 || (int)ids.size() > component_limit)
            return;
        vector<int> leaves = h56_union_component_leaves(best, ids);
        if ((int)leaves.size() < min_leaves ||
            (int)leaves.size() > max_leaves) {
            return;
        }
        string key = h45_component_key(leaves);
        int score = 4096 * (int)ids.size() - 17 * (int)leaves.size() +
                    (node >= 0 ? guide.nodes[node].depth : 0) +
                    conflict_bonus;
        auto it = by_key.find(key);
        if (it == by_key.end() || score > it->second.score) {
            H56Region region;
            region.leaves = std::move(leaves);
            region.score = score;
            by_key[key] = std::move(region);
        }
    };

    const Tree *guides[2] = {&T1, &T2};
    for (const Tree *guide : guides) {
        for (int node = 0; node < (int)guide->nodes.size(); node++) {
            vector<int> ids = local_components(*guide, best, node,
                                               component_limit + 1);
            consider(std::move(ids), *guide, node);
        }
    }

    // Also inspect a few ancestors of every active component root. This is
    // where the 424-style medium components tend to sit: not a whole local
    // clade, but a compact neighborhood around one troublesome component.
    const int ancestor_steps = h56_env_int("H56_REGION_ANCESTOR_STEPS", 4);
    for (const Component &component : best.comp) {
        if (!component.active || component.leaves.empty()) continue;
        const Tree *trees[2] = {&T1, &T2};
        for (const Tree *tree : trees) {
            int node = lca_of_leaf_set(*tree, component.leaves);
            for (int step = 0;
                 node >= 0 && step <= ancestor_steps;
                 step++, node = tree->nodes[node].parent) {
                vector<int> ids = local_components(*tree, best, node,
                                                   component_limit + 1);
                consider(std::move(ids), *tree, node);
            }
        }
    }

    if (h56_env_int("H56_ENABLE_CONFLICT_REGION_RANK", 0)) {
        // Conflict-ranked regions: if two active components are a cherry in one
        // tree but not the other, the useful replacement usually lives around
        // their LCA in the opposite tree. Keep this optional because changing
        // the region order can displace the known-good 424 pool on 000.
        auto add_conflict_regions =
            [&](const Tree &other_tree,
                const vector<pair<int, int>> &cherries,
                const unordered_set<uint64_t> &other_set) {
            int used = 0;
            const int limit =
                h56_env_int("H56_CONFLICT_REGION_PAIR_LIMIT", 500);
            for (auto [A, B] : cherries) {
                if (used >= limit) break;
                if (other_set.count(comp_pair_key(A, B))) continue;
                if (A < 0 || B < 0 || A >= (int)best.comp.size() ||
                    B >= (int)best.comp.size() || !best.comp[A].active ||
                    !best.comp[B].active) {
                    continue;
                }
                int ra = lca_of_leaf_set(other_tree, best.comp[A].leaves);
                int rb = lca_of_leaf_set(other_tree, best.comp[B].leaves);
                int node = lca_node(other_tree, ra, rb);
                if (node < 0) continue;
                vector<int> ids =
                    local_components(other_tree, best, node,
                                     component_limit + 1);
                int size_bonus =
                    min(2048, 64 * ((int)best.comp[A].leaves.size() +
                                    (int)best.comp[B].leaves.size()));
                consider(std::move(ids), other_tree, node,
                         20000 + size_bonus);
                used++;
            }
        };

        vector<pair<int, int>> cherries1 = active_cherry_pairs(T1, best);
        vector<pair<int, int>> cherries2 = active_cherry_pairs(T2, best);
        unordered_set<uint64_t> cherry_set1, cherry_set2;
        cherry_set1.reserve(cherries1.size() * 2 + 1);
        cherry_set2.reserve(cherries2.size() * 2 + 1);
        for (auto [A, B] : cherries1) cherry_set1.insert(comp_pair_key(A, B));
        for (auto [A, B] : cherries2) cherry_set2.insert(comp_pair_key(A, B));
        add_conflict_regions(T2, cherries1, cherry_set2);
        add_conflict_regions(T1, cherries2, cherry_set1);
    }

    vector<H56Region> regions;
    regions.reserve(by_key.size());
    for (auto &entry : by_key) regions.push_back(std::move(entry.second));
    sort(regions.begin(), regions.end(),
         [](const H56Region &a, const H56Region &b) {
             if (a.score != b.score) return a.score > b.score;
             if (a.leaves.size() != b.leaves.size())
                 return a.leaves.size() > b.leaves.size();
             return a.leaves < b.leaves;
         });
    if ((int)regions.size() > region_limit) regions.resize(region_limit);
    cerr << "h56 regions=" << regions.size()
         << " raw=" << by_key.size() << '\n';
    return regions;
}

static void h56_generate_region_components(
    const H56Region &region,
    unordered_map<string, vector<int>> &pool,
    const Tree &T1, const Tree &T2,
    int max_pool) {
    const int max_component_size = h56_env_int("H56_MAX_COMPONENT_SIZE", 16);
    const int pair_distance = h56_env_int("H56_PAIR_DISTANCE", 24);
    const int add_distance = h56_env_int("H56_ADD_DISTANCE", 24);
    const int per_region_limit = h56_env_int("H56_PER_REGION_LIMIT", 140);
    const int frontier_limit = h56_env_int("H56_FRONTIER_LIMIT", 48);
    const vector<int> &region_leaves = region.leaves;
    unordered_set<string> local_seen;
    vector<vector<int>> frontier;
    int added = 0;

    auto try_add = [&](vector<int> leaves) {
        sort(leaves.begin(), leaves.end());
        leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());
        string key = h45_component_key(leaves);
        if (!local_seen.insert(key).second) return false;
        size_t before = pool.size();
        bool ok = h56_add_candidate(pool, std::move(leaves), T1, T2, max_pool);
        if (ok && pool.size() != before) added++;
        return ok;
    };

    for (int i = 0; i < (int)region_leaves.size(); i++) {
        for (int j = i + 1; j < (int)region_leaves.size(); j++) {
            int a = region_leaves[i];
            int b = region_leaves[j];
            int distance = max(h56_leaf_distance(T1, a, b),
                               h56_leaf_distance(T2, a, b));
            if (distance > pair_distance) continue;
            vector<int> pair{a, b};
            if (try_add(pair)) frontier.push_back(std::move(pair));
            if (added >= per_region_limit || (int)pool.size() >= max_pool)
                return;
        }
    }

    sort(frontier.begin(), frontier.end(),
         [&](const vector<int> &a, const vector<int> &b) {
             int sa = h56_component_score(a, T1, T2);
             int sb = h56_component_score(b, T1, T2);
             if (sa != sb) return sa > sb;
             return a < b;
         });
    if ((int)frontier.size() > frontier_limit)
        frontier.resize(frontier_limit);

    for (int target_size = 3;
         target_size <= max_component_size && !frontier.empty() &&
         added < per_region_limit && (int)pool.size() < max_pool;
         target_size++) {
        vector<vector<int>> next_frontier;
        unordered_set<string> next_seen;

        for (const vector<int> &state : frontier) {
            vector<char> in_region(T1.n + 1, 0);
            for (int leaf : state) in_region[leaf] = 1;

            vector<pair<int, int>> options;
            for (int leaf : region_leaves) {
                if (in_region[leaf]) continue;
                int nearest = 1 << 28;
                for (int existing : state) {
                    nearest = min(nearest, max(
                        h56_leaf_distance(T1, existing, leaf),
                        h56_leaf_distance(T2, existing, leaf)));
                }
                if (nearest <= add_distance) options.push_back({nearest, leaf});
            }
            sort(options.begin(), options.end(),
                 [](const pair<int, int> &a, const pair<int, int> &b) {
                     if (a.first != b.first) return a.first < b.first;
                     return a.second < b.second;
                 });
            if ((int)options.size() > 24) options.resize(24);

            for (const auto &option : options) {
                vector<int> grown = state;
                grown.push_back(option.second);
                sort(grown.begin(), grown.end());
                grown.erase(unique(grown.begin(), grown.end()), grown.end());
                if ((int)grown.size() != target_size) continue;
                string key = h45_component_key(grown);
                if (!next_seen.insert(key).second) continue;
                if (try_add(grown)) next_frontier.push_back(std::move(grown));
                if (added >= per_region_limit || (int)pool.size() >= max_pool)
                    return;
            }
        }

        sort(next_frontier.begin(), next_frontier.end(),
             [&](const vector<int> &a, const vector<int> &b) {
                 int sa = h56_component_score(a, T1, T2);
                 int sb = h56_component_score(b, T1, T2);
                 if (sa != sb) return sa > sb;
                 if (a.size() != b.size()) return a.size() > b.size();
                 return a < b;
             });
        next_frontier.erase(unique(next_frontier.begin(),
                                   next_frontier.end()),
                            next_frontier.end());
        if ((int)next_frontier.size() > frontier_limit)
            next_frontier.resize(frontier_limit);
        frontier = std::move(next_frontier);
    }
}

static vector<vector<int>> h56_region_restricted_node_sets(
    const Tree &T,
    const vector<int> &region_leaves) {
    const int max_node_set = h56_env_int("H56_NODE_SET_MAX_LEAVES", 42);
    const int keep_limit = h56_env_int("H56_NODE_SET_KEEP", 120);
    vector<char> mark(T.n + 1, 0);
    for (int leaf : region_leaves) mark[leaf] = 1;

    vector<int> order;
    order.reserve(T.nodes.size());
    vector<int> stack;
    if (T.root >= 0) stack.push_back(T.root);
    while (!stack.empty()) {
        int node = stack.back();
        stack.pop_back();
        order.push_back(node);
        for (int child : T.nodes[node].child) stack.push_back(child);
    }

    vector<vector<int>> below(T.nodes.size());
    vector<vector<int>> sets;
    for (int oi = (int)order.size() - 1; oi >= 0; oi--) {
        int node = order[oi];
        const TreeNode &tn = T.nodes[node];
        vector<int> leaves;
        if (tn.label >= 1) {
            if (tn.label <= T.n && mark[tn.label]) leaves.push_back(tn.label);
        } else {
            for (int child : tn.child) {
                leaves.insert(leaves.end(),
                              below[child].begin(), below[child].end());
                if ((int)leaves.size() > max_node_set) break;
            }
            sort(leaves.begin(), leaves.end());
            leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());
        }
        if ((int)leaves.size() <= max_node_set)
            below[node] = leaves;
        if ((int)leaves.size() >= 2 &&
            (int)leaves.size() <= max_node_set) {
            sets.push_back(std::move(leaves));
        }
    }

    sort(sets.begin(), sets.end(),
         [&](const vector<int> &a, const vector<int> &b) {
             if (a.size() != b.size()) return a.size() > b.size();
             int sa = h56_component_score(a, T, T);
             int sb = h56_component_score(b, T, T);
             if (sa != sb) return sa > sb;
             return a < b;
         });
    sets.erase(unique(sets.begin(), sets.end()), sets.end());
    if ((int)sets.size() > keep_limit) sets.resize(keep_limit);
    return sets;
}

static vector<int> h56_intersection(const vector<int> &a,
                                    const vector<int> &b) {
    vector<int> out;
    out.reserve(min(a.size(), b.size()));
    int i = 0;
    int j = 0;
    while (i < (int)a.size() && j < (int)b.size()) {
        if (a[i] == b[j]) {
            out.push_back(a[i]);
            i++;
            j++;
        } else if (a[i] < b[j]) {
            i++;
        } else {
            j++;
        }
    }
    return out;
}

static void h56_generate_clade_intersection_components(
    const H56Region &region,
    unordered_map<string, vector<int>> &pool,
    const Tree &T1, const Tree &T2,
    int max_pool) {
    const int medium_min =
        h56_env_int("H56_MEDIUM_INSTANCE_MIN_LEAVES",
                    H56_MEDIUM_INSTANCE_MIN_LEAVES_VALUE);
    const int default_clade_cap =
        T1.n >= medium_min
            ? h56_env_int("H56_REGION_MAX_LEAVES", 42)
            : H56_MAX_COMPONENT_SIZE_VALUE;
    const int max_component_size =
        h56_env_int("H56_CLADE_MAX_COMPONENT_SIZE", default_clade_cap);
    const int per_region_limit = h56_env_int("H56_CLADE_PER_REGION_LIMIT", 80);
    int added = 0;

    vector<vector<int>> first_sets =
        h56_region_restricted_node_sets(T1, region.leaves);
    vector<vector<int>> second_sets =
        h56_region_restricted_node_sets(T2, region.leaves);

    auto try_add = [&](vector<int> leaves) {
        if ((int)leaves.size() < 2 ||
            (int)leaves.size() > max_component_size) {
            return;
        }
        size_t before = pool.size();
        h56_add_candidate(pool, std::move(leaves), T1, T2, max_pool);
        if (pool.size() != before) added++;
    };

    for (const vector<int> &leaves : first_sets) {
        try_add(leaves);
        if (added >= per_region_limit || (int)pool.size() >= max_pool) return;
    }
    for (const vector<int> &leaves : second_sets) {
        try_add(leaves);
        if (added >= per_region_limit || (int)pool.size() >= max_pool) return;
    }

    for (const vector<int> &a : first_sets) {
        for (const vector<int> &b : second_sets) {
            vector<int> common = h56_intersection(a, b);
            try_add(std::move(common));
            if (added >= per_region_limit || (int)pool.size() >= max_pool)
                return;
        }
    }
}

static vector<int> h56_difference(const vector<int> &a,
                                  const vector<int> &b) {
    vector<int> out;
    out.reserve(a.size());
    int i = 0;
    int j = 0;
    while (i < (int)a.size()) {
        while (j < (int)b.size() && b[j] < a[i]) j++;
        if (j >= (int)b.size() || a[i] != b[j]) out.push_back(a[i]);
        i++;
    }
    return out;
}

static bool h56_is_subset_sorted(const vector<int> &small,
                                 const vector<int> &large) {
    int i = 0;
    int j = 0;
    while (i < (int)small.size() && j < (int)large.size()) {
        if (small[i] == large[j]) {
            i++;
            j++;
        } else if (small[i] > large[j]) {
            j++;
        } else {
            return false;
        }
    }
    return i == (int)small.size();
}

static void h56_add_pruned_clade_components(
    unordered_map<string, vector<int>> &pool,
    const vector<H56Region> &regions,
    const Tree &T1,
    const Tree &T2,
    int max_pool,
    const chrono::steady_clock::time_point &global_start,
    double pool_cutoff_seconds) {
    if (!h56_env_int("H56_ENABLE_PRUNED_CLADES", 0)) return;

    const int medium_min =
        h56_env_int("H56_MEDIUM_INSTANCE_MIN_LEAVES",
                    H56_MEDIUM_INSTANCE_MIN_LEAVES_VALUE);
    const bool medium_instance = T1.n >= medium_min;
    const int region_limit =
        h56_env_int("H56_PRUNED_REGION_LIMIT",
                    medium_instance ? 60 : 36);
    const int base_limit =
        h56_env_int("H56_PRUNED_BASE_LIMIT",
                    medium_instance ? 160 : 90);
    const int remove_limit =
        h56_env_int("H56_PRUNED_REMOVE_LIMIT",
                    medium_instance ? 80 : 50);
    const int max_component_size =
        h56_env_int("H56_PRUNED_MAX_COMPONENT_SIZE",
                    medium_instance ? 42 : 18);
    const int min_component_size =
        h56_env_int("H56_PRUNED_MIN_COMPONENT_SIZE", 2);
    const int max_remove_size =
        h56_env_int("H56_PRUNED_MAX_REMOVE_SIZE",
                    medium_instance ? 18 : 8);
    const int add_limit =
        h56_env_int("H56_PRUNED_ADD_LIMIT",
                    medium_instance ? 1800 : 900);
    const double seconds_limit =
        h56_env_double("H56_PRUNED_SECONDS",
                       medium_instance ? 8.0 : 4.0);

    int before = (int)pool.size();
    int accepted = 0;
    int used_regions = 0;
    long long tried = 0;
    auto start = chrono::steady_clock::now();

    auto still_time = [&]() {
        double local_elapsed = chrono::duration<double>(
            chrono::steady_clock::now() - start).count();
        return local_elapsed < seconds_limit &&
               h56_elapsed_seconds(global_start) < pool_cutoff_seconds &&
               !g_stop_requested;
    };

    auto normalize = [](vector<int> leaves) {
        sort(leaves.begin(), leaves.end());
        leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());
        return leaves;
    };

    auto try_add = [&](vector<int> leaves) {
        leaves = normalize(std::move(leaves));
        if ((int)leaves.size() < min_component_size ||
            (int)leaves.size() > max_component_size) {
            return;
        }
        tried++;
        size_t old_size = pool.size();
        if (h56_add_candidate(pool, std::move(leaves), T1, T2, max_pool))
            accepted += pool.size() != old_size;
    };

    for (const H56Region &region : regions) {
        if ((int)pool.size() >= max_pool || accepted >= add_limit ||
            !still_time()) {
            break;
        }
        used_regions++;

        vector<vector<int>> bases =
            h56_region_restricted_node_sets(T1, region.leaves);
        vector<vector<int>> second =
            h56_region_restricted_node_sets(T2, region.leaves);
        bases.insert(bases.end(), second.begin(), second.end());
        for (vector<int> &base : bases) base = normalize(std::move(base));
        sort(bases.begin(), bases.end(),
             [&](const vector<int> &a, const vector<int> &b) {
                 if (a.size() != b.size()) return a.size() > b.size();
                 int sa = h56_component_score(a, T1, T2);
                 int sb = h56_component_score(b, T1, T2);
                 if (sa != sb) return sa > sb;
                 return a < b;
             });
        bases.erase(unique(bases.begin(), bases.end()), bases.end());
        if ((int)bases.size() > base_limit) bases.resize(base_limit);

        vector<vector<int>> removals = bases;
        sort(removals.begin(), removals.end(),
             [&](const vector<int> &a, const vector<int> &b) {
                 if (a.size() != b.size()) return a.size() < b.size();
                 int sa = h56_component_score(a, T1, T2);
                 int sb = h56_component_score(b, T1, T2);
                 if (sa != sb) return sa > sb;
                 return a < b;
             });
        if ((int)removals.size() > remove_limit)
            removals.resize(remove_limit);

        unordered_set<string> local_seen;
        for (const vector<int> &base : bases) {
            if ((int)pool.size() >= max_pool || accepted >= add_limit ||
                !still_time()) {
                break;
            }
            if ((int)base.size() < min_component_size + 1 ||
                (int)base.size() > max_component_size + max_remove_size) {
                continue;
            }

            string base_key = h45_component_key(base);
            if (local_seen.insert(base_key).second) try_add(base);

            vector<vector<int>> usable_removals;
            usable_removals.reserve(removals.size());
            for (const vector<int> &rem : removals) {
                if (rem.empty() || rem.size() >= base.size() ||
                    (int)rem.size() > max_remove_size) {
                    continue;
                }
                if (!h56_is_subset_sorted(rem, base)) continue;
                usable_removals.push_back(rem);
            }

            for (int i = 0; i < (int)usable_removals.size(); i++) {
                vector<int> candidate =
                    h56_difference(base, usable_removals[i]);
                string key = h45_component_key(candidate);
                if (local_seen.insert(key).second) try_add(candidate);
                if ((int)pool.size() >= max_pool || accepted >= add_limit ||
                    !still_time()) {
                    break;
                }
            }

            const int pair_remove_limit =
                h56_env_int("H56_PRUNED_PAIR_REMOVE_LIMIT",
                            medium_instance ? 30 : 16);
            int pair_cap =
                min((int)usable_removals.size(), pair_remove_limit);
            for (int i = 0; i < pair_cap; i++) {
                for (int j = i + 1; j < pair_cap; j++) {
                    vector<int> removed = usable_removals[i];
                    removed.insert(removed.end(),
                                   usable_removals[j].begin(),
                                   usable_removals[j].end());
                    removed = normalize(std::move(removed));
                    if ((int)removed.size() > max_remove_size * 2 ||
                        !h56_is_subset_sorted(removed, base)) {
                        continue;
                    }
                    vector<int> candidate = h56_difference(base, removed);
                    string key = h45_component_key(candidate);
                    if (local_seen.insert(key).second) try_add(candidate);
                    if ((int)pool.size() >= max_pool ||
                        accepted >= add_limit || !still_time()) {
                        break;
                    }
                }
                if ((int)pool.size() >= max_pool || accepted >= add_limit ||
                    !still_time()) {
                    break;
                }
            }
        }

        if (used_regions >= region_limit) break;
    }

    cerr << "h56 pruned_clades before=" << before
         << " after=" << pool.size()
         << " regions=" << used_regions
         << " accepted=" << accepted
         << " tried=" << tried
         << " elapsed=" << h56_elapsed_seconds(global_start) << '\n';
}

static void h56_add_repartition_growth_components(
    unordered_map<string, vector<int>> &pool,
    const vector<H56Region> &regions,
    const Tree &T1,
    const Tree &T2,
    int max_pool,
    const chrono::steady_clock::time_point &global_start,
    double pool_cutoff_seconds) {
    if (!h56_env_int("H56_ENABLE_REPARTITION_GROWTH", 0)) return;

    const int medium_min =
        h56_env_int("H56_MEDIUM_INSTANCE_MIN_LEAVES",
                    H56_MEDIUM_INSTANCE_MIN_LEAVES_VALUE);
    const bool medium_instance = T1.n >= medium_min;
    const int region_limit =
        h56_env_int("H56_REPARTITION_REGION_LIMIT",
                    medium_instance ? 44 : 28);
    const int seed_limit =
        h56_env_int("H56_REPARTITION_SEED_LIMIT",
                    medium_instance ? 220 : 140);
    const int max_component_size =
        h56_env_int("H56_REPARTITION_MAX_COMPONENT_SIZE",
                    medium_instance ? 42 : 18);
    const int max_extra =
        h56_env_int("H56_REPARTITION_MAX_EXTRA",
                    medium_instance ? 9 : 7);
    const int beam_width =
        h56_env_int("H56_REPARTITION_BEAM_WIDTH",
                    medium_instance ? 24 : 18);
    const int option_limit =
        h56_env_int("H56_REPARTITION_OPTION_LIMIT",
                    medium_instance ? 56 : 36);
    const int add_limit =
        h56_env_int("H56_REPARTITION_ADD_LIMIT",
                    medium_instance ? 1700 : 900);
    const double seconds_limit =
        h56_env_double("H56_REPARTITION_SECONDS",
                       medium_instance ? 9.0 : 5.0);
    const int preferred_size =
        h56_env_int("H56_REPARTITION_PREFERRED_SIZE",
                    medium_instance ? 18 : 10);

    int before = (int)pool.size();
    int accepted = 0;
    int used_regions = 0;
    long long expanded = 0;
    auto start = chrono::steady_clock::now();

    auto still_time = [&]() {
        double local_elapsed = chrono::duration<double>(
            chrono::steady_clock::now() - start).count();
        return local_elapsed < seconds_limit &&
               h56_elapsed_seconds(global_start) < pool_cutoff_seconds &&
               !g_stop_requested;
    };

    auto normalize = [](vector<int> leaves) {
        sort(leaves.begin(), leaves.end());
        leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());
        return leaves;
    };

    for (const H56Region &region : regions) {
        if ((int)pool.size() >= max_pool || accepted >= add_limit ||
            !still_time()) {
            break;
        }
        if (region.leaves.size() < 4) continue;
        used_regions++;

        vector<vector<int>> first_sets =
            h56_region_restricted_node_sets(T1, region.leaves);
        vector<vector<int>> second_sets =
            h56_region_restricted_node_sets(T2, region.leaves);
        vector<vector<int>> seeds;
        unordered_set<string> seed_seen;

        auto add_seed = [&](vector<int> leaves) {
            leaves = normalize(std::move(leaves));
            if ((int)leaves.size() < 2 ||
                (int)leaves.size() > max_component_size) {
                return;
            }
            string key = h45_component_key(leaves);
            if (!seed_seen.insert(key).second) return;
            seeds.push_back(std::move(leaves));
        };

        for (const vector<int> &s : first_sets) add_seed(s);
        for (const vector<int> &s : second_sets) add_seed(s);
        for (const vector<int> &a : first_sets) {
            for (const vector<int> &b : second_sets) {
                add_seed(h56_intersection(a, b));
                add_seed(h56_difference(a, b));
                add_seed(h56_difference(b, a));
            }
        }

        sort(seeds.begin(), seeds.end(),
             [&](const vector<int> &a, const vector<int> &b) {
                 int pa = abs((int)a.size() - preferred_size);
                 int pb = abs((int)b.size() - preferred_size);
                 if (pa != pb) return pa < pb;
                 int sa = h56_component_score(a, T1, T2);
                 int sb = h56_component_score(b, T1, T2);
                 if (sa != sb) return sa > sb;
                 if (a.size() != b.size()) return a.size() > b.size();
                 return a < b;
             });
        if ((int)seeds.size() > seed_limit) seeds.resize(seed_limit);

        vector<char> in_region(T1.n + 1, 0);
        for (int leaf : region.leaves) in_region[leaf] = 1;
        unordered_set<string> region_seen;

        auto try_add = [&](const vector<int> &leaves) {
            size_t old_size = pool.size();
            if (h56_add_candidate(pool, leaves, T1, T2, max_pool)) {
                accepted += pool.size() != old_size;
            }
        };

        for (const vector<int> &seed : seeds) {
            if ((int)pool.size() >= max_pool || accepted >= add_limit ||
                !still_time()) {
                break;
            }
            try_add(seed);

            vector<vector<int>> beam{seed};
            region_seen.insert(h45_component_key(seed));
            int target_cap = min(max_component_size,
                                 (int)seed.size() +
                                     max(1, max_extra));
            if ((int)seed.size() <= 4)
                target_cap = min(target_cap, 16);

            while (!beam.empty() && (int)pool.size() < max_pool &&
                   accepted < add_limit && still_time()) {
                vector<vector<int>> next_states;
                for (const vector<int> &state : beam) {
                    if ((int)state.size() >= target_cap) continue;
                    vector<char> in_state(T1.n + 1, 0);
                    for (int leaf : state) in_state[leaf] = 1;

                    vector<pair<int, int>> options;
                    for (int leaf : region.leaves) {
                        if (in_state[leaf]) continue;
                        int nearest = 1 << 28;
                        for (int existing : state) {
                            nearest = min(nearest, max(
                                h56_leaf_distance(T1, existing, leaf),
                                h56_leaf_distance(T2, existing, leaf)));
                        }
                        options.push_back({nearest, leaf});
                    }
                    sort(options.begin(), options.end(),
                         [](const pair<int, int> &a,
                            const pair<int, int> &b) {
                             if (a.first != b.first) return a.first < b.first;
                             return a.second < b.second;
                         });
                    if ((int)options.size() > option_limit)
                        options.resize(option_limit);

                    for (const auto &option : options) {
                        vector<int> grown = state;
                        grown.push_back(option.second);
                        grown = normalize(std::move(grown));
                        if ((int)grown.size() != (int)state.size() + 1)
                            continue;
                        string key = h45_component_key(grown);
                        if (!region_seen.insert(key).second) continue;
                        expanded++;
                        try_add(grown);
                        next_states.push_back(std::move(grown));
                        if ((int)pool.size() >= max_pool ||
                            accepted >= add_limit || !still_time()) {
                            break;
                        }
                    }
                    if ((int)pool.size() >= max_pool ||
                        accepted >= add_limit || !still_time()) {
                        break;
                    }
                }

                sort(next_states.begin(), next_states.end(),
                     [&](const vector<int> &a, const vector<int> &b) {
                         int pa = abs((int)a.size() - preferred_size);
                         int pb = abs((int)b.size() - preferred_size);
                         if (pa != pb) return pa < pb;
                         int sa = h56_component_score(a, T1, T2);
                         int sb = h56_component_score(b, T1, T2);
                         if (sa != sb) return sa > sb;
                         if (a.size() != b.size())
                             return a.size() > b.size();
                         return a < b;
                     });
                next_states.erase(unique(next_states.begin(),
                                         next_states.end()),
                                  next_states.end());
                if ((int)next_states.size() > beam_width)
                    next_states.resize(beam_width);
                beam = std::move(next_states);
            }
        }

        if (used_regions >= region_limit) break;
    }

    cerr << "h56 repartition_growth before=" << before
         << " after=" << pool.size()
         << " regions=" << used_regions
         << " accepted=" << accepted
         << " expanded=" << expanded
         << " elapsed=" << h56_elapsed_seconds(global_start) << '\n';
}

static void h56_small_subset_dfs(
    const vector<int> &leaves,
    int target_size,
    int start_index,
    vector<int> &current,
    unordered_map<string, vector<int>> &pool,
    unordered_set<string> &local_seen,
    const Tree &T1,
    const Tree &T2,
    int max_pool,
    int added_limit,
    long long visited_limit,
    int &added,
    long long &visited,
    const chrono::steady_clock::time_point &global_start,
    double pool_cutoff_seconds) {
    if ((int)pool.size() >= max_pool || added >= added_limit ||
        visited >= visited_limit ||
        h56_elapsed_seconds(global_start) >= pool_cutoff_seconds) {
        return;
    }
    if ((int)current.size() == target_size) {
        visited++;
        string key = h45_component_key(current);
        if (!local_seen.insert(key).second) return;
        size_t before = pool.size();
        h56_add_candidate(pool, current, T1, T2, max_pool);
        if (pool.size() != before) added++;
        return;
    }

    int need = target_size - (int)current.size();
    int last = (int)leaves.size() - need;
    for (int index = start_index; index <= last; index++) {
        current.push_back(leaves[index]);
        h56_small_subset_dfs(leaves, target_size, index + 1, current,
                             pool, local_seen, T1, T2, max_pool,
                             added_limit, visited_limit, added, visited,
                             global_start, pool_cutoff_seconds);
        current.pop_back();
        if ((int)pool.size() >= max_pool || added >= added_limit ||
            visited >= visited_limit ||
            h56_elapsed_seconds(global_start) >= pool_cutoff_seconds) {
            return;
        }
    }
}

static void h56_add_small_region_subset_components(
    unordered_map<string, vector<int>> &pool,
    const vector<H56Region> &regions,
    const Tree &T1,
    const Tree &T2,
    int max_pool,
    const chrono::steady_clock::time_point &global_start,
    double pool_cutoff_seconds) {
    const int region_limit = h56_env_int("H56_SMALL_SUBSET_REGION_LIMIT", 48);
    const int enum_leaf_limit =
        h56_env_int("H56_SMALL_SUBSET_LEAF_LIMIT", 18);
    const int max_subset_size =
        h56_env_int("H56_SMALL_SUBSET_MAX_SIZE", enum_leaf_limit);
    const int min_subset_size =
        h56_env_int("H56_SMALL_SUBSET_MIN_SIZE", 2);
    const int per_region_limit =
        h56_env_int("H56_SMALL_SUBSET_PER_REGION", 120);
    const long long mask_limit =
        (long long)h56_env_int("H56_SMALL_SUBSET_MASK_LIMIT", 45000);

    int before = (int)pool.size();
    int used_regions = 0;
    int added = 0;
    long long visited = 0;

    for (const H56Region &region : regions) {
        if ((int)pool.size() >= max_pool ||
            h56_elapsed_seconds(global_start) >= pool_cutoff_seconds) {
            break;
        }
        vector<int> leaves = region.leaves;
        sort(leaves.begin(), leaves.end());
        leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());
        if ((int)leaves.size() < min_subset_size ||
            (int)leaves.size() > enum_leaf_limit) {
            continue;
        }

        used_regions++;
        unordered_set<string> local_seen;
        local_seen.reserve((size_t)per_region_limit * 4 + 16);
        vector<int> current;
        current.reserve(leaves.size());
        int local_added_before = added;
        long long local_visited_before = visited;

        int top_size = min((int)leaves.size(), max_subset_size);
        for (int target_size = top_size;
             target_size >= min_subset_size &&
             (int)pool.size() < max_pool &&
             added - local_added_before < per_region_limit &&
             visited - local_visited_before < mask_limit &&
             h56_elapsed_seconds(global_start) < pool_cutoff_seconds;
             target_size--) {
            h56_small_subset_dfs(leaves, target_size, 0, current, pool,
                                 local_seen, T1, T2, max_pool,
                                 local_added_before + per_region_limit,
                                 local_visited_before + mask_limit,
                                 added, visited, global_start,
                                 pool_cutoff_seconds);
        }

        if (used_regions >= region_limit) break;
    }

    cerr << "h56 small_region_subsets before=" << before
         << " after=" << pool.size()
         << " regions=" << used_regions
         << " added=" << added
         << " visited=" << visited
         << " elapsed=" << h56_elapsed_seconds(global_start) << '\n';
}

static bool h56_collect_leaves_limited(const Tree &T, int node, int limit,
                                       vector<int> &leaves) {
    leaves.clear();
    vector<int> stack{node};
    while (!stack.empty()) {
        int current = stack.back();
        stack.pop_back();
        const TreeNode &tn = T.nodes[current];
        if (tn.label >= 1) {
            leaves.push_back(tn.label);
            if ((int)leaves.size() > limit) return false;
        } else {
            for (int child : tn.child) stack.push_back(child);
        }
    }
    sort(leaves.begin(), leaves.end());
    return !leaves.empty();
}

static vector<vector<int>> h56_global_seed_pairs(const Tree &T1,
                                                 const Tree &T2,
                                                 int leaf_count) {
    const int clade_cap = h56_env_int("H56_GROWTH_SEED_CLADE_CAP", 28);
    const int seed_limit = h56_env_int("H56_GROWTH_SEED_LIMIT", 12000);
    unordered_set<string> seen;
    vector<vector<int>> seeds;
    vector<int> leaves;

    auto scan_tree = [&](const Tree &T) {
        for (int node = 0; node < (int)T.nodes.size(); node++) {
            if ((int)seeds.size() >= seed_limit * 4) return;
            if (!h56_collect_leaves_limited(T, node, clade_cap, leaves))
                continue;
            if ((int)leaves.size() < 2 || (int)leaves.size() > clade_cap)
                continue;
            for (int i = 0; i < (int)leaves.size(); i++) {
                for (int j = i + 1; j < (int)leaves.size(); j++) {
                    vector<int> pair{leaves[i], leaves[j]};
                    if (!h56_component_valid(pair, T1, T2)) continue;
                    string key = h45_component_key(pair);
                    if (!seen.insert(key).second) continue;
                    seeds.push_back(std::move(pair));
                }
            }
        }
    };

    scan_tree(T1);
    scan_tree(T2);
    sort(seeds.begin(), seeds.end(),
         [&](const vector<int> &a, const vector<int> &b) {
             int sa = h56_component_score(a, T1, T2);
             int sb = h56_component_score(b, T1, T2);
             if (sa != sb) return sa > sb;
             return a < b;
         });
    if ((int)seeds.size() > seed_limit) seeds.resize(seed_limit);
    return seeds;
}

static vector<int> h56_global_neighbor_leaves(const vector<int> &state,
                                              const Tree &T1,
                                              const Tree &T2,
                                              int leaf_count) {
    const int ancestor_steps = h56_env_int("H56_GROWTH_ANCESTOR_STEPS", 5);
    const int neighborhood_cap =
        h56_env_int("H56_GROWTH_NEIGHBORHOOD_CAP", 96);
    const int neighbor_limit = h56_env_int("H56_GROWTH_NEIGHBOR_LIMIT", 120);
    vector<char> in_state(leaf_count + 1, 0);
    vector<char> seen(leaf_count + 1, 0);
    vector<int> neighbors;
    vector<int> leaves;
    for (int leaf : state) in_state[leaf] = 1;

    auto scan_tree = [&](const Tree &T) {
        int node = lca_of_leaf_set(T, state);
        for (int step = 0;
             node >= 0 && step <= ancestor_steps &&
             (int)neighbors.size() < neighbor_limit;
             step++, node = T.nodes[node].parent) {
            if (!h56_collect_leaves_limited(T, node, neighborhood_cap, leaves))
                continue;
            for (int leaf : leaves) {
                if (in_state[leaf] || seen[leaf]) continue;
                seen[leaf] = 1;
                neighbors.push_back(leaf);
                if ((int)neighbors.size() >= neighbor_limit) break;
            }
        }
    };

    scan_tree(T1);
    scan_tree(T2);
    sort(neighbors.begin(), neighbors.end(), [&](int a, int b) {
        vector<int> aa = state;
        vector<int> bb = state;
        aa.push_back(a);
        bb.push_back(b);
        sort(aa.begin(), aa.end());
        sort(bb.begin(), bb.end());
        int sa = h56_component_score(aa, T1, T2);
        int sb = h56_component_score(bb, T1, T2);
        if (sa != sb) return sa > sb;
        return a < b;
    });
    return neighbors;
}

static void h56_add_global_growth_components(
    unordered_map<string, vector<int>> &pool,
    const Tree &T1, const Tree &T2,
    int max_pool) {
    const int leaf_count = T1.n;
    const int max_component_size =
        h56_env_int("H56_GROWTH_MAX_COMPONENT_SIZE", 16);
    const int beam_width = h56_env_int("H56_GROWTH_BEAM_WIDTH", 20);
    const double seconds_limit =
        h56_env_double("H56_GROWTH_SECONDS", 18.0);
    auto start = chrono::steady_clock::now();
    vector<vector<int>> seeds = h56_global_seed_pairs(T1, T2, leaf_count);
    int before = (int)pool.size();
    int accepted = 0;
    long long expansions = 0;

    for (const vector<int> &seed : seeds) {
        if ((int)pool.size() >= max_pool) break;
        double elapsed = chrono::duration<double>(
            chrono::steady_clock::now() - start).count();
        if (elapsed >= seconds_limit) break;

        h56_add_candidate(pool, seed, T1, T2, max_pool);
        vector<vector<int>> beam{seed};
        unordered_set<string> local_seen;
        local_seen.insert(h45_component_key(seed));

        while (!beam.empty() && (int)pool.size() < max_pool) {
            vector<vector<int>> next_states;
            for (const vector<int> &state : beam) {
                if ((int)state.size() >= max_component_size) continue;
                vector<int> neighbors =
                    h56_global_neighbor_leaves(state, T1, T2, leaf_count);
                for (int leaf : neighbors) {
                    vector<int> grown = state;
                    grown.push_back(leaf);
                    sort(grown.begin(), grown.end());
                    grown.erase(unique(grown.begin(), grown.end()),
                                grown.end());
                    if ((int)grown.size() != (int)state.size() + 1)
                        continue;
                    string key = h45_component_key(grown);
                    if (!local_seen.insert(key).second) continue;
                    expansions++;
                    size_t old_size = pool.size();
                    if (h56_add_candidate(pool, grown, T1, T2, max_pool)) {
                        accepted += pool.size() != old_size;
                        next_states.push_back(std::move(grown));
                    }
                    if ((int)pool.size() >= max_pool) break;
                }
                if ((int)pool.size() >= max_pool) break;
            }
            sort(next_states.begin(), next_states.end(),
                 [&](const vector<int> &a, const vector<int> &b) {
                     int sa = h56_component_score(a, T1, T2);
                     int sb = h56_component_score(b, T1, T2);
                     if (sa != sb) return sa > sb;
                     if (a.size() != b.size()) return a.size() > b.size();
                     return a < b;
                 });
            next_states.erase(unique(next_states.begin(),
                                     next_states.end()),
                              next_states.end());
            if ((int)next_states.size() > beam_width)
                next_states.resize(beam_width);
            beam = std::move(next_states);

            elapsed = chrono::duration<double>(
                chrono::steady_clock::now() - start).count();
            if (elapsed >= seconds_limit) break;
        }
    }

    cerr << "h56 global_growth before=" << before
         << " after=" << pool.size()
         << " seeds=" << seeds.size()
         << " accepted=" << accepted
         << " expansions=" << expansions << '\n';
}

static void h56_add_distance_clique_components(
    unordered_map<string, vector<int>> &pool,
    const Tree &T1, const Tree &T2,
    int max_pool) {
    const int n = T1.n;
    const int max_n = h56_env_int("H56_CLIQUE_MAX_N", 1500);
    if (n > max_n) {
        cerr << "h56 clique_skip leaves=" << n << " max_n=" << max_n << '\n';
        return;
    }
    const int distance_limit = h56_env_int("H56_CLIQUE_DISTANCE", 22);
    const int max_component_size =
        h56_env_int("H56_CLIQUE_MAX_COMPONENT_SIZE", 16);
    const int seed_limit = h56_env_int("H56_CLIQUE_SEED_LIMIT", 24000);
    const int beam_width = h56_env_int("H56_CLIQUE_BEAM_WIDTH", 20);
    const double seconds_limit =
        h56_env_double("H56_CLIQUE_SECONDS", 18.0);
    auto start = chrono::steady_clock::now();

    vector<vector<int>> adjacency(n + 1);
    vector<vector<int>> seeds;
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            int distance = max(h56_leaf_distance(T1, a, b),
                               h56_leaf_distance(T2, a, b));
            if (distance > distance_limit) continue;
            adjacency[a].push_back(b);
            adjacency[b].push_back(a);
            seeds.push_back(vector<int>{a, b});
        }
    }
    for (int leaf = 1; leaf <= n; leaf++)
        sort(adjacency[leaf].begin(), adjacency[leaf].end());

    sort(seeds.begin(), seeds.end(),
         [&](const vector<int> &a, const vector<int> &b) {
             int sa = h56_component_score(a, T1, T2);
             int sb = h56_component_score(b, T1, T2);
             if (sa != sb) return sa > sb;
             return a < b;
         });
    if ((int)seeds.size() > seed_limit) seeds.resize(seed_limit);
    {
        const int stride = max(1, h56_env_int("H56_CLIQUE_SEED_STRIDE", 257));
        vector<vector<int>> scheduled;
        scheduled.reserve(seeds.size());
        for (int offset = 0; offset < stride; offset++) {
            for (int index = offset; index < (int)seeds.size();
                 index += stride) {
                scheduled.push_back(std::move(seeds[index]));
            }
        }
        seeds = std::move(scheduled);
    }

    unordered_set<string> global_seen;
    int before = (int)pool.size();
    int accepted = 0;
    long long expanded = 0;

    auto adjacent_to_all = [&](int leaf, const vector<int> &state) {
        const vector<int> &neighbors = adjacency[leaf];
        for (int existing : state) {
            if (!binary_search(neighbors.begin(), neighbors.end(), existing))
                return false;
        }
        return true;
    };

    auto run_closure_order = [&](const vector<int> &seed,
                                 vector<int> options) {
        vector<int> state = seed;
        h56_add_candidate(pool, state, T1, T2, max_pool);
        for (int leaf : options) {
            if ((int)state.size() >= max_component_size ||
                (int)pool.size() >= max_pool) {
                break;
            }
            if (binary_search(state.begin(), state.end(), leaf))
                continue;
            if (!adjacent_to_all(leaf, state)) continue;
            vector<int> grown = state;
            grown.push_back(leaf);
            sort(grown.begin(), grown.end());
            string key = h45_component_key(grown);
            if (!global_seen.insert(key).second) continue;
            expanded++;
            size_t old_size = pool.size();
            if (h56_add_candidate(pool, grown, T1, T2, max_pool)) {
                accepted += pool.size() != old_size;
                state = std::move(grown);
            }
        }
    };

    for (const vector<int> &seed : seeds) {
        if ((int)pool.size() >= max_pool) break;
        double elapsed = chrono::duration<double>(
            chrono::steady_clock::now() - start).count();
        if (elapsed >= seconds_limit) break;

        vector<int> common;
        const vector<int> &first_neighbors = adjacency[seed[0]];
        const vector<int> &second_neighbors = adjacency[seed[1]];
        set_intersection(first_neighbors.begin(), first_neighbors.end(),
                         second_neighbors.begin(), second_neighbors.end(),
                         back_inserter(common));
        common.erase(remove(common.begin(), common.end(), seed[0]),
                     common.end());
        common.erase(remove(common.begin(), common.end(), seed[1]),
                     common.end());
        vector<pair<int, int>> ranked_common;
        ranked_common.reserve(common.size());
        for (int leaf : common) {
            int degree = 0;
            const vector<int> &neighbors = adjacency[leaf];
            for (int other : common)
                degree += binary_search(neighbors.begin(), neighbors.end(),
                                        other);
            ranked_common.push_back({-degree, leaf});
        }
        sort(ranked_common.begin(), ranked_common.end());
        vector<int> degree_order;
        degree_order.reserve(ranked_common.size());
        for (const auto &entry : ranked_common)
            degree_order.push_back(entry.second);
        vector<int> leaf_order = common;
        sort(leaf_order.begin(), leaf_order.end());
        run_closure_order(seed, leaf_order);
        reverse(leaf_order.begin(), leaf_order.end());
        run_closure_order(seed, leaf_order);
        run_closure_order(seed, degree_order);
        reverse(degree_order.begin(), degree_order.end());
        run_closure_order(seed, degree_order);
        if ((int)pool.size() >= max_pool) break;

        vector<vector<int>> beam{seed};
        global_seen.insert(h45_component_key(seed));
        h56_add_candidate(pool, seed, T1, T2, max_pool);

        while (!beam.empty() && (int)pool.size() < max_pool) {
            vector<vector<int>> next_states;
            for (const vector<int> &state : beam) {
                if ((int)state.size() >= max_component_size) continue;
                vector<pair<int, int>> options;
                for (int leaf = 1; leaf <= n; leaf++) {
                    if (binary_search(state.begin(), state.end(), leaf))
                        continue;
                    if (!adjacent_to_all(leaf, state)) continue;
                    vector<int> grown = state;
                    grown.push_back(leaf);
                    sort(grown.begin(), grown.end());
                    int score = h56_component_score(grown, T1, T2);
                    options.push_back({-score, leaf});
                }
                sort(options.begin(), options.end());
                const int option_limit =
                    h56_env_int("H56_CLIQUE_OPTION_LIMIT", 96);
                if ((int)options.size() > option_limit)
                    options.resize(option_limit);

                for (const auto &option : options) {
                    vector<int> grown = state;
                    grown.push_back(option.second);
                    sort(grown.begin(), grown.end());
                    grown.erase(unique(grown.begin(), grown.end()),
                                grown.end());
                    string key = h45_component_key(grown);
                    if (!global_seen.insert(key).second) continue;
                    expanded++;
                    size_t old_size = pool.size();
                    if (h56_add_candidate(pool, grown, T1, T2, max_pool)) {
                        accepted += pool.size() != old_size;
                        next_states.push_back(std::move(grown));
                    }
                    if ((int)pool.size() >= max_pool) break;
                }
                if ((int)pool.size() >= max_pool) break;
            }

            sort(next_states.begin(), next_states.end(),
                 [&](const vector<int> &a, const vector<int> &b) {
                     int sa = h56_component_score(a, T1, T2);
                     int sb = h56_component_score(b, T1, T2);
                     if (sa != sb) return sa > sb;
                     if (a.size() != b.size()) return a.size() > b.size();
                     return a < b;
                 });
            next_states.erase(unique(next_states.begin(),
                                     next_states.end()),
                              next_states.end());
            if ((int)next_states.size() > beam_width)
                next_states.resize(beam_width);
            beam = std::move(next_states);

            elapsed = chrono::duration<double>(
                chrono::steady_clock::now() - start).count();
            if (elapsed >= seconds_limit) break;
        }
    }
    cerr << "h56 clique_growth before=" << before
         << " after=" << pool.size()
         << " seeds=" << seeds.size()
         << " accepted=" << accepted
         << " expanded=" << expanded << '\n';
}

static void h56_add_rank_bucket_beam_components(
    unordered_map<string, vector<int>> &pool,
    const Tree &T1, const Tree &T2,
    int max_pool) {
    const int n = T1.n;
    const int max_n = h56_env_int("H56_BUCKET_MAX_N", 1500);
    if (n > max_n) {
        cerr << "h56 bucket_beam_skip leaves=" << n
             << " max_n=" << max_n << '\n';
        return;
    }

    const int distance_limit = h56_env_int("H56_BUCKET_DISTANCE", 22);
    const int max_component_size =
        h56_env_int("H56_BUCKET_MAX_COMPONENT_SIZE", 16);
    const int beam_width = h56_env_int("H56_BUCKET_BEAM_WIDTH", 20);
    const int option_limit = h56_env_int("H56_BUCKET_OPTION_LIMIT", 40);
    const int bucket_step = max(1, h56_env_int("H56_BUCKET_STEP", 250));
    const int bucket_half = max(0, h56_env_int("H56_BUCKET_HALF_WIDTH", 20));
    const int selected_limit =
        h56_env_int("H56_BUCKET_SEED_LIMIT", 6500);
    const int add_limit = h56_env_int("H56_BUCKET_ADD_LIMIT", 2800);
    const int min_add_size = h56_env_int("H56_BUCKET_MIN_ADD_SIZE", 4);
    const int per_seed_add_limit =
        h56_env_int("H56_BUCKET_PER_SEED_ADD_LIMIT", 8);
    const int per_seed_expand_limit =
        h56_env_int("H56_BUCKET_PER_SEED_EXPAND_LIMIT", 2200);
    const double seconds_limit =
        h56_env_double("H56_BUCKET_SECONDS", 24.0);
    auto start = chrono::steady_clock::now();

    vector<vector<int>> adjacency(n + 1);
    vector<vector<int>> seeds;
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            int distance = max(h56_leaf_distance(T1, a, b),
                               h56_leaf_distance(T2, a, b));
            if (distance > distance_limit) continue;
            adjacency[a].push_back(b);
            adjacency[b].push_back(a);
            seeds.push_back(vector<int>{a, b});
        }
    }
    for (int leaf = 1; leaf <= n; leaf++)
        sort(adjacency[leaf].begin(), adjacency[leaf].end());

    sort(seeds.begin(), seeds.end(),
         [&](const vector<int> &a, const vector<int> &b) {
             int sa = h56_component_score(a, T1, T2);
             int sb = h56_component_score(b, T1, T2);
             if (sa != sb) return sa > sb;
             return a < b;
         });

    vector<vector<int>> selected;
    selected.reserve(min((int)seeds.size(), selected_limit));
    vector<char> selected_index(seeds.size(), 0);
    auto add_seed_index = [&](int index) {
        if (index < 0 || index >= (int)seeds.size()) return;
        if (selected_index[index]) return;
        selected_index[index] = 1;
        selected.push_back(seeds[index]);
    };
    int bucket_count = ((int)seeds.size() + bucket_step - 1) / bucket_step;
    for (int delta = 0;
         delta <= bucket_half && (int)selected.size() < selected_limit;
         delta++) {
        for (int bucket = 0;
             bucket < bucket_count && (int)selected.size() < selected_limit;
             bucket++) {
            int center = bucket * bucket_step;
            add_seed_index(center + delta);
            if (delta > 0) add_seed_index(center - delta);
        }
    }

    auto adjacent_to_all = [&](int leaf, const vector<int> &state) {
        const vector<int> &neighbors = adjacency[leaf];
        for (int existing : state) {
            if (!binary_search(neighbors.begin(), neighbors.end(), existing))
                return false;
        }
        return true;
    };

    int before = (int)pool.size();
    int accepted = 0;
    long long expanded = 0;
    unordered_set<string> global_seen;

    for (const vector<int> &seed : selected) {
        if ((int)pool.size() >= max_pool || accepted >= add_limit) break;
        double elapsed = chrono::duration<double>(
            chrono::steady_clock::now() - start).count();
        if (elapsed >= seconds_limit) break;

        vector<vector<int>> beam{seed};
        unordered_set<string> local_seen;
        local_seen.insert(h45_component_key(seed));
        vector<vector<int>> seed_candidates;
        long long seed_expanded = 0;

        while (!beam.empty() && (int)pool.size() < max_pool &&
               accepted < add_limit) {
            vector<vector<int>> next_states;
            for (const vector<int> &state : beam) {
                if ((int)state.size() >= max_component_size) continue;

                vector<pair<int, int>> options;
                for (int leaf = 1; leaf <= n; leaf++) {
                    if (binary_search(state.begin(), state.end(), leaf))
                        continue;
                    if (!adjacent_to_all(leaf, state)) continue;
                    vector<int> grown = state;
                    grown.push_back(leaf);
                    sort(grown.begin(), grown.end());
                    if (!h56_component_valid(grown, T1, T2)) continue;
                    options.push_back(
                        {-h56_component_score(grown, T1, T2), leaf});
                }
                sort(options.begin(), options.end());
                if ((int)options.size() > option_limit)
                    options.resize(option_limit);

                for (const auto &option : options) {
                    vector<int> grown = state;
                    grown.push_back(option.second);
                    sort(grown.begin(), grown.end());
                    grown.erase(unique(grown.begin(), grown.end()),
                                grown.end());
                    string key = h45_component_key(grown);
                    if (!local_seen.insert(key).second) continue;
                    expanded++;
                    seed_expanded++;
                    if ((int)grown.size() >= min_add_size)
                        seed_candidates.push_back(grown);
                    next_states.push_back(std::move(grown));
                    if ((int)pool.size() >= max_pool ||
                        accepted >= add_limit ||
                        seed_expanded >= per_seed_expand_limit) {
                        break;
                    }
                }
                if ((int)pool.size() >= max_pool || accepted >= add_limit ||
                    seed_expanded >= per_seed_expand_limit)
                    break;
            }

            sort(next_states.begin(), next_states.end(),
                 [&](const vector<int> &a, const vector<int> &b) {
                     int sa = h56_component_score(a, T1, T2);
                     int sb = h56_component_score(b, T1, T2);
                     if (sa != sb) return sa > sb;
                     if (a.size() != b.size()) return a.size() > b.size();
                     return a < b;
                 });
            next_states.erase(unique(next_states.begin(),
                                     next_states.end()),
                              next_states.end());
            if ((int)next_states.size() > beam_width)
                next_states.resize(beam_width);
            beam = std::move(next_states);

            elapsed = chrono::duration<double>(
                chrono::steady_clock::now() - start).count();
            if (elapsed >= seconds_limit ||
                seed_expanded >= per_seed_expand_limit) {
                break;
            }
        }

        sort(seed_candidates.begin(), seed_candidates.end(),
             [&](const vector<int> &a, const vector<int> &b) {
                 int da = abs((int)a.size() -
                              h56_env_int("H56_BUCKET_PREFERRED_SIZE", 10));
                 int db = abs((int)b.size() -
                              h56_env_int("H56_BUCKET_PREFERRED_SIZE", 10));
                 if (da != db) return da < db;
                 int sa = h56_component_score(a, T1, T2);
                 int sb = h56_component_score(b, T1, T2);
                 if (sa != sb) return sa > sb;
                 return a < b;
             });
        seed_candidates.erase(unique(seed_candidates.begin(),
                                     seed_candidates.end()),
                              seed_candidates.end());

        int seed_added = 0;
        for (const vector<int> &candidate : seed_candidates) {
            if (seed_added >= per_seed_add_limit ||
                accepted >= add_limit ||
                (int)pool.size() >= max_pool) {
                break;
            }
            string key = h45_component_key(candidate);
            if (!global_seen.insert(key).second) continue;
            size_t old_size = pool.size();
            if (h56_add_candidate(pool, candidate, T1, T2, max_pool)) {
                int delta = pool.size() != old_size;
                accepted += delta;
                seed_added += delta;
            }
        }
    }

    cerr << "h56 bucket_beam before=" << before
         << " after=" << pool.size()
         << " seeds=" << seeds.size()
         << " selected=" << selected.size()
         << " accepted=" << accepted
         << " expanded=" << expanded << '\n';
}

static bool h56_choose_rec(const vector<int> &items,
                           int choose,
                           int start,
                           vector<int> &current,
                           const function<bool(const vector<int>&)> &visit) {
    if ((int)current.size() == choose) return visit(current);
    int need = choose - (int)current.size();
    for (int i = start; i <= (int)items.size() - need; i++) {
        current.push_back(items[i]);
        if (!h56_choose_rec(items, choose, i + 1, current, visit))
            return false;
        current.pop_back();
    }
    return true;
}

static bool h56_adjacent_to_all_leaves(const vector<vector<int>> &adjacency,
                                       int leaf,
                                       const vector<int> &state) {
    const vector<int> &neighbors = adjacency[leaf];
    for (int existing : state) {
        if (!binary_search(neighbors.begin(), neighbors.end(), existing))
            return false;
    }
    return true;
}

static pair<int, int> h56_component_badness(const vector<int> &leaves,
                                            const Tree &T1,
                                            const Tree &T2) {
    if (leaves.size() <= 2) return {0, 0};
    int min_total = INT_MAX;
    int max_total = 0;
    int max_single = 0;
    for (int leaf : leaves) {
        int total = 0;
        int local_max = 0;
        for (int other : leaves) {
            if (leaf == other) continue;
            int distance = max(h56_leaf_distance(T1, leaf, other),
                               h56_leaf_distance(T2, leaf, other));
            total += distance;
            local_max = max(local_max, distance);
        }
        min_total = min(min_total, total);
        max_total = max(max_total, total);
        max_single = max(max_single, local_max);
    }
    return {max_total - min_total, max_single};
}

static void h56_keep_top_edits_for_remove_set(
    vector<vector<int>> &base_candidates,
    unordered_set<string> &base_seen,
    vector<vector<int>> remove_candidates,
    const Tree &T1,
    const Tree &T2,
    int per_remove_limit) {
    sort(remove_candidates.begin(), remove_candidates.end(),
         [&](const vector<int> &a, const vector<int> &b) {
             int sa = h56_component_score(a, T1, T2);
             int sb = h56_component_score(b, T1, T2);
             if (sa != sb) return sa > sb;
             if (a.size() != b.size()) return a.size() > b.size();
             return a < b;
         });
    remove_candidates.erase(unique(remove_candidates.begin(),
                                   remove_candidates.end()),
                            remove_candidates.end());
    int kept = 0;
    for (const vector<int> &candidate : remove_candidates) {
        if (kept >= per_remove_limit) break;
        string key = h45_component_key(candidate);
        if (!base_seen.insert(key).second) continue;
        base_candidates.push_back(candidate);
        kept++;
    }
}

static void h56_add_swap_refinement_components(
    unordered_map<string, vector<int>> &pool,
    const Tree &T1, const Tree &T2,
    int max_pool) {
    const int n = T1.n;
    const int max_n = h56_env_int("H56_SWAP_MAX_N", 1500);
    if (n > max_n) {
        cerr << "h56 swap_refine_skip leaves=" << n
             << " max_n=" << max_n << '\n';
        return;
    }

    const int distance_limit = h56_env_int("H56_SWAP_DISTANCE", 22);
    const int max_remove = h56_env_int("H56_SWAP_MAX_REMOVE", 3);
    const int min_base_size = h56_env_int("H56_SWAP_MIN_BASE_SIZE", 4);
    const int max_base_size = h56_env_int("H56_SWAP_MAX_BASE_SIZE", 13);
    const int option_limit = h56_env_int("H56_SWAP_OPTION_LIMIT", 80);
    const int base_limit = h56_env_int("H56_SWAP_BASE_LIMIT", 4400);
    const int per_base_limit = h56_env_int("H56_SWAP_PER_BASE_LIMIT", 2);
    const int add_beam_width =
        h56_env_int("H56_SWAP_ADD_BEAM_WIDTH", 20);
    const int add_option_limit =
        h56_env_int("H56_SWAP_ADD_OPTION_LIMIT", 24);
    const int add_limit = h56_env_int("H56_SWAP_ADD_LIMIT", 2200);
    const int preferred_size =
        h56_env_int("H56_SWAP_PREFERRED_SIZE", 10);
    const double seconds_limit =
        h56_env_double("H56_SWAP_SECONDS", 16.0);
    auto start = chrono::steady_clock::now();

    vector<vector<int>> adjacency(n + 1);
    for (int a = 1; a <= n; a++) {
        for (int b = a + 1; b <= n; b++) {
            int distance = max(h56_leaf_distance(T1, a, b),
                               h56_leaf_distance(T2, a, b));
            if (distance > distance_limit) continue;
            adjacency[a].push_back(b);
            adjacency[b].push_back(a);
        }
    }
    for (int leaf = 1; leaf <= n; leaf++)
        sort(adjacency[leaf].begin(), adjacency[leaf].end());

    vector<vector<int>> bases;
    bases.reserve(pool.size());
    for (const auto &entry : pool) {
        const vector<int> &component = entry.second;
        if ((int)component.size() < min_base_size ||
            (int)component.size() > max_base_size) {
            continue;
        }
        bases.push_back(component);
    }
    sort(bases.begin(), bases.end());
    bases.erase(unique(bases.begin(), bases.end()), bases.end());

    vector<vector<int>> preferred_order = bases;
    sort(preferred_order.begin(), preferred_order.end(),
         [&](const vector<int> &a, const vector<int> &b) {
             int da = abs((int)a.size() - preferred_size);
             int db = abs((int)b.size() - preferred_size);
             if (da != db) return da < db;
             int sa = h56_component_score(a, T1, T2);
             int sb = h56_component_score(b, T1, T2);
             if (sa != sb) return sa > sb;
             return a < b;
         });

    vector<vector<int>> badness_order = bases;
    sort(badness_order.begin(), badness_order.end(),
         [&](const vector<int> &a, const vector<int> &b) {
             pair<int, int> ba = h56_component_badness(a, T1, T2);
             pair<int, int> bb = h56_component_badness(b, T1, T2);
             if (ba.first != bb.first) return ba.first > bb.first;
             if (ba.second != bb.second) return ba.second > bb.second;
             int da = abs((int)a.size() - preferred_size);
             int db = abs((int)b.size() - preferred_size);
             if (da != db) return da < db;
             return a < b;
         });

    vector<vector<int>> small_order = bases;
    sort(small_order.begin(), small_order.end(),
         [&](const vector<int> &a, const vector<int> &b) {
             if (a.size() != b.size()) return a.size() < b.size();
             pair<int, int> ba = h56_component_badness(a, T1, T2);
             pair<int, int> bb = h56_component_badness(b, T1, T2);
             if (ba.first != bb.first) return ba.first > bb.first;
             return a < b;
         });

    vector<vector<int>> scheduled_bases;
    scheduled_bases.reserve(min(base_limit, (int)bases.size()));
    unordered_set<string> scheduled_seen;
    auto append_order = [&](const vector<vector<int>> &ordered,
                            int limit) {
        for (const vector<int> &base : ordered) {
            if ((int)scheduled_bases.size() >= base_limit ||
                limit <= 0) {
                break;
            }
            string key = h45_component_key(base);
            if (!scheduled_seen.insert(key).second) continue;
            scheduled_bases.push_back(base);
            limit--;
        }
    };
    append_order(preferred_order,
                 h56_env_int("H56_SWAP_PREFERRED_BASES", 1400));
    append_order(badness_order, h56_env_int("H56_SWAP_BAD_BASES", 2200));
    append_order(small_order, h56_env_int("H56_SWAP_SMALL_BASES", 800));
    bases = std::move(scheduled_bases);

    int before = (int)pool.size();
    int accepted = 0;
    long long tried = 0;
    unordered_set<string> generated;

    for (const vector<int> &base : bases) {
        if ((int)pool.size() >= max_pool || accepted >= add_limit) break;
        double elapsed = chrono::duration<double>(
            chrono::steady_clock::now() - start).count();
        if (elapsed >= seconds_limit) break;

        vector<char> in_base(n + 1, 0);
        for (int leaf : base) in_base[leaf] = 1;

        vector<pair<pair<int, int>, int>> ranked_options;
        for (int leaf = 1; leaf <= n; leaf++) {
            if (in_base[leaf]) continue;
            int conflicts = 0;
            int total_distance = 0;
            for (int existing : base) {
                int distance = max(h56_leaf_distance(T1, leaf, existing),
                                   h56_leaf_distance(T2, leaf, existing));
                total_distance += distance;
                if (distance > distance_limit) conflicts++;
            }
            if (conflicts > max_remove) continue;
            ranked_options.push_back({{conflicts, total_distance}, leaf});
        }
        sort(ranked_options.begin(), ranked_options.end());
        vector<int> options;
        for (const auto &entry : ranked_options) {
            options.push_back(entry.second);
            if ((int)options.size() >= option_limit) break;
        }

        int base_added = 0;
        int base_tried = 0;
        vector<vector<int>> base_candidates;
        unordered_set<string> base_candidate_seen;
        vector<int> removal_items = base;
        sort(removal_items.begin(), removal_items.end(),
             [&](int a, int b) {
                 int total_a = 0;
                 int total_b = 0;
                 int max_a = 0;
                 int max_b = 0;
                 for (int other : base) {
                     if (other != a) {
                         int distance = max(h56_leaf_distance(T1, a, other),
                                            h56_leaf_distance(T2, a, other));
                         total_a += distance;
                         max_a = max(max_a, distance);
                     }
                     if (other != b) {
                         int distance = max(h56_leaf_distance(T1, b, other),
                                            h56_leaf_distance(T2, b, other));
                         total_b += distance;
                         max_b = max(max_b, distance);
                     }
                 }
                 if (total_a != total_b) return total_a > total_b;
                 if (max_a != max_b) return max_a > max_b;
                 return a < b;
             });
        for (int remove_count = 1;
             remove_count <= max_remove &&
             base_tried < per_base_limit * 24;
             remove_count++) {
            vector<int> removed;
            bool keep_base = h56_choose_rec(
                removal_items, remove_count, 0, removed,
                [&](const vector<int> &remove_set) {
                    vector<char> removed_mark(n + 1, 0);
                    for (int leaf : remove_set) removed_mark[leaf] = 1;
                    vector<int> kept;
                    kept.reserve(base.size());
                    for (int leaf : base) {
                        if (!removed_mark[leaf]) kept.push_back(leaf);
                    }
                    vector<int> compatible_options;
                    for (int option : options) {
                        if (h56_adjacent_to_all_leaves(adjacency, option,
                                                       kept)) {
                            compatible_options.push_back(option);
                        }
                    }
                    vector<vector<int>> add_beam{kept};
                    unordered_set<string> add_seen;
                    add_seen.insert(h45_component_key(kept));
                    for (int add_step = 1;
                         add_step <= remove_count &&
                         base_tried < per_base_limit * 24;
                         add_step++) {
                        vector<vector<int>> next_beam;
                        for (const vector<int> &state : add_beam) {
                            vector<pair<int, int>> add_options;
                            for (int option : compatible_options) {
                                if (binary_search(state.begin(), state.end(),
                                                  option)) {
                                    continue;
                                }
                                if (!h56_adjacent_to_all_leaves(
                                        adjacency, option, state)) {
                                    continue;
                                }
                                vector<int> grown = state;
                                grown.push_back(option);
                                sort(grown.begin(), grown.end());
                                if (!h56_component_valid(grown, T1, T2))
                                    continue;
                                add_options.push_back(
                                    {-h56_component_score(grown, T1, T2),
                                     option});
                            }
                            sort(add_options.begin(), add_options.end());
                            if ((int)add_options.size() > add_option_limit)
                                add_options.resize(add_option_limit);

                            for (const auto &entry : add_options) {
                                vector<int> candidate = state;
                                candidate.push_back(entry.second);
                                sort(candidate.begin(), candidate.end());
                                candidate.erase(unique(candidate.begin(),
                                                       candidate.end()),
                                                candidate.end());
                                string key = h45_component_key(candidate);
                                if (!add_seen.insert(key).second) continue;
                                tried++;
                                base_tried++;
                                if (base_candidate_seen.insert(key).second)
                                    base_candidates.push_back(candidate);
                                next_beam.push_back(std::move(candidate));
                                if ((int)pool.size() >= max_pool ||
                                    accepted >= add_limit ||
                                    base_tried >= per_base_limit * 24) {
                                    break;
                                }
                            }
                            if ((int)pool.size() >= max_pool ||
                                accepted >= add_limit ||
                                base_tried >= per_base_limit * 24) {
                                break;
                            }
                        }
                        sort(next_beam.begin(), next_beam.end(),
                             [&](const vector<int> &a,
                                 const vector<int> &b) {
                                 int sa = h56_component_score(a, T1, T2);
                                 int sb = h56_component_score(b, T1, T2);
                                 if (sa != sb) return sa > sb;
                                 if (a.size() != b.size())
                                     return a.size() > b.size();
                                 return a < b;
                             });
                        next_beam.erase(unique(next_beam.begin(),
                                               next_beam.end()),
                                        next_beam.end());
                        if ((int)next_beam.size() > add_beam_width)
                            next_beam.resize(add_beam_width);
                        add_beam = std::move(next_beam);
                        if (add_beam.empty()) break;
                    }
                    if ((int)pool.size() >= max_pool ||
                        accepted >= add_limit ||
                        base_tried >= per_base_limit * 24) {
                        return false;
                    }
                    return true;
                });
            if (!keep_base) break;
        }

        sort(base_candidates.begin(), base_candidates.end(),
             [&](const vector<int> &a, const vector<int> &b) {
                 int da = abs((int)a.size() - preferred_size);
                 int db = abs((int)b.size() - preferred_size);
                 if (da != db) return da < db;
                 int sa = h56_component_score(a, T1, T2);
                 int sb = h56_component_score(b, T1, T2);
                 if (sa != sb) return sa > sb;
                 return a < b;
             });
        base_candidates.erase(unique(base_candidates.begin(),
                                     base_candidates.end()),
                              base_candidates.end());
        for (const vector<int> &candidate : base_candidates) {
            if (base_added >= per_base_limit ||
                accepted >= add_limit ||
                (int)pool.size() >= max_pool) {
                break;
            }
            string key = h45_component_key(candidate);
            if (!generated.insert(key).second) continue;
            size_t old_size = pool.size();
            if (h56_add_candidate(pool, candidate, T1, T2, max_pool)) {
                int delta = pool.size() != old_size;
                accepted += delta;
                base_added += delta;
            }
        }
    }

    cerr << "h56 swap_refine before=" << before
         << " after=" << pool.size()
         << " bases=" << bases.size()
         << " accepted=" << accepted
         << " tried=" << tried << '\n';
}

static void h56_add_local_edit_components(
    unordered_map<string, vector<int>> &pool,
    const Tree &T1, const Tree &T2,
    int max_pool) {
    const int n = T1.n;
    const int max_n = h56_env_int("H56_EDIT_MAX_N", 1500);
    if (n > max_n) {
        cerr << "h56 local_edit_skip leaves=" << n
             << " max_n=" << max_n << '\n';
        return;
    }

    const int distance_limit = h56_env_int("H56_EDIT_DISTANCE", 24);
    const int max_remove = h56_env_int("H56_EDIT_MAX_REMOVE", 3);
    const int max_add = h56_env_int("H56_EDIT_MAX_ADD", 3);
    const int min_base_size = h56_env_int("H56_EDIT_MIN_BASE_SIZE", 4);
    const int max_base_size = h56_env_int("H56_EDIT_MAX_BASE_SIZE", 13);
    const int extra_limit = h56_env_int("H56_EDIT_EXTRA_LIMIT", 48);
    const int extra_conflicts =
        h56_env_int("H56_EDIT_EXTRA_CONFLICTS", max_remove);
    const int base_limit = h56_env_int("H56_EDIT_BASE_LIMIT", 2600);
    const int add_limit = h56_env_int("H56_EDIT_ADD_LIMIT", 2200);
    const int per_base_limit = h56_env_int("H56_EDIT_PER_BASE_LIMIT", 32);
    const int per_remove_limit =
        h56_env_int("H56_EDIT_PER_REMOVE_LIMIT", 1);
    const int per_base_try_limit =
        h56_env_int("H56_EDIT_PER_BASE_TRY_LIMIT", 1500);
    const int remove_set_limit =
        h56_env_int("H56_EDIT_REMOVE_SET_LIMIT", 4);
    const int beam_width = h56_env_int("H56_EDIT_BEAM_WIDTH", 96);
    const int beam_option_limit =
        h56_env_int("H56_EDIT_BEAM_OPTION_LIMIT", 48);
    const int exact_base_limit =
        h56_env_int("H56_EDIT_EXACT_BASES", 700);
    const int pair_beam = h56_env_int("H56_EDIT_PAIR_BEAM", 256);
    const int preferred_size =
        h56_env_int("H56_EDIT_PREFERRED_SIZE", 10);
    const double seconds_limit =
        h56_env_double("H56_EDIT_SECONDS", 20.0);
    auto start = chrono::steady_clock::now();

    vector<vector<int>> bases;
    bases.reserve(pool.size());
    for (const auto &entry : pool) {
        const vector<int> &component = entry.second;
        if ((int)component.size() < min_base_size ||
            (int)component.size() > max_base_size) {
            continue;
        }
        bases.push_back(component);
    }
    sort(bases.begin(), bases.end());
    bases.erase(unique(bases.begin(), bases.end()), bases.end());
    sort(bases.begin(), bases.end(),
         [&](const vector<int> &a, const vector<int> &b) {
             int da = abs((int)a.size() - preferred_size);
             int db = abs((int)b.size() - preferred_size);
             if (da != db) return da < db;
             pair<int, int> ba = h56_component_badness(a, T1, T2);
             pair<int, int> bb = h56_component_badness(b, T1, T2);
             if (ba.first != bb.first) return ba.first > bb.first;
             int sa = h56_component_score(a, T1, T2);
             int sb = h56_component_score(b, T1, T2);
             if (sa != sb) return sa > sb;
             return a < b;
         });
    if ((int)bases.size() > base_limit) bases.resize(base_limit);

    int before = (int)pool.size();
    int accepted = 0;
    long long tried = 0;
    unordered_set<string> generated;

    for (int base_index = 0; base_index < (int)bases.size(); base_index++) {
        const vector<int> &base = bases[base_index];
        if ((int)pool.size() >= max_pool || accepted >= add_limit) break;
        double elapsed = chrono::duration<double>(
            chrono::steady_clock::now() - start).count();
        if (elapsed >= seconds_limit) break;

        vector<char> in_base(n + 1, 0);
        for (int leaf : base) in_base[leaf] = 1;
        vector<int> removal_items = base;
        sort(removal_items.begin(), removal_items.end(),
             [&](int a, int b) {
                 int total_a = 0;
                 int total_b = 0;
                 int max_a = 0;
                 int max_b = 0;
                 for (int other : base) {
                     if (other != a) {
                         int distance = max(h56_leaf_distance(T1, a, other),
                                            h56_leaf_distance(T2, a, other));
                         total_a += distance;
                         max_a = max(max_a, distance);
                     }
                     if (other != b) {
                         int distance = max(h56_leaf_distance(T1, b, other),
                                            h56_leaf_distance(T2, b, other));
                         total_b += distance;
                         max_b = max(max_b, distance);
                     }
                 }
                 if (total_a != total_b) return total_a > total_b;
                 if (max_a != max_b) return max_a > max_b;
                 return a < b;
             });

        vector<vector<int>> base_candidates;
        unordered_set<string> base_seen;
        int base_tried = 0;

        for (int remove_count = max_remove;
             remove_count >= 1 &&
             base_tried < per_base_try_limit;
             remove_count--) {
            vector<int> removed;
            int remove_sets_seen = 0;
            bool keep_base = h56_choose_rec(
                removal_items, remove_count, 0, removed,
                [&](const vector<int> &remove_set) {
                    remove_sets_seen++;
                    if (remove_sets_seen > remove_set_limit)
                        return false;
                    vector<char> removed_mark(n + 1, 0);
                    for (int leaf : remove_set) removed_mark[leaf] = 1;
                    vector<int> kept;
                    kept.reserve(base.size());
                    for (int leaf : base)
                        if (!removed_mark[leaf]) kept.push_back(leaf);
                    if (kept.empty()) return true;

                    vector<pair<int, int>> ranked_extras;
                    for (int leaf = 1; leaf <= n; leaf++) {
                        if (in_base[leaf]) continue;
                        int conflict_count = 0;
                        int total_distance = 0;
                        int worst_distance = 0;
                        for (int existing : kept) {
                            int distance = max(
                                h56_leaf_distance(T1, leaf, existing),
                                h56_leaf_distance(T2, leaf, existing));
                            total_distance += distance;
                            worst_distance = max(worst_distance, distance);
                            if (distance > distance_limit) conflict_count++;
                        }
                        if (conflict_count > extra_conflicts) continue;
                        vector<int> one = kept;
                        one.push_back(leaf);
                        sort(one.begin(), one.end());
                        int score = h56_component_score(one, T1, T2) -
                                    90000 * conflict_count -
                                    25 * total_distance -
                                    80 * worst_distance;
                        ranked_extras.push_back({-score, leaf});
                    }
                    sort(ranked_extras.begin(), ranked_extras.end());
                    vector<int> extras;
                    for (const auto &entry : ranked_extras) {
                        extras.push_back(entry.second);
                        if ((int)extras.size() >= extra_limit) break;
                    }

                    vector<vector<int>> remove_candidates;
                    auto consider_fast_complete = [&](vector<int> candidate) {
                        sort(candidate.begin(), candidate.end());
                        candidate.erase(unique(candidate.begin(),
                                               candidate.end()),
                                        candidate.end());
                        base_tried++;
                        tried++;
                        if (h56_component_valid(candidate, T1, T2))
                            remove_candidates.push_back(std::move(candidate));
                        return base_tried < per_base_try_limit;
                    };

                    if (base_index < exact_base_limit &&
                        max_add >= 3 && (int)extras.size() >= 3 &&
                        base_tried < per_base_try_limit) {
                        int top_count = min(24, (int)extras.size());
                        for (int i = 0;
                             i < top_count &&
                             base_tried < per_base_try_limit;
                             i++) {
                            for (int j = 0;
                                 j + 1 < (int)extras.size() &&
                                 base_tried < per_base_try_limit;
                                 j++) {
                                if (j == i || j + 1 == i) continue;
                                vector<int> candidate = kept;
                                candidate.push_back(extras[i]);
                                candidate.push_back(extras[j]);
                                candidate.push_back(extras[j + 1]);
                                if (!consider_fast_complete(
                                        std::move(candidate))) {
                                    break;
                                }
                            }
                        }
                    }

                    vector<vector<int>> beam{kept};
                    unordered_set<string> beam_seen;
                    beam_seen.insert(h45_component_key(kept));
                    for (int add_step = 1;
                         add_step <= max_add &&
                         base_tried < per_base_try_limit &&
                         !beam.empty();
                         add_step++) {
                        vector<vector<int>> next_beam;
                        for (const vector<int> &state : beam) {
                            vector<pair<int, int>> add_options;
                            for (int option : extras) {
                                if (binary_search(state.begin(), state.end(),
                                                  option)) {
                                    continue;
                                }
                                vector<int> grown = state;
                                grown.push_back(option);
                                sort(grown.begin(), grown.end());
                                int score = h56_component_score(grown, T1, T2);
                                add_options.push_back({-score, option});
                            }
                            sort(add_options.begin(), add_options.end());
                            if ((int)add_options.size() > beam_option_limit)
                                add_options.resize(beam_option_limit);

                            for (const auto &entry : add_options) {
                                vector<int> candidate = state;
                                candidate.push_back(entry.second);
                                sort(candidate.begin(), candidate.end());
                                candidate.erase(unique(candidate.begin(),
                                                       candidate.end()),
                                                candidate.end());
                                string key = h45_component_key(candidate);
                                if (!beam_seen.insert(key).second) continue;
                                base_tried++;
                                tried++;
                                if (h56_component_valid(candidate, T1, T2))
                                    remove_candidates.push_back(candidate);
                                next_beam.push_back(std::move(candidate));
                                if (base_tried >= per_base_try_limit) break;
                            }
                            if (base_tried >= per_base_try_limit) break;
                        }
                        sort(next_beam.begin(), next_beam.end(),
                             [&](const vector<int> &a,
                                 const vector<int> &b) {
                                 int sa = h56_component_score(a, T1, T2);
                                 int sb = h56_component_score(b, T1, T2);
                                 if (sa != sb) return sa > sb;
                                 if (a.size() != b.size())
                                     return a.size() > b.size();
                                 return a < b;
                             });
                        next_beam.erase(unique(next_beam.begin(),
                                               next_beam.end()),
                                        next_beam.end());
                        if ((int)next_beam.size() > beam_width)
                            next_beam.resize(beam_width);
                        beam = std::move(next_beam);
                    }

                    if (base_index < exact_base_limit &&
                        base_tried < per_base_try_limit) {
                        auto consider_complete = [&](vector<int> candidate) {
                            sort(candidate.begin(), candidate.end());
                            candidate.erase(unique(candidate.begin(),
                                                   candidate.end()),
                                            candidate.end());
                            string key = h45_component_key(candidate);
                            base_tried++;
                            tried++;
                            if (h56_component_valid(candidate, T1, T2))
                                remove_candidates.push_back(std::move(candidate));
                            return base_tried < per_base_try_limit;
                        };

                        if (max_add >= 3 && (int)extras.size() >= 3) {
                            int top_count = min(24, (int)extras.size());
                            for (int i = 0;
                                 i < top_count &&
                                 base_tried < per_base_try_limit;
                                 i++) {
                                for (int j = 0;
                                     j + 1 < (int)extras.size() &&
                                     base_tried < per_base_try_limit;
                                     j++) {
                                    if (j == i || j + 1 == i) continue;
                                    vector<int> candidate = kept;
                                    candidate.push_back(extras[i]);
                                    candidate.push_back(extras[j]);
                                    candidate.push_back(extras[j + 1]);
                                    if (!consider_complete(
                                            std::move(candidate))) {
                                        break;
                                    }
                                }
                            }
                        }

                        if (max_add >= 2) {
                            vector<pair<int, pair<int, int>>> ranked_pairs;
                            for (int i = 0; i < (int)extras.size(); i++) {
                                for (int j = i + 1;
                                     j < (int)extras.size(); j++) {
                                    vector<int> candidate = kept;
                                    candidate.push_back(extras[i]);
                                    candidate.push_back(extras[j]);
                                    sort(candidate.begin(), candidate.end());
                                    int score = h56_component_score(candidate,
                                                                    T1, T2);
                                    ranked_pairs.push_back(
                                        {-score, {extras[i], extras[j]}});
                                }
                            }
                            sort(ranked_pairs.begin(), ranked_pairs.end());

                            int pair_limit = min(pair_beam,
                                                 (int)ranked_pairs.size());
                            for (int pi = 0;
                                 pi < pair_limit &&
                                 base_tried < per_base_try_limit;
                                 pi++) {
                                vector<int> candidate = kept;
                                candidate.push_back(ranked_pairs[pi].second.first);
                                candidate.push_back(ranked_pairs[pi].second.second);
                                if (!consider_complete(std::move(candidate)))
                                    break;
                            }

                            if (max_add >= 3) {
                                for (int pi = 0;
                                     pi < pair_limit &&
                                     base_tried < per_base_try_limit;
                                     pi++) {
                                    int a = ranked_pairs[pi].second.first;
                                    int b = ranked_pairs[pi].second.second;
                                    vector<pair<int, int>> third_options;
                                    for (int leaf : extras) {
                                        if (leaf == a || leaf == b) continue;
                                        vector<int> candidate = kept;
                                        candidate.push_back(a);
                                        candidate.push_back(b);
                                        candidate.push_back(leaf);
                                        sort(candidate.begin(),
                                             candidate.end());
                                        int score = h56_component_score(
                                            candidate, T1, T2);
                                        third_options.push_back({-score, leaf});
                                    }
                                    sort(third_options.begin(),
                                         third_options.end());
                                    int third_limit =
                                        min(beam_option_limit,
                                            (int)third_options.size());
                                    for (int ti = 0;
                                         ti < third_limit &&
                                         base_tried < per_base_try_limit;
                                         ti++) {
                                        vector<int> candidate = kept;
                                        candidate.push_back(a);
                                        candidate.push_back(b);
                                        candidate.push_back(
                                            third_options[ti].second);
                                        if (!consider_complete(
                                                std::move(candidate))) {
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    h56_keep_top_edits_for_remove_set(
                        base_candidates, base_seen,
                        std::move(remove_candidates), T1, T2,
                        per_remove_limit);

                    return base_tried < per_base_try_limit &&
                           (int)pool.size() < max_pool &&
                           accepted < add_limit;
                });
            if (!keep_base) break;
        }

        sort(base_candidates.begin(), base_candidates.end(),
             [&](const vector<int> &a, const vector<int> &b) {
                 int da = abs((int)a.size() - preferred_size);
                 int db = abs((int)b.size() - preferred_size);
                 if (da != db) return da < db;
                 int sa = h56_component_score(a, T1, T2);
                 int sb = h56_component_score(b, T1, T2);
                 if (sa != sb) return sa > sb;
                 return a < b;
             });
        base_candidates.erase(unique(base_candidates.begin(),
                                     base_candidates.end()),
                              base_candidates.end());

        int base_added = 0;
        for (const vector<int> &candidate : base_candidates) {
            if (base_added >= per_base_limit ||
                accepted >= add_limit ||
                (int)pool.size() >= max_pool) {
                break;
            }
            string key = h45_component_key(candidate);
            if (!generated.insert(key).second) continue;
            size_t old_size = pool.size();
            if (h56_add_candidate(pool, candidate, T1, T2, max_pool)) {
                int delta = pool.size() != old_size;
                accepted += delta;
                base_added += delta;
            }
        }
    }

    cerr << "h56 local_edit before=" << before
         << " after=" << pool.size()
         << " bases=" << bases.size()
         << " accepted=" << accepted
         << " tried=" << tried << '\n';
}

static void h56_add_external_local_triples(
    unordered_map<string, vector<int>> &pool,
    const Tree &T1, const Tree &T2,
    int max_pool,
    const chrono::steady_clock::time_point &global_start,
    double pool_cutoff_seconds) {
    const int n = T1.n;
    const int distance_limit =
        h56_env_int("H56_EXTERNAL_TRIPLE_DISTANCE", 8);
    const double seconds_limit =
        h56_env_double("H56_EXTERNAL_TRIPLE_SECONDS", 18.0);
    auto local_start = chrono::steady_clock::now();
    auto still_time = [&]() {
        return h56_elapsed_seconds(global_start) < pool_cutoff_seconds &&
               h56_elapsed_seconds(local_start) < seconds_limit;
    };

    int generated_pairs = 0;
    int generated_triples = 0;
    vector<vector<int>> local_after(n + 1);
    unordered_set<uint64_t> local_pair;
    local_pair.reserve((size_t)n * 16);

    for (int a = 1; a <= n && still_time(); a++) {
        for (int b = a + 1; b <= n; b++) {
            int distance = max(h56_leaf_distance(T1, a, b),
                               h56_leaf_distance(T2, a, b));
            if (distance > distance_limit) continue;
            local_after[a].push_back(b);
            local_pair.insert(representative_pair_key(a, b));
            if (h56_add_candidate(pool, vector<int>{a, b}, T1, T2,
                                  max_pool)) {
                generated_pairs++;
            }
            if ((int)pool.size() >= max_pool || !still_time()) {
                cerr << "h56 external_pairs_triples pairs="
                     << generated_pairs
                     << " triples=" << generated_triples
                     << " pool=" << pool.size() << '\n';
                return;
            }
        }
    }

    for (int a = 1; a <= n && still_time(); a++) {
        const vector<int> &neighbors = local_after[a];
        for (int left = 0; left < (int)neighbors.size(); left++) {
            int b = neighbors[left];
            for (int right = left + 1; right < (int)neighbors.size();
                 right++) {
                int c = neighbors[right];
                if (!local_pair.count(representative_pair_key(b, c)))
                    continue;
                if (h56_add_candidate(pool, vector<int>{a, b, c}, T1, T2,
                                      max_pool)) {
                    generated_triples++;
                }
                if ((int)pool.size() >= max_pool || !still_time()) {
                    cerr << "h56 external_pairs_triples pairs="
                         << generated_pairs
                         << " triples=" << generated_triples
                         << " pool=" << pool.size() << '\n';
                    return;
                }
            }
        }
    }
    cerr << "h56 external_pairs_triples pairs=" << generated_pairs
         << " triples=" << generated_triples
         << " pool=" << pool.size() << '\n';
}

static void h56_add_external_component_unions(
    unordered_map<string, vector<int>> &pool,
    const Tree &T1, const Tree &T2,
    int max_pool,
    const chrono::steady_clock::time_point &global_start,
    double pool_cutoff_seconds) {
    const int distance_limit =
        h56_env_int("H56_EXTERNAL_MERGE_DISTANCE", 12);
    const int max_leaves =
        h56_env_int("H56_EXTERNAL_MERGE_MAX_LEAVES", 40);
    const int component_limit =
        h56_env_int("H56_EXTERNAL_COMPONENT_LIMIT", 18000);
    const double seconds_limit =
        h56_env_double("H56_EXTERNAL_UNION_SECONDS", 28.0);
    auto local_start = chrono::steady_clock::now();
    auto still_time = [&]() {
        return h56_elapsed_seconds(global_start) < pool_cutoff_seconds &&
               h56_elapsed_seconds(local_start) < seconds_limit;
    };

    vector<vector<int>> components =
        h56_valid_pool_components(pool, T1, T2);
    if ((int)components.size() > component_limit)
        components.resize(component_limit);

    vector<int> root1(components.size());
    vector<int> root2(components.size());
    for (int i = 0; i < (int)components.size(); i++) {
        root1[i] = lca_of_leaf_set(T1, components[i]);
        root2[i] = lca_of_leaf_set(T2, components[i]);
    }

    int checked = 0;
    int generated = 0;
    for (int i = 0; i < (int)components.size() && still_time(); i++) {
        if ((int)components[i].size() >= max_leaves) continue;
        for (int j = i + 1; j < (int)components.size(); j++) {
            int combined_size =
                (int)components[i].size() + (int)components[j].size();
            if (combined_size > max_leaves) break;
            if (!h56_components_disjoint(components[i], components[j]))
                continue;
            if (h56_node_distance(T1, root1[i], root1[j]) >
                distance_limit) {
                continue;
            }
            if (h56_node_distance(T2, root2[i], root2[j]) >
                distance_limit) {
                continue;
            }
            checked++;
            vector<int> leaves = components[i];
            leaves.insert(leaves.end(),
                          components[j].begin(), components[j].end());
            if (h56_add_candidate(pool, std::move(leaves), T1, T2,
                                  max_pool)) {
                generated++;
            }
            if ((int)pool.size() >= max_pool || !still_time()) {
                cerr << "h56 external_component_unions generated="
                     << generated
                     << " checked=" << checked
                     << " components=" << components.size()
                     << " pool=" << pool.size() << '\n';
                return;
            }
        }
    }
    cerr << "h56 external_component_unions generated=" << generated
         << " checked=" << checked
         << " components=" << components.size()
         << " pool=" << pool.size() << '\n';
}

static void h56_add_external_like_pool(
    unordered_map<string, vector<int>> &pool,
    const Tree &T1, const Tree &T2,
    int max_pool,
    const chrono::steady_clock::time_point &global_start,
    double pool_cutoff_seconds) {
    const int max_n = h56_env_int("H56_EXTERNAL_MAX_N", 1800);
    if (T1.n > max_n) {
        cerr << "h56 external_pool_skip leaves=" << T1.n
             << " max_n=" << max_n << '\n';
        return;
    }
    int before = (int)pool.size();
    for (int leaf = 1; leaf <= T1.n && (int)pool.size() < max_pool; leaf++)
        h56_add_candidate(pool, vector<int>{leaf}, T1, T2, max_pool);
    h56_add_external_local_triples(pool, T1, T2, max_pool, global_start,
                                   pool_cutoff_seconds);
    if ((int)pool.size() < max_pool &&
        h56_elapsed_seconds(global_start) < pool_cutoff_seconds) {
        h56_add_external_component_unions(pool, T1, T2, max_pool,
                                          global_start,
                                          pool_cutoff_seconds);
    }
    cerr << "h56 external_pool before=" << before
         << " after=" << pool.size()
         << " cap=" << max_pool
         << " elapsed=" << h56_elapsed_seconds(global_start) << '\n';
}

static int h56_leaf_to_component_distance(const Tree &T,
                                          int leaf,
                                          const vector<int> &component) {
    int leaf_node = T.leaf_node[leaf];
    int best = 1 << 28;
    for (int other : component) {
        int other_node = T.leaf_node[other];
        best = min(best, h56_node_distance(T, leaf_node, other_node));
    }
    return best;
}

static vector<int> h56_rank_removal_items(const vector<int> &base,
                                          const Tree &T1,
                                          const Tree &T2,
                                          int limit) {
    vector<pair<pair<int, int>, int>> ranked;
    ranked.reserve(base.size());
    for (int leaf : base) {
        int total = 0;
        int worst = 0;
        for (int other : base) {
            if (leaf == other) continue;
            int distance = max(h56_leaf_distance(T1, leaf, other),
                               h56_leaf_distance(T2, leaf, other));
            total += distance;
            worst = max(worst, distance);
        }
        ranked.push_back({{-total, -worst}, leaf});
    }
    sort(ranked.begin(), ranked.end());
    vector<int> out;
    for (const auto &entry : ranked) {
        out.push_back(entry.second);
        if ((int)out.size() >= limit) break;
    }
    sort(out.begin(), out.end());
    return out;
}

static void h56_add_region_complement_components(
    unordered_map<string, vector<int>> &pool,
    const vector<H56Region> &regions,
    const Tree &T1,
    const Tree &T2,
    int max_pool,
    const chrono::steady_clock::time_point &global_start,
    double pool_cutoff_seconds) {
    const int region_limit =
        h56_env_int("H56_COMPLEMENT_REGION_LIMIT", 90);
    const int min_leaves =
        h56_env_int("H56_COMPLEMENT_MIN_LEAVES", 8);
    const int max_leaves =
        h56_env_int("H56_COMPLEMENT_MAX_LEAVES", 42);
    const int removal_limit =
        h56_env_int("H56_COMPLEMENT_MAX_REMOVE", 6);
    const int removal_item_limit =
        h56_env_int("H56_COMPLEMENT_REMOVAL_ITEM_LIMIT", 16);
    const int per_region_limit =
        h56_env_int("H56_COMPLEMENT_PER_REGION", 80);
    const int add_limit =
        h56_env_int("H56_COMPLEMENT_ADD_LIMIT", 1400);
    const double seconds_limit =
        h56_env_double("H56_COMPLEMENT_SECONDS", 4.0);

    auto local_start = chrono::steady_clock::now();
    auto still_time = [&]() {
        return h56_elapsed_seconds(global_start) < pool_cutoff_seconds &&
               h56_elapsed_seconds(local_start) < seconds_limit;
    };

    int before = (int)pool.size();
    int used_regions = 0;
    int accepted = 0;
    long long tried = 0;
    unordered_set<string> seen;
    seen.reserve((size_t)add_limit * 4 + 16);

    for (const H56Region &region : regions) {
        if ((int)pool.size() >= max_pool || accepted >= add_limit ||
            !still_time()) {
            break;
        }
        vector<int> leaves = region.leaves;
        sort(leaves.begin(), leaves.end());
        leaves.erase(unique(leaves.begin(), leaves.end()), leaves.end());
        if ((int)leaves.size() < min_leaves ||
            (int)leaves.size() > max_leaves) {
            continue;
        }
        used_regions++;
        vector<int> removal_items = h56_rank_removal_items(
            leaves, T1, T2, min(removal_item_limit, (int)leaves.size()));
        int local_before = accepted;

        for (int remove_count = 1;
             remove_count <= removal_limit &&
             remove_count < (int)leaves.size() &&
             accepted - local_before < per_region_limit &&
             (int)pool.size() < max_pool &&
             accepted < add_limit &&
             still_time();
             remove_count++) {
            vector<int> remove_set;
            bool keep_going = h56_choose_rec(
                removal_items, remove_count, 0, remove_set,
                [&](const vector<int> &chosen) {
                    if ((int)pool.size() >= max_pool ||
                        accepted - local_before >= per_region_limit ||
                        accepted >= add_limit ||
                        !still_time()) {
                        return false;
                    }
                    tried++;
                    vector<char> removed(T1.n + 1, 0);
                    for (int leaf : chosen) removed[leaf] = 1;
                    vector<int> candidate;
                    candidate.reserve(leaves.size() - chosen.size());
                    for (int leaf : leaves)
                        if (!removed[leaf]) candidate.push_back(leaf);
                    if ((int)candidate.size() < 2) return true;
                    string key = h45_component_key(candidate);
                    if (!seen.insert(key).second) return true;
                    size_t before_pool = pool.size();
                    h56_add_candidate(pool, std::move(candidate), T1, T2,
                                      max_pool);
                    if (pool.size() != before_pool) accepted++;
                    return true;
                });
            if (!keep_going) break;
        }
        if (used_regions >= region_limit) break;
    }

    cerr << "h56 region_complements before=" << before
         << " after=" << pool.size()
         << " regions=" << used_regions
         << " accepted=" << accepted
         << " tried=" << tried
         << " elapsed=" << h56_elapsed_seconds(global_start) << '\n';
}

static void h56_add_completion_growth_components(
    unordered_map<string, vector<int>> &pool,
    const Tree &T1, const Tree &T2,
    int max_pool,
    const chrono::steady_clock::time_point &global_start,
    double pool_cutoff_seconds) {
    if (!h56_env_int("H56_ENABLE_COMPLETION",
                     H56_ENABLE_COMPLETION_VALUE)) {
        return;
    }
    const int n = T1.n;
    const int max_n = h56_env_int("H56_COMPLETION_MAX_N", 1800);
    if (n > max_n) {
        cerr << "h56 completion_skip leaves=" << n
             << " max_n=" << max_n << '\n';
        return;
    }

    const int min_base_size =
        h56_env_int("H56_COMPLETION_MIN_BASE_SIZE", 3);
    const int max_base_size =
        h56_env_int("H56_COMPLETION_MAX_BASE_SIZE", 34);
    const int max_component_size =
        h56_env_int("H56_COMPLETION_MAX_COMPONENT_SIZE", 40);
    const int max_add = h56_env_int("H56_COMPLETION_MAX_ADD", 12);
    const int max_remove = h56_env_int("H56_COMPLETION_MAX_REMOVE", 2);
    const int removal_item_limit =
        h56_env_int("H56_COMPLETION_REMOVAL_ITEM_LIMIT", 10);
    const int base_limit =
        h56_env_int("H56_COMPLETION_BASE_LIMIT", 1800);
    const int option_limit =
        h56_env_int("H56_COMPLETION_OPTION_LIMIT", 96);
    const int beam_width =
        h56_env_int("H56_COMPLETION_BEAM_WIDTH", 24);
    const int add_limit =
        h56_env_int("H56_COMPLETION_ADD_LIMIT", 2800);
    const double seconds_limit =
        h56_env_double("H56_COMPLETION_SECONDS", 18.0);
    auto local_start = chrono::steady_clock::now();
    auto still_time = [&]() {
        return h56_elapsed_seconds(global_start) < pool_cutoff_seconds &&
               h56_elapsed_seconds(local_start) < seconds_limit;
    };

    vector<vector<int>> bases = h56_valid_pool_components(pool, T1, T2);
    vector<vector<int>> filtered;
    filtered.reserve(bases.size());
    for (vector<int> &base : bases) {
        if ((int)base.size() < min_base_size ||
            (int)base.size() > max_base_size) {
            continue;
        }
        filtered.push_back(std::move(base));
    }
    sort(filtered.begin(), filtered.end(),
         [&](const vector<int> &a, const vector<int> &b) {
             int da = abs((int)a.size() -
                          h56_env_int("H56_COMPLETION_PREFERRED_SIZE", 24));
             int db = abs((int)b.size() -
                          h56_env_int("H56_COMPLETION_PREFERRED_SIZE", 24));
             if (da != db) return da < db;
             if (a.size() != b.size()) return a.size() > b.size();
             int sa = h56_component_score(a, T1, T2);
             int sb = h56_component_score(b, T1, T2);
             if (sa != sb) return sa > sb;
             return a < b;
         });
    if ((int)filtered.size() > base_limit)
        filtered.resize(base_limit);

    int before = (int)pool.size();
    int accepted = 0;
    long long tried = 0;
    unordered_set<string> generated;

    auto add_seed_variants = [&](const vector<int> &base) {
        vector<vector<int>> seeds{base};
        vector<int> removal_items = h56_rank_removal_items(
            base, T1, T2, removal_item_limit);
        for (int remove_count = 1; remove_count <= max_remove;
             remove_count++) {
            vector<int> removed;
            h56_choose_rec(
                removal_items, remove_count, 0, removed,
                [&](const vector<int> &remove_set) {
                    vector<char> mark(n + 1, 0);
                    for (int leaf : remove_set) mark[leaf] = 1;
                    vector<int> seed;
                    seed.reserve(base.size() - remove_set.size());
                    for (int leaf : base)
                        if (!mark[leaf]) seed.push_back(leaf);
                    if ((int)seed.size() >= min_base_size &&
                        h56_component_valid(seed, T1, T2)) {
                        seeds.push_back(std::move(seed));
                    }
                    return (int)seeds.size() < 16;
                });
            if ((int)seeds.size() >= 16) break;
        }
        return seeds;
    };

    for (const vector<int> &base : filtered) {
        if ((int)pool.size() >= max_pool || accepted >= add_limit ||
            !still_time()) {
            break;
        }

        vector<vector<int>> seeds = add_seed_variants(base);
        for (const vector<int> &seed : seeds) {
            if ((int)pool.size() >= max_pool || accepted >= add_limit ||
                !still_time()) {
                break;
            }
            vector<char> in_seed(n + 1, 0);
            for (int leaf : seed) in_seed[leaf] = 1;

            vector<pair<pair<int, int>, int>> ranked_options;
            ranked_options.reserve(n);
            int root1 = lca_of_leaf_set(T1, seed);
            int root2 = lca_of_leaf_set(T2, seed);
            for (int leaf = 1; leaf <= n; leaf++) {
                if (in_seed[leaf]) continue;
                int pair_distance =
                    max(h56_leaf_to_component_distance(T1, leaf, seed),
                        h56_leaf_to_component_distance(T2, leaf, seed));
                int root_distance =
                    h56_node_distance(T1, T1.leaf_node[leaf], root1) +
                    h56_node_distance(T2, T2.leaf_node[leaf], root2);
                ranked_options.push_back(
                    {{pair_distance, root_distance}, leaf});
            }
            sort(ranked_options.begin(), ranked_options.end());
            if ((int)ranked_options.size() > option_limit)
                ranked_options.resize(option_limit);
            vector<int> options;
            options.reserve(ranked_options.size());
            for (const auto &entry : ranked_options)
                options.push_back(entry.second);

            vector<vector<int>> beam{seed};
            unordered_set<string> local_seen;
            local_seen.insert(h45_component_key(seed));
            for (int step = 1; step <= max_add && !beam.empty(); step++) {
                vector<vector<int>> next_beam;
                for (const vector<int> &state : beam) {
                    if ((int)state.size() >= max_component_size) continue;
                    for (int option : options) {
                        if (binary_search(state.begin(), state.end(),
                                          option)) {
                            continue;
                        }
                        vector<int> grown = state;
                        grown.push_back(option);
                        sort(grown.begin(), grown.end());
                        grown.erase(unique(grown.begin(), grown.end()),
                                    grown.end());
                        if ((int)grown.size() > max_component_size)
                            continue;
                        string key = h45_component_key(grown);
                        if (!local_seen.insert(key).second) continue;
                        tried++;
                        if (!h56_component_valid(grown, T1, T2))
                            continue;
                        next_beam.push_back(grown);
                        if (generated.insert(key).second) {
                            size_t old_size = pool.size();
                            if (h56_add_candidate(pool, std::move(grown),
                                                  T1, T2, max_pool)) {
                                accepted += pool.size() != old_size;
                            }
                        }
                        if ((int)pool.size() >= max_pool ||
                            accepted >= add_limit || !still_time()) {
                            break;
                        }
                    }
                    if ((int)pool.size() >= max_pool ||
                        accepted >= add_limit || !still_time()) {
                        break;
                    }
                }
                sort(next_beam.begin(), next_beam.end(),
                     [&](const vector<int> &a, const vector<int> &b) {
                         int sa = h56_component_score(a, T1, T2);
                         int sb = h56_component_score(b, T1, T2);
                         if (sa != sb) return sa > sb;
                         if (a.size() != b.size())
                             return a.size() > b.size();
                         return a < b;
                     });
                next_beam.erase(unique(next_beam.begin(), next_beam.end()),
                                next_beam.end());
                if ((int)next_beam.size() > beam_width)
                    next_beam.resize(beam_width);
                beam = std::move(next_beam);
                if ((int)pool.size() >= max_pool || accepted >= add_limit ||
                    !still_time()) {
                    break;
                }
            }
        }
    }
    cerr << "h56 completion_growth before=" << before
         << " after=" << pool.size()
         << " bases=" << filtered.size()
         << " accepted=" << accepted
         << " tried=" << tried
         << " elapsed=" << h56_elapsed_seconds(global_start) << '\n';
}

struct H56BatchState {
    vector<int> add;
    int score = 0;
};

static void h56_add_batch_completion_components(
    unordered_map<string, vector<int>> &pool,
    const Tree &T1, const Tree &T2,
    int max_pool,
    const chrono::steady_clock::time_point &global_start,
    double pool_cutoff_seconds) {
    if (!h56_env_int("H56_ENABLE_BATCH_COMPLETION",
                     H56_ENABLE_BATCH_COMPLETION_VALUE)) {
        return;
    }
    const int n = T1.n;
    const int max_n = h56_env_int("H56_BATCH_MAX_N", 1800);
    if (n > max_n) {
        cerr << "h56 batch_completion_skip leaves=" << n
             << " max_n=" << max_n << '\n';
        return;
    }

    const int min_base_size = h56_env_int("H56_BATCH_MIN_BASE_SIZE", 4);
    const int max_base_size = h56_env_int("H56_BATCH_MAX_BASE_SIZE", 34);
    const int max_component_size =
        h56_env_int("H56_BATCH_MAX_COMPONENT_SIZE", 42);
    const int max_remove = h56_env_int("H56_BATCH_MAX_REMOVE", 2);
    const int max_add = h56_env_int("H56_BATCH_MAX_ADD", 10);
    const int option_limit = h56_env_int("H56_BATCH_OPTION_LIMIT", 34);
    const int base_limit = h56_env_int("H56_BATCH_BASE_LIMIT", 700);
    const int beam_width = h56_env_int("H56_BATCH_BEAM_WIDTH", 1800);
    const int add_limit = h56_env_int("H56_BATCH_ADD_LIMIT", 1800);
    const int preferred_size = h56_env_int("H56_BATCH_PREFERRED_SIZE", 24);
    const int removal_item_limit =
        h56_env_int("H56_BATCH_REMOVAL_ITEM_LIMIT", 10);
    const double seconds_limit =
        h56_env_double("H56_BATCH_SECONDS", 18.0);
    auto local_start = chrono::steady_clock::now();
    auto still_time = [&]() {
        return h56_elapsed_seconds(global_start) < pool_cutoff_seconds &&
               h56_elapsed_seconds(local_start) < seconds_limit;
    };

    vector<vector<int>> bases = h56_valid_pool_components(pool, T1, T2);
    vector<vector<int>> filtered;
    filtered.reserve(bases.size());
    for (vector<int> &base : bases) {
        if ((int)base.size() < min_base_size ||
            (int)base.size() > max_base_size) {
            continue;
        }
        filtered.push_back(std::move(base));
    }
    sort(filtered.begin(), filtered.end(),
         [&](const vector<int> &a, const vector<int> &b) {
             int da = abs((int)a.size() - preferred_size);
             int db = abs((int)b.size() - preferred_size);
             if (da != db) return da < db;
             if (a.size() != b.size()) return a.size() > b.size();
             pair<int, int> ba = h56_component_badness(a, T1, T2);
             pair<int, int> bb = h56_component_badness(b, T1, T2);
             if (ba.first != bb.first) return ba.first > bb.first;
             int sa = h56_component_score(a, T1, T2);
             int sb = h56_component_score(b, T1, T2);
             if (sa != sb) return sa > sb;
             return a < b;
         });
    if ((int)filtered.size() > base_limit)
        filtered.resize(base_limit);

    int before = (int)pool.size();
    int accepted = 0;
    long long tried = 0;
    unordered_set<string> generated;

    auto add_seed_variants = [&](const vector<int> &base) {
        vector<vector<int>> seeds{base};
        vector<int> removal_items = h56_rank_removal_items(
            base, T1, T2, removal_item_limit);
        for (int remove_count = 1; remove_count <= max_remove;
             remove_count++) {
            vector<int> removed;
            h56_choose_rec(
                removal_items, remove_count, 0, removed,
                [&](const vector<int> &remove_set) {
                    vector<char> mark(n + 1, 0);
                    for (int leaf : remove_set) mark[leaf] = 1;
                    vector<int> seed;
                    seed.reserve(base.size() - remove_set.size());
                    for (int leaf : base)
                        if (!mark[leaf]) seed.push_back(leaf);
                    if ((int)seed.size() >= min_base_size)
                        seeds.push_back(std::move(seed));
                    return (int)seeds.size() < 12;
                });
            if ((int)seeds.size() >= 12) break;
        }
        return seeds;
    };

    auto rank_options = [&](const vector<int> &seed) {
        vector<char> in_seed(n + 1, 0);
        for (int leaf : seed) in_seed[leaf] = 1;
        int root1 = lca_of_leaf_set(T1, seed);
        int root2 = lca_of_leaf_set(T2, seed);
        vector<pair<pair<int, int>, int>> ranked;
        ranked.reserve(n);
        for (int leaf = 1; leaf <= n; leaf++) {
            if (in_seed[leaf]) continue;
            int pair_distance =
                max(h56_leaf_to_component_distance(T1, leaf, seed),
                    h56_leaf_to_component_distance(T2, leaf, seed));
            int root_distance =
                h56_node_distance(T1, T1.leaf_node[leaf], root1) +
                h56_node_distance(T2, T2.leaf_node[leaf], root2);
            ranked.push_back({{pair_distance, root_distance}, leaf});
        }
        sort(ranked.begin(), ranked.end());
        if ((int)ranked.size() > option_limit) ranked.resize(option_limit);
        vector<int> options;
        options.reserve(ranked.size());
        for (const auto &entry : ranked) options.push_back(entry.second);
        sort(options.begin(), options.end());
        return options;
    };

    auto state_score = [&](const vector<int> &seed,
                           const vector<int> &add) {
        int score = 0;
        for (int leaf : add) {
            score += max(h56_leaf_to_component_distance(T1, leaf, seed),
                         h56_leaf_to_component_distance(T2, leaf, seed));
        }
        for (int i = 0; i < (int)add.size(); i++) {
            for (int j = i + 1; j < (int)add.size(); j++) {
                score += max(h56_leaf_distance(T1, add[i], add[j]),
                             h56_leaf_distance(T2, add[i], add[j])) / 2;
            }
        }
        return score;
    };

    for (const vector<int> &base : filtered) {
        if ((int)pool.size() >= max_pool || accepted >= add_limit ||
            !still_time()) {
            break;
        }
        vector<vector<int>> seeds = add_seed_variants(base);
        for (const vector<int> &seed : seeds) {
            if ((int)pool.size() >= max_pool || accepted >= add_limit ||
                !still_time()) {
                break;
            }
            if ((int)seed.size() >= max_component_size) continue;
            vector<int> options = rank_options(seed);
            vector<H56BatchState> beam{{{}, 0}};
            unordered_set<string> local_seen;
            for (int step = 1; step <= max_add && !beam.empty(); step++) {
                vector<H56BatchState> next_beam;
                for (const H56BatchState &state : beam) {
                    int start_index = 0;
                    if (!state.add.empty()) {
                        auto it = upper_bound(options.begin(), options.end(),
                                              state.add.back());
                        start_index = (int)(it - options.begin());
                    }
                    for (int oi = start_index; oi < (int)options.size();
                         oi++) {
                        vector<int> add = state.add;
                        add.push_back(options[oi]);
                        if ((int)seed.size() + (int)add.size() >
                            max_component_size) {
                            continue;
                        }
                        string add_key = h45_component_key(add);
                        if (!local_seen.insert(add_key).second) continue;
                        vector<int> candidate = seed;
                        candidate.insert(candidate.end(),
                                         add.begin(), add.end());
                        sort(candidate.begin(), candidate.end());
                        candidate.erase(unique(candidate.begin(),
                                               candidate.end()),
                                        candidate.end());
                        tried++;
                        if (h56_component_valid(candidate, T1, T2)) {
                            string key = h45_component_key(candidate);
                            if (generated.insert(key).second) {
                                size_t old_size = pool.size();
                                if (h56_add_candidate(pool, candidate, T1, T2,
                                                      max_pool)) {
                                    accepted += pool.size() != old_size;
                                }
                            }
                        }
                        int next_score = state_score(seed, add);
                        next_beam.push_back({std::move(add), next_score});
                        if ((int)pool.size() >= max_pool ||
                            accepted >= add_limit || !still_time()) {
                            break;
                        }
                    }
                    if ((int)pool.size() >= max_pool ||
                        accepted >= add_limit || !still_time()) {
                        break;
                    }
                }
                sort(next_beam.begin(), next_beam.end(),
                     [](const H56BatchState &a,
                        const H56BatchState &b) {
                         if (a.score != b.score)
                             return a.score < b.score;
                         if (a.add.size() != b.add.size())
                             return a.add.size() > b.add.size();
                         return a.add < b.add;
                     });
                next_beam.erase(
                    unique(next_beam.begin(), next_beam.end(),
                           [](const H56BatchState &a,
                              const H56BatchState &b) {
                               return a.add == b.add;
                           }),
                    next_beam.end());
                if ((int)next_beam.size() > beam_width)
                    next_beam.resize(beam_width);
                beam = std::move(next_beam);
            }
        }
    }

    cerr << "h56 batch_completion before=" << before
         << " after=" << pool.size()
         << " bases=" << filtered.size()
         << " accepted=" << accepted
         << " tried=" << tried
         << " elapsed=" << h56_elapsed_seconds(global_start) << '\n';
}

static void h56_enrich_targeted_medium(
    unordered_map<string, vector<int>> &pool,
    const AgreementForest &best,
    const Tree &T1, const Tree &T2,
    const chrono::steady_clock::time_point &global_start,
    double pool_cutoff_seconds) {
    const int medium_min =
        h56_env_int("H56_MEDIUM_INSTANCE_MIN_LEAVES",
                    H56_MEDIUM_INSTANCE_MIN_LEAVES_VALUE);
    const bool medium_instance = T1.n >= medium_min;
    const int default_max_pool =
        medium_instance ? H56_MEDIUM_MAX_POOL_VALUE : H56_MAX_POOL_VALUE;
    const int max_pool = h56_env_int("H56_MAX_POOL", default_max_pool);
    const int default_core_pool =
        medium_instance ? H56_MEDIUM_CORE_POOL_VALUE : max_pool;
    const int core_pool = min(max_pool,
                              h56_env_int("H56_CORE_MAX_POOL",
                                          default_core_pool));
    const int local_pool = h56_env_int("H56_LOCAL_POOL_LIMIT",
                                       H56_LOCAL_POOL_LIMIT_VALUE);
    int before = (int)pool.size();
    vector<H56Region> regions = h56_collect_regions(best, T1, T2);
    auto still_have_pool_time = [&]() {
        return h56_elapsed_seconds(global_start) < pool_cutoff_seconds;
    };
    for (const H56Region &region : regions) {
        if ((int)pool.size() >= local_pool || !still_have_pool_time()) break;
        h56_generate_region_components(region, pool, T1, T2, local_pool);
        h56_generate_clade_intersection_components(region, pool, T1, T2,
                                                   local_pool);
    }
    if ((int)pool.size() < core_pool && still_have_pool_time())
        h56_add_pruned_clade_components(pool, regions, T1, T2, core_pool,
                                        global_start,
                                        pool_cutoff_seconds);
    if ((int)pool.size() < core_pool && still_have_pool_time())
        h56_add_repartition_growth_components(pool, regions, T1, T2,
                                              core_pool, global_start,
                                              pool_cutoff_seconds);
    if ((int)pool.size() < core_pool && still_have_pool_time())
        h56_add_completion_growth_components(pool, T1, T2, core_pool,
                                             global_start,
                                             pool_cutoff_seconds);
    if ((int)pool.size() < core_pool && still_have_pool_time())
        h56_add_batch_completion_components(pool, T1, T2, core_pool,
                                            global_start,
                                            pool_cutoff_seconds);
    if ((int)pool.size() < core_pool && still_have_pool_time())
        h56_add_rank_bucket_beam_components(pool, T1, T2, core_pool);
    if ((int)pool.size() < core_pool && still_have_pool_time()) {
        int reserve_for_swap =
            h56_env_int("H56_SWAP_RESERVED_POOL", 1800);
        int clique_cap = max((int)pool.size(),
                             core_pool - reserve_for_swap);
        h56_add_distance_clique_components(pool, T1, T2, clique_cap);
    }
    if (const char *dump_path = getenv("H56_DUMP_PRE_SWAP_POOL")) {
        ofstream dump(dump_path);
        vector<string> keys;
        keys.reserve(pool.size());
        for (const auto &entry : pool)
            keys.push_back(h45_component_key(entry.second));
        sort(keys.begin(), keys.end());
        keys.erase(unique(keys.begin(), keys.end()), keys.end());
        for (const string &key : keys) dump << key << '\n';
        cerr << "h56 dumped_pre_swap_pool path=" << dump_path
             << " size=" << keys.size() << '\n';
    }
    if ((int)pool.size() < core_pool && still_have_pool_time())
        h56_add_local_edit_components(pool, T1, T2, core_pool);
    if (h56_env_int("H56_ENABLE_SWAP_REFINEMENT",
                    H56_ENABLE_SWAP_REFINEMENT_VALUE) &&
        (int)pool.size() < core_pool && still_have_pool_time())
        h56_add_swap_refinement_components(pool, T1, T2, core_pool);
    if (h56_env_int("H56_ENABLE_GLOBAL_GROWTH", 0) &&
        (int)pool.size() < core_pool && still_have_pool_time()) {
        h56_add_global_growth_components(pool, T1, T2, core_pool);
    }
    if (h56_env_int("H56_ENABLE_SMALL_REGION_SUBSETS",
                    H56_ENABLE_SMALL_REGION_SUBSETS_VALUE) &&
        (int)pool.size() < core_pool && still_have_pool_time()) {
        h56_add_small_region_subset_components(pool, regions, T1, T2,
                                               core_pool, global_start,
                                               pool_cutoff_seconds);
    }
    const bool enable_external =
        h56_env_int("H56_ENABLE_EXTERNAL_POOL",
                    medium_instance
                        ? H56_MEDIUM_ENABLE_EXTERNAL_POOL_VALUE
                        : H56_ENABLE_EXTERNAL_POOL_VALUE) != 0;
    if (enable_external &&
        (int)pool.size() < max_pool && still_have_pool_time()) {
        int external_cap = h56_env_int("H56_EXTERNAL_POOL_LIMIT",
                                       medium_instance
                                           ? max_pool
                                           : H56_EXTERNAL_POOL_LIMIT_VALUE);
        external_cap = max((int)pool.size(), min(max_pool, external_cap));
        h56_add_external_like_pool(pool, T1, T2, external_cap,
                                   global_start, pool_cutoff_seconds);
    }
    cerr << "h56 targeted_pool before=" << before
         << " after=" << pool.size()
         << " regions=" << regions.size()
         << " elapsed=" << h56_elapsed_seconds(global_start)
         << " pool_cutoff=" << pool_cutoff_seconds << '\n';
}

#ifndef H56_NO_MAIN
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    stringstream input_buffer;
    input_buffer << cin.rdbuf();
    string input_text = input_buffer.str();
    if (input_text.empty()) return 1;

    TwoTreeInstance instance;
    stringstream instance_input(input_text);
    if (!read_two_tree_instance(instance_input, instance)) return 1;
    if (instance.leaf_count > 500) {
        while (true) pause();
    }
    Tree T1 = parse_newick(instance.trees[0], instance.leaf_count);
    Tree T2 = parse_newick(instance.trees[1], instance.leaf_count);
    if (T1.root < 0 || T2.root < 0) return 1;
    build_lca(T1);
    build_lca(T2);

    signal(SIGTERM, handle_signal);
    signal(SIGINT, handle_signal);

    auto global_start = chrono::steady_clock::now();
    const double total_seconds =
        h56_env_double("H56_TOTAL_SECONDS", H56_TOTAL_SECONDS_VALUE);
    const double output_safety_seconds =
        h56_env_double("H56_OUTPUT_SAFETY_SECONDS",
                       H56_OUTPUT_SAFETY_SECONDS_VALUE);
    const double highs_setup_seconds =
        h56_env_double("H56_HIGHS_SETUP_SECONDS",
                       H56_HIGHS_SETUP_SECONDS_VALUE);
    const double configured_highs_seconds = h46_highs_seconds();
    const double pool_cutoff_seconds =
        max(0.0, total_seconds - output_safety_seconds -
                     highs_setup_seconds - configured_highs_seconds);

    string prefix = "/tmp/pace_h56_" + to_string((long long)getpid());
    H56RouteFiles route;
    if (!h56_run_h45_route(
            input_text, prefix, h56_env_int("H56_SCORE_MODE", 0),
            h56_env_int("H56_PRIMARY_SEED", 0),
            h56_env_int("H56_BEST_SEED", 193),
            h56_env_double("H56_ROUTE_SECONDS", H56_ROUTE_SECONDS_VALUE),
            route)) {
        return 1;
    }

    AgreementForest best;
    bool have_best = false;
    AgreementForest previous_context;
    bool have_previous_context = false;
    vector<AgreementForest> milestone_contexts;
    unordered_map<string, vector<int>> component_pool;

    vector<string> candidate_paths;
    candidate_paths.push_back(route.output_path);
    candidate_paths.push_back(route.snapshot_path);
    for (const string &path : candidate_paths) {
        AgreementForest candidate;
        string output = h46_read_file(path);
        if (h46_parse_forest(output, instance.leaf_count, T1, T2, candidate)) {
            if (!have_best || active_count(candidate) < active_count(best)) {
                best = std::move(candidate);
                have_best = true;
            }
        }
    }
    {
        AgreementForest candidate;
        string output = h46_read_file(route.previous_snapshot_path);
        if (h46_parse_forest(output, instance.leaf_count, T1, T2,
                             candidate)) {
            previous_context = std::move(candidate);
            have_previous_context = true;
            cerr << "h56 previous_context components="
                 << active_count(previous_context) << '\n';
        }
    }
    for (const string &path : route.milestone_paths) {
        AgreementForest candidate;
        string output = h46_read_file(path);
        if (h46_parse_forest(output, instance.leaf_count, T1, T2,
                             candidate)) {
            cerr << "h56 milestone_context components="
                 << active_count(candidate)
                 << " path=" << path << '\n';
            milestone_contexts.push_back(std::move(candidate));
        }
    }
    h46_parse_pool(route.pool_path, instance.leaf_count, component_pool);
    if (!have_best) return 1;
    update_signal_best_output(forest_to_output(best, T1));

    for (const Component &component : best.comp) {
        if (!component.active || component.leaves.empty()) continue;
        component_pool.emplace(h45_component_key(component.leaves),
                               component.leaves);
    }
    h56_add_fast_greedy_pool(component_pool, best, T1, T2);
    for (const AgreementForest &milestone : milestone_contexts)
        h56_add_forest_components(component_pool, milestone);
    cerr << "h56 pre_target incumbent=" << active_count(best)
         << " raw_pool=" << component_pool.size() << '\n';

    if (h56_env_int("H56_ENABLE_MILESTONE_CONTEXTS", 1) &&
        !milestone_contexts.empty()) {
        sort(milestone_contexts.begin(), milestone_contexts.end(),
             [](const AgreementForest &a, const AgreementForest &b) {
                 int ka = active_count(a);
                 int kb = active_count(b);
                 if (ka != kb) return ka < kb;
                 return a.comp.size() < b.comp.size();
             });
        const char *old_core = getenv("H56_CORE_MAX_POOL");
        string old_core_value = old_core ? string(old_core) : string();
        bool had_old_core = old_core != nullptr;
        int milestone_cap = h56_env_int("H56_MILESTONE_POOL_CAP", 5200);
        milestone_cap = max((int)component_pool.size(), milestone_cap);
        string milestone_cap_text = to_string(milestone_cap);
        setenv("H56_CORE_MAX_POOL", milestone_cap_text.c_str(), 1);
        double milestone_seconds =
            h56_env_double("H56_MILESTONE_SECONDS", 24.0);
        double milestone_cutoff = min(
            pool_cutoff_seconds,
            h56_elapsed_seconds(global_start) + milestone_seconds);
        cerr << "h56 milestone_enrich start cap=" << milestone_cap
             << " contexts=" << milestone_contexts.size()
             << " cutoff=" << milestone_cutoff << '\n';
        for (const AgreementForest &milestone : milestone_contexts) {
            if (h56_elapsed_seconds(global_start) >= milestone_cutoff) break;
            if (active_count(milestone) >
                h56_env_int("H56_MILESTONE_MAX_COMPONENTS", 430)) {
                continue;
            }
            h56_enrich_targeted_medium(component_pool, milestone, T1, T2,
                                       global_start, milestone_cutoff);
            if ((int)component_pool.size() >= milestone_cap) break;
        }
        if (had_old_core)
            setenv("H56_CORE_MAX_POOL", old_core_value.c_str(), 1);
        else
            unsetenv("H56_CORE_MAX_POOL");
        cerr << "h56 milestone_enrich end pool="
             << component_pool.size()
             << " elapsed=" << h56_elapsed_seconds(global_start) << '\n';
    }

    if (h56_env_int("H56_ENABLE_EARLY_FLOOR_ILP", 1) &&
        active_count(best) > 425 &&
        h56_elapsed_seconds(global_start) <
            pool_cutoff_seconds - h56_env_double("H56_EARLY_HIGHS_SECONDS",
                                                 12.0) - 1.0) {
        const char *old_highs = getenv("H46_HIGHS_SECONDS");
        string old_highs_value = old_highs ? string(old_highs) : string();
        bool had_old_highs = old_highs != nullptr;
        string early_highs = to_string(
            h56_env_double("H56_EARLY_HIGHS_SECONDS", 12.0));
        setenv("H46_HIGHS_SECONDS", early_highs.c_str(), 1);
        AgreementForest floor_forest;
        int before_floor = active_count(best);
        cerr << "h56 early_floor_highs start incumbent=" << before_floor
             << " seconds=" << early_highs
             << " pool=" << component_pool.size()
             << " elapsed=" << h56_elapsed_seconds(global_start) << '\n';
        if (h46_solve_pool(component_pool, instance.leaf_count, T1, T2,
                           before_floor, floor_forest)) {
            best = std::move(floor_forest);
            update_signal_best_output(forest_to_output(best, T1));
            h56_add_forest_components(component_pool, best);
            cerr << "h56 early_floor_highs improved "
                 << before_floor << " -> " << active_count(best)
                 << " pool=" << component_pool.size()
                 << " elapsed=" << h56_elapsed_seconds(global_start)
                 << '\n';
        } else {
            cerr << "h56 early_floor_highs no_improvement incumbent="
                 << before_floor
                 << " elapsed=" << h56_elapsed_seconds(global_start)
                 << '\n';
        }
        if (had_old_highs)
            setenv("H46_HIGHS_SECONDS", old_highs_value.c_str(), 1);
        else
            unsetenv("H46_HIGHS_SECONDS");
    }

    if (have_previous_context &&
        h56_env_int("H56_ENABLE_PREVIOUS_CONTEXT", 0) &&
        active_count(best) <= 425 &&
        active_count(previous_context) <= active_count(best) + 3) {
        const char *old_core = getenv("H56_CORE_MAX_POOL");
        string old_core_value = old_core ? string(old_core) : string();
        bool had_old_core = old_core != nullptr;
        int context_cap =
            h56_env_int("H56_PREVIOUS_CONTEXT_POOL_CAP", 9000);
        context_cap = max((int)component_pool.size(), context_cap);
        string context_cap_text = to_string(context_cap);
        setenv("H56_CORE_MAX_POOL", context_cap_text.c_str(), 1);
        cerr << "h56 previous_context_enrich start cap="
             << context_cap
             << " context=" << active_count(previous_context)
             << " incumbent=" << active_count(best) << '\n';
        h56_enrich_targeted_medium(component_pool, previous_context,
                                   T1, T2, global_start,
                                   pool_cutoff_seconds);
        if (had_old_core)
            setenv("H56_CORE_MAX_POOL", old_core_value.c_str(), 1);
        else
            unsetenv("H56_CORE_MAX_POOL");
        cerr << "h56 previous_context_enrich end pool="
             << component_pool.size()
             << " elapsed=" << h56_elapsed_seconds(global_start) << '\n';
    }

    h56_enrich_targeted_medium(component_pool, best, T1, T2, global_start,
                               pool_cutoff_seconds);
    if (const char *dump_path = getenv("H56_DUMP_POOL")) {
        ofstream dump(dump_path);
        vector<string> keys;
        keys.reserve(component_pool.size());
        for (const auto &entry : component_pool)
            keys.push_back(h45_component_key(entry.second));
        sort(keys.begin(), keys.end());
        keys.erase(unique(keys.begin(), keys.end()), keys.end());
        for (const string &key : keys) dump << key << '\n';
        cerr << "h56 dumped_pool path=" << dump_path
             << " size=" << keys.size() << '\n';
    }

    const int ilp_passes = max(
        1, h56_env_int("H56_ILP_PASSES", H56_ILP_PASSES_VALUE));
    const double min_highs_seconds =
        h56_env_double("H56_MIN_HIGHS_SECONDS",
                       H56_MIN_HIGHS_SECONDS_VALUE);
    const bool medium_instance =
        T1.n >= h56_env_int("H56_MEDIUM_INSTANCE_MIN_LEAVES",
                            H56_MEDIUM_INSTANCE_MIN_LEAVES_VALUE);
    for (int pass = 0; pass < ilp_passes && !g_stop_requested; pass++) {
        AgreementForest optimized;
        double elapsed = h56_elapsed_seconds(global_start);
        double available_for_highs =
            total_seconds - elapsed - output_safety_seconds -
            highs_setup_seconds;
        double highs_seconds =
            min(configured_highs_seconds, available_for_highs);
        if (highs_seconds >= min_highs_seconds) {
            string highs_value = to_string(highs_seconds);
            setenv("H46_HIGHS_SECONDS", highs_value.c_str(), 1);
            int before_k = active_count(best);
            cerr << "h56 highs_budget pass=" << pass
                 << " seconds=" << highs_seconds
                 << " elapsed=" << elapsed
                 << " total=" << total_seconds
                 << " incumbent=" << before_k << '\n';
            if (h46_solve_pool(component_pool, instance.leaf_count, T1, T2,
                               active_count(best), optimized)) {
                best = std::move(optimized);
                update_signal_best_output(forest_to_output(best, T1));
                int after_k = active_count(best);
                h56_add_forest_components(component_pool, best);
                cerr << "h56 ilp_pass_improved pass=" << pass
                     << " " << before_k << " -> " << after_k
                     << " pool=" << component_pool.size()
                     << " elapsed=" << h56_elapsed_seconds(global_start)
                     << '\n';
                if (pass + 1 < ilp_passes) {
                    int second_max_pool = h56_env_int(
                        "H56_SECOND_MAX_POOL", H56_SECOND_MAX_POOL_VALUE);
                    if (medium_instance && second_max_pool <= 0)
                        second_max_pool = H56_MEDIUM_MAX_POOL_VALUE;
                    if (second_max_pool > 0) {
                        string value = to_string(second_max_pool);
                        setenv("H56_MAX_POOL", value.c_str(), 1);
                        cerr << "h56 second_max_pool=" << second_max_pool
                             << '\n';
                    }
                    int second_external_cap = h56_env_int(
                        "H56_SECOND_EXTERNAL_POOL_LIMIT",
                        H56_SECOND_EXTERNAL_POOL_LIMIT_VALUE);
                    if (medium_instance && second_external_cap <= 0)
                        second_external_cap = H56_MEDIUM_MAX_POOL_VALUE;
                    if (second_external_cap > 0) {
                        string value = to_string(second_external_cap);
                        setenv("H56_EXTERNAL_POOL_LIMIT", value.c_str(), 1);
                        cerr << "h56 second_external_pool_limit="
                             << second_external_cap << '\n';
                    }
                    if (medium_instance &&
                        h56_env_int("H56_SECOND_ENABLE_EXTERNAL_POOL", 1)) {
                        setenv("H56_ENABLE_EXTERNAL_POOL", "1", 1);
                        cerr << "h56 second_external_pool_enabled=1\n";
                    }
                    double reenrich_seconds = h56_env_double(
                        "H56_REENRICH_SECONDS", H56_REENRICH_SECONDS_VALUE);
                    double reenrich_cutoff = min(
                        total_seconds - output_safety_seconds -
                            highs_setup_seconds - min_highs_seconds,
                        h56_elapsed_seconds(global_start) + reenrich_seconds);
                    if (medium_instance &&
                        h56_env_int("H56_SECOND_ENABLE_COMPLEMENTS", 0) &&
                        second_max_pool > (int)component_pool.size() &&
                        reenrich_cutoff > h56_elapsed_seconds(global_start)) {
                        vector<H56Region> improved_regions =
                            h56_collect_regions(best, T1, T2);
                        h56_add_region_complement_components(
                            component_pool, improved_regions, T1, T2,
                            second_max_pool, global_start, reenrich_cutoff);
                    }
                    if (medium_instance &&
                        h56_env_int("H56_SECOND_ENABLE_COMPLETION", 0) &&
                        second_max_pool > (int)component_pool.size() &&
                        reenrich_cutoff > h56_elapsed_seconds(global_start)) {
                        setenv("H56_ENABLE_COMPLETION", "1", 0);
                        setenv("H56_COMPLETION_SECONDS", "4.0", 0);
                        setenv("H56_COMPLETION_MAX_REMOVE", "6", 0);
                        setenv("H56_COMPLETION_MAX_ADD", "14", 0);
                        setenv("H56_COMPLETION_BASE_LIMIT", "1400", 0);
                        setenv("H56_COMPLETION_ADD_LIMIT", "1600", 0);
                        setenv("H56_COMPLETION_OPTION_LIMIT", "80", 0);
                        setenv("H56_COMPLETION_BEAM_WIDTH", "16", 0);
                        setenv("H56_COMPLETION_MAX_COMPONENT_SIZE", "40", 0);
                        int before_completion = (int)component_pool.size();
                        h56_add_completion_growth_components(
                            component_pool, T1, T2, second_max_pool,
                            global_start, reenrich_cutoff);
                        cerr << "h56 second_completion before="
                             << before_completion
                             << " after=" << component_pool.size()
                             << " elapsed="
                             << h56_elapsed_seconds(global_start) << '\n';
                    }
                    if (medium_instance &&
                        h56_env_int("H56_SECOND_ENABLE_EXTERNAL_POOL", 1)) {
                        string core_lock =
                            to_string((int)component_pool.size());
                        setenv("H56_CORE_MAX_POOL", core_lock.c_str(), 1);
                        cerr << "h56 second_core_pool_locked="
                             << core_lock << '\n';
                    }
                    if (reenrich_cutoff > h56_elapsed_seconds(global_start)) {
                        h56_enrich_targeted_medium(
                            component_pool, best, T1, T2, global_start,
                            reenrich_cutoff);
                    }
                }
            } else {
                cerr << "h56 ilp_pass_no_improvement pass=" << pass
                     << " incumbent=" << before_k
                     << " elapsed=" << h56_elapsed_seconds(global_start)
                     << '\n';
                break;
            }
        } else {
            cerr << "h56 highs_skip pass=" << pass
                 << " elapsed=" << elapsed
                 << " available=" << available_for_highs
                 << " min=" << min_highs_seconds << '\n';
            break;
        }
    }

    if (const char *dump_path = getenv("H56_DUMP_FINAL_POOL")) {
        ofstream dump(dump_path);
        vector<string> keys;
        keys.reserve(component_pool.size());
        for (const auto &entry : component_pool)
            keys.push_back(h45_component_key(entry.second));
        sort(keys.begin(), keys.end());
        keys.erase(unique(keys.begin(), keys.end()), keys.end());
        for (const string &key : keys) dump << key << '\n';
        cerr << "h56 dumped_final_pool path=" << dump_path
             << " size=" << keys.size() << '\n';
    }

    AgreementForest output_forest = exact_output_repair(best, T1, T2);
    update_signal_best_output(forest_to_output(output_forest, T1));
    cerr << "h56 final components=" << active_count(output_forest) << '\n';
    cout << forest_to_output(output_forest, T1);
    cout.flush();

    unlink(route.output_path.c_str());
    unlink(route.pool_path.c_str());
    unlink(route.snapshot_path.c_str());
    unlink(route.previous_snapshot_path.c_str());
    for (const string &path : route.milestone_paths) unlink(path.c_str());
    return 0;
}
#endif

// ===== END heuristic/clean_refactor/h56_targeted_medium.cpp =====
