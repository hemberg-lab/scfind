#include <Rcpp.h>
#include <iostream>

class FPNode;
class Pattern;

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
std::vector<bool> int2bin(int n) {
  std::vector<bool> b;
  while (n>=0) {
    if (n==0) {
      b.insert(b.begin(), false);
      break;
    }
    else if (n==1) {
      b.insert(b.begin(), true);
      break;
    }
    else {
      b.insert(b.begin(), (n % 2)==1);
      n /= 2;
    }
  }
  return b;
}

// [[Rcpp::export]]
LogicalVector quantizeCpp(NumericVector ids, int n) {
  std::vector<bool> tf;
  //Quantize each of the values to n bits
  for (int i=0; i<ids.size(); i++) {
    std::vector<bool> b = int2bin(ids[i]*pow(2.0, n));
    if (b.size()<n) { //pre-pend zeros to make length n
      std::vector<bool> tmp (n-b.size(), false);
      b.insert(b.begin(), tmp.begin(), tmp.end());
    }
    tf.insert(tf.end(), b.begin()+b.size()-n, b.end());
  }
  LogicalVector v(tf.size());
  for (int i=0; i<tf.size(); i++)
    v[i] = tf[i];
  return v;
}

// [[Rcpp::export]]
NumericVector dequantizeCpp(LogicalVector bits, int n) {
  std::vector<int> e(bits.size()/n);
  //Quantize each of the values to n bits
  for (int i=0; i<bits.size()/n; i++) {
    for (int j=0; j<n; j++) {
      if (bits[i*n+j])
	e[i] += pow(2.0, n-j-1);
    }
  }
  NumericVector v(e.size());
  for (int i=0; i<e.size(); i++)
    v[i] = e[i];
  return v;
}

// [[Rcpp::export]]
Rcpp::List eliasFanoCodingCpp(NumericVector ids, int l) {
  std::vector<bool> L, H;
  //Append 0 to the start of the list
  for (int i=0; i<ids.size(); i++) {
    std::vector<bool> c = int2bin(ids[i]);
    if (c.size()<l) { //pre-pend zeros to make length l
      std::vector<bool> tmp (l-c.size(), false);
      c.insert(c.begin(), tmp.begin(), tmp.end());
    }
    L.insert(L.end(), c.begin()+c.size()-l, c.end()); //append the lower bits
    //Use a unary code for the high bits
    int m = (int)(ids[i]/pow(2, l)) - (int)(ids[i-1]/pow(2, l));
    std::vector<bool> h (m+1, false);
    h[m] = true;
    H.insert(H.end(), h.begin(), h.end());
  }
  return List::create(_["H"] = H, _["L"] = L);
}


// [[Rcpp::export]]
NumericVector eliasFanoDecodingCpp(LogicalVector H, LogicalVector L, int m, int l, int ih, int il) {
  NumericVector ids(m);
  int nZeros = 0, nCellsFound = 0, i = ih, j = 0, prevH = 0;
  while (nCellsFound<m) {
    if (!H[i])
      nZeros++;
    else {
      for (int k=0; k<l; k++) {//Calculate the value stored in the low bits
	if (L[il+j*l+k])
	  ids[nCellsFound] += pow(2, l-k-1);
      }
      //Calculate the value stored in the high bits
      std::vector<bool> h = int2bin(nZeros + prevH);
      for (int k=0; k<h.size(); k++) {
	if (h[k])
	  ids[nCellsFound] += pow(2, h.size()-k+l-1);
      }
      prevH = (int)(ids[nCellsFound]/pow(2, l));
      j++;
      nZeros = 0;
      nCellsFound++;
    }
    i++;
  }
  return ids;
}

// [[Rcpp::export]]
NumericVector eliasFanoDecodingGeneCpp(LogicalVector Hs, LogicalVector Ls, NumericVector ms,  NumericVector ihs, NumericVector ns, NumericVector permutation) {
  //Need to figure out how to decode H and L since the parameter l varies. We know that there are ls[i]*ms[i] bits used in L for the ith entry. Similarly, we know that we need to read until we have encountered ms[i] 1s in H.
  std::vector<int> idsOriginal;
  int il = 0, idsOffset = 0; 
  for (int i=0; i<ns.size(); i++) {
    if (ms[i]>0) {
      int l = floor(log2(ns[i]/ms[i]));
      NumericVector tmp = eliasFanoDecodingCpp(Hs, Ls, ms[i], l, ihs[i]-1, il);
      //Now we need to look in the permutation vector to figure out the original ID for the cell
      for (int j=0; j<tmp.size(); j++) 
	idsOriginal.push_back(permutation[idsOffset + tmp[j] - 1]);
      //Figure out the index where we should start in the L vector
      il += l*ms[i];
    }
    idsOffset += ns[i];
  }
  NumericVector v(idsOriginal.size());
  for (int i=0; i<idsOriginal.size(); i++)
    v[i] = idsOriginal[i];
  return v;
}

// [[Rcpp::export]]
Rcpp::List eliasFanoDecodingGeneByClusterCpp(LogicalVector Hs, LogicalVector Ls, NumericVector ms,  NumericVector ihs, NumericVector ns, NumericVector permutation) {
  //Similar to the function above, but returns a list where each entry corresponds to a vector with indexes for that experiment
  Rcpp::List ret = List::create();
  std::vector<int> idsOriginal;
  int il = 0, idsOffset = 0, experimentInd = 0; 
  for (int i=0; i<ns.size(); i++) {
    if (ms[i]>0) {
      int l = floor(log2(ns[i]/ms[i]));
      NumericVector tmp = eliasFanoDecodingCpp(Hs, Ls, ms[i], l, ihs[i]-1, il);
      //Now we need to look in the permutation vector to figure out the original ID for the cell
      for (int j=0; j<tmp.size(); j++) 
	idsOriginal.push_back(permutation[idsOffset + tmp[j] - 1]);
      //Figure out the index where we should start in the L vector
      il += l*ms[i];
      NumericVector v(idsOriginal.size());
      for (int j=0; j<idsOriginal.size(); j++)
	v[j] = idsOriginal[j];
      ret.push_back(v);
      idsOriginal.clear();
    }
    else
      ret.push_back(NumericVector());
    idsOffset += ns[i];
  }
  return ret;
}

// [[Rcpp::export]]
Rcpp::List eliasFanoDecodingGeneByExperimentCpp(LogicalVector Hs, LogicalVector Ls, NumericVector ms,  NumericVector ihs, NumericVector ns, NumericVector permutation, NumericVector experimentStarts) {
  //Similar to the function above, but returns a list where each entry corresponds to a vector with indexes for that experiment
  Rcpp::List ret = List::create();
  std::vector<int> idsOriginal;
  int il = 0, idsOffset = 0, experimentInd = 0; 
  for (int i=0; i<ns.size(); i++) {
    if (ms[i]>0) {
      int l = floor(log2(ns[i]/ms[i]));
      NumericVector tmp = eliasFanoDecodingCpp(Hs, Ls, ms[i], l, ihs[i]-1, il);
      //Now we need to look in the permutation vector to figure out the original ID for the cell
      for (int j=0; j<tmp.size(); j++) 
	idsOriginal.push_back(permutation[idsOffset + tmp[j] - 1]);
      //Figure out the index where we should start in the L vector
      il += l*ms[i];
    }
    idsOffset += ns[i];
    if (idsOffset>=experimentStarts[experimentInd]) {
      NumericVector v(idsOriginal.size());
      for (int j=0; j<idsOriginal.size(); j++)
	v[j] = idsOriginal[j];
      ret.push_back(v);
      idsOriginal.clear();
      experimentInd++;
    }
  }
  return ret;
}



//////////////////////////////////////////////////////////////////

class Pattern {
public:
  Pattern();
  Pattern(int id, bool itemsGenes) { 
    if(itemsGenes) 
      cellIDs.insert(id); 
    else
      geneIDs.insert(id);
  }

  void addGene(int geneID) {
    geneIDs.insert(geneID);
  }

  void addCell(int cellID) {
    cellIDs.insert(cellID);
  }

  int utility() {
    return (geneIDs.size()-1)*(cellIDs.size()-1);
  }

  void print() {
    Rprintf("Cells: ");
    for (std::set<int>::iterator it=cellIDs.begin(); it!=cellIDs.end(); ++it)
      Rprintf("%d, ", *it);
    Rprintf("Genes: ");
    for (std::set<int>::iterator it=geneIDs.begin(); it!=geneIDs.end(); ++it)
      Rprintf("%d, ", *it);
    Rprintf("\n");
  }

  Rcpp::List asList() {
    return List::create(_["cells"] = cellIDs, _["genes"] = geneIDs);
  }

  std::set<int> cellIDs, geneIDs;
};

/*RCPP_EXPOSED_CLASS(Pattern)

RCPP_MODULE(PatternMod) {
  
  class_<FPNode>("Pattern")

    .constructor<int>("Creates the root")

    .method("getChild", &FPNode::getChild, "Returns the child if it exists")
    .method("addChild", &FPNode::addChild, "Adds a child")
    .method("inc", &FPNode::inc, "Increases the counter")
    .method("addItem", &FPNode::addItem, "Adds as an item")
    .method("getID", &FPNode::getID, "Returns ID")
    .method("print", &FPNode::print, "prints")
    .method("generatePotentialItemsetsList", &FPNode::generatePotentialItemsetsList, "Finds combinations of genes")
    ;
    }*/


class FPNode {
public: 
  FPNode() : id(-1), depth(0) { count = 0; marked = false; parent = NULL; }
  FPNode(int id_) : id(id_), depth(0) { count = 1; marked = false; parent = NULL; }
  FPNode(int id_, FPNode* parent_) : id(id_), parent(parent_) { 
    count = 1; 
    marked = false; 
    depth = parent->depth + 1;
  }
  
  bool hasChild(int j) {
    for (int i=0; i<children.size(); i++) {
      if (j==children[i]->id)
	return true;
    }
    return false;
  }

  FPNode* getChild(int j) {
    for (int i=0; i<children.size(); i++) {
      if (j==children[i]->id)
	return children[i];
    }
    return (new FPNode());
  }
  
  FPNode* addChild(int j) {
    FPNode* n = new FPNode(j, this);
    children.push_back(n);
    return n;
  }
  
  void inc() {
    count++;
  }

  int getID() { return id; }

  void addItem(int i) {
    items.push_back(i);
  }

  int numberOfNodes() {
    int n = 1;
    for (int i=0; i<children.size(); i++) 
      n += children[i]->numberOfNodes();
    return n;
  }

  int maxDepth() {
    int d = depth;
    for (int i=0; i<children.size(); i++) {
      int d2 = children[i]->maxDepth();
      if (d2>d)
	d = d2;
    }
    return d;
  }

  int maxNumberOfItems() {
    int n = items.size();
    for (int i=0; i<children.size(); i++) {
      int n2 = children[i]->maxNumberOfItems();
      if (n2>n)
	n = n2;
    }
    return n;
  }

  void print() {
    printChild(id);
  }

  void printChild(int parentID) {
    Rprintf("%d\t%d\t%d\t%d\t%d\t", id, parentID, count, children.size(), items.size());
//  if (items.size()>0) {
//    for (int i=0; i<items.size(); i++)
//      Rprintf("%d,", items[i]);
//  }
    Rprintf("\n");
    for (int i=0; i<children.size(); i++) 
      children[i]->printChild(id);
  }

  void prune(int minCount) { //Recursively remove all nodes that fall below minCount
    for (std::vector<FPNode*>::iterator it=children.begin(); it!=children.end(); it++) {
      FPNode* f = *it;
      if (f->count<minCount) 
	children.erase(it--);
      else 
	f->prune(minCount);
    }
  }

  void pruneItems(int minItems) { //Recursively remove all nodes where the number of items fall below minItems
    for (std::vector<FPNode*>::iterator it=children.begin(); it!=children.end(); it++) {
      FPNode* f = *it;
      if (f->items.size()<minItems) 
	children.erase(it--);
      else 
	f->pruneItems(minItems);
    }
  }

  std::vector<Pattern*> generatePotentialItemsets(int minGenes, int minCells, std::vector<Pattern*> l, bool itemsGenes) {
    //Traverse the tree to the leaves and as we go back up to the root, we create potential itemsets
    for (int i=0; i<children.size(); i++) { 
      FPNode* c = children[i];
      if (c->items.size()>0) 
	l = c->generatePotentialItemsets(minGenes, minCells, l, itemsGenes);
    }
    if (id!=-1 && children.size()==0 && ((itemsGenes && depth>=minCells && items.size()>=minGenes) || (~itemsGenes && depth>=minGenes && items.size()>=minCells))) {
      Pattern* pattern = new Pattern(id, itemsGenes);
      Rprintf("Creating pattern from node %d with %d items at depth %d\n", id, items.size(), depth);
      for (int i=0; i<items.size(); i++) {
	if (itemsGenes) 
	  pattern->addGene(items[i]);
	else
	  pattern->addCell(items[i]);
	//Rprintf("%d,", items[i]);
      }
      //Rprintf("\n");
      pattern = this->markNode(pattern, minGenes, minCells, itemsGenes);
      Rprintf("Found pattern with utility %d, %d genes and %d cells.\n", pattern->utility(), pattern->geneIDs.size(), pattern->cellIDs.size());
      if (pattern->utility()>(minGenes-1)*(minCells-1))
	l.push_back(pattern);
    }
    return l;
  }

  Rcpp::List generatePotentialItemsetsList(int minGenes, int minCells, bool itemsGenes) {
    std::vector<Pattern*> potentialItemsets;
    potentialItemsets = generatePotentialItemsets(minGenes, minCells, potentialItemsets, itemsGenes);
    //Convert to a list and return to R
    Rcpp::List itemsets = List::create();
    for (int i=0; i<potentialItemsets.size(); i++)
      itemsets.push_back(potentialItemsets[i]->asList());
    return itemsets;
  }

  Pattern* markNode(Pattern* pattern, int minGenes, int minCells, bool itemsGenes) {
    int n = items.size();
    if (!marked && n>0 && id!=-1) {
      marked = true;
      //Consider either adding the gene/cell at this level to the pattern, or discard the items (genes/cells) collected so far and use the ones found at this node. This works since we know that the sets of items get progressively smaller as we climb down the tree. Thus, we can calculate the maximum possible size of a pattern and compare
      //bool createNewPattern = false;
      //int newUtility = depth*items.size(), oldUtility = 0;
      if (itemsGenes) {
	//oldUtility = (depth+pattern->cellIDs.size())*pattern->geneIDs.size();
	//if (newUtility<oldUtility) 
	  pattern->addCell(id);
	  //else 
	  //createNewPattern = true;	
      }
      else {
	//oldUtility = (depth+pattern->geneIDs.size())*pattern->cellIDs.size();
	//if (newUtility<oldUtility)
	  pattern->addGene(id);
	  //else 
	  //createNewPattern = true;
      }
      //if (createNewPattern) {
	//Rprintf("Found better pattern with depth %d with utility %d instead of %d\n", depth, newUtility, oldUtility);
	/*pattern = new Pattern(id, itemsGenes);
	for (int i=0; i<items.size(); i++) {
	  if (itemsGenes) 
	    pattern->addGene(items[i]);
	  else 
	    pattern->addCell(items[i]);
	    }*/
      //}
      parent->markNode(pattern, minGenes, minCells, itemsGenes);
    }
    return pattern;
  }

  Rcpp::List asList() {
    //Recursively convert all children to lists first
    List l = List::create();
    for (int i=0; i<children.size(); i++) 
      l.push_back(children[i]->asList());
    IntegerVector v;
    for (int i=0; i<items.size(); i++) 
      v.push_back(items[i]);
    List list = List::create(_["count"] = count, _["id"] = id, _["children"] = l, _["items"] = v);
    return list;
  }
  
  //private:
  int id, count, depth;
  std::vector<FPNode*> children;
  FPNode* parent;
  std::vector<int> items;
  bool marked;
};

RCPP_EXPOSED_CLASS(FPNode)

RCPP_MODULE(FPNodeMod) {
  
  class_<FPNode>("FPNode")

    .constructor<int>("Creates the root")

    .method("getChild", &FPNode::getChild, "Returns the child if it exists")
    .method("addChild", &FPNode::addChild, "Adds a child")
    .method("inc", &FPNode::inc, "Increases the counter")
    .method("addItem", &FPNode::addItem, "Adds as an item")
    .method("getID", &FPNode::getID, "Returns ID")
    .method("print", &FPNode::print, "prints")
    .method("generatePotentialItemsetsList", &FPNode::generatePotentialItemsetsList, "Finds combinations of genes")
    ;
}

// [[Rcpp::export]]
NumericVector findUniqueIndexCpp(NumericVector cellInds, NumericVector permutation, int permutationStartInd, int permutationEndInd) {
  //This function is equivalent to "which(cellInds %in% permutation)" operation in R
  NumericVector inds;
  //Use double for-loop, not efficient
  for (int i=0; i<cellInds.length(); i++) {
    for (int j=permutationStartInd; j<permutationEndInd; j++) {
      if (cellInds[i]==permutation[j]) {
	inds.push_back(j);
	break;
      }
    }
  }
  return inds;
}

// [[Rcpp::export]]
Rcpp::List findAssociationRulesCellsFPgrowthCpp(int minGenes, int minCells, Rcpp::List &Hss, Rcpp::List &Lss, NumericMatrix mss,  NumericMatrix ihss, NumericVector ns, NumericVector permutation, NumericVector cellsPerExperimentStart, bool transpose) { //For some reason minGenes and minCells are swapped, otherwise this seems to be workign very well now...
  FPNode* fproot = new FPNode(-1);
  if (transpose) {
    for (int i=0; i<Hss.size(); i++) {
      FPNode* n = fproot;
      List cellInds = eliasFanoDecodingGeneByExperimentCpp(Hss[i], Lss[i], mss(i, _ ), ihss(i, _ ), ns, permutation, cellsPerExperimentStart);
      int start = 0;
      for (int j=0; j<cellInds.size(); j++) {
	//Find out the unique ID from this gene by comparing to its index in the permutation vector
	NumericVector ids = findUniqueIndexCpp(cellInds[j], permutation, start, cellsPerExperimentStart[j]);
	start = cellsPerExperimentStart[j];
	for (int k=0; k<ids.size(); k++) {
	  FPNode* child = n->getChild(ids[k]); //nodes correspond to cells
	  if (child->getID()==-1) {
	    n->addChild(ids[k]);
	    child = n->getChild(ids[k]);
	  }
	  child->inc();
	  child->addItem(i); //Items correspond to genes
	  n = child;
	}
      }
    }
  }
  else {
    //We need to first extract all of the information and put in a transposed matrix.
    int nCells = sum(ns), nGenes = Hss.size();
    bool mb[nGenes][nCells];
    LogicalMatrix m(nGenes, nCells);
    for (int i=0; i<nGenes; i++) {
      for (int k=0; k<nCells; k++) { m[i, k] = false; mb[i][k] = false; }
      List cellInds = eliasFanoDecodingGeneByExperimentCpp(Hss[i], Lss[i], mss(i, _ ), ihss(i, _ ), ns, permutation, cellsPerExperimentStart);
      int start = 0;
      for (int j=0; j<cellInds.size(); j++) {
	//Find out the unique ID from this gene by comparing to its index in the permutation vector
	NumericVector ids = findUniqueIndexCpp(cellInds[j], permutation, start, cellsPerExperimentStart[j]);
	start = cellsPerExperimentStart[j];
	//Rprintf("ids: %.2f, %.2f, %.2f, %.2f\n", ids[0], ids[1], ids[2], ids[3]);
	for (int k=0; k<ids.size(); k++) {
	  m[i, ids[k]] = true;	  
	  mb[i][(int)ids[k]] = true;
	}
      }
      //Rprintf("m%d: %i, %i, %i, %i\n", i, m[i,0], m[i,1], m[i,2], m[i,3]);
      //Rprintf("mb%d: %i, %i, %i, %i\n", i, mb[i][0], mb[i][1], mb[i][2], mb[i][3]);
    }
    for (int i=0; i<nGenes; i++) {
      FPNode* n = fproot;
      //Rprintf("m%d: %i, %i, %i, %i\n", i, m[i,0], m[i,1], m[i,2], m[i,3]);
      //Rprintf("mb%d: %i, %i, %i, %i\n", i, mb[i][0], mb[i][1], mb[i][2], mb[i][3]);
      for (int j=0; j<nCells; j++) {
	if (mb[i][j]) {
	  FPNode* child = n->getChild(j); //nodes correspond to genes
	  if (child->getID()==-1) {
	    n->addChild(j);
	    child = n->getChild(j);
	  }
	  child->inc();
	  child->addItem(i); //Items correspond to cells
	  n = child;
	}
      }
    }
  }
  Rprintf("Built tree with %d nodes and depth %d with max items %d before pruning\n", fproot->numberOfNodes(), fproot->maxDepth(), fproot->maxNumberOfItems());
  fproot->prune(transpose ? minCells : minGenes);
  Rprintf("Kept %d nodes and depth %d with max items %d after pruning\n", fproot->numberOfNodes(), fproot->maxDepth(), fproot->maxNumberOfItems());
  fproot->pruneItems(transpose ? minGenes : minCells);
  Rprintf("Kept %d nodes and depth %d with max items %d after pruning items\n", fproot->numberOfNodes(), fproot->maxDepth(), fproot->maxNumberOfItems());
  Rcpp::List potentialItemsets = fproot->generatePotentialItemsetsList(transpose ? minGenes : minCells, transpose ? minCells : minGenes, true);
  return List::create(/*_["fp"] = fproot,*/ _["potential.itemsets"] = potentialItemsets);
}
 
// [[Rcpp::export]]
Rcpp::List buildFPTree(IntegerMatrix m, int minCells, int minGenes) {
  int nCells = m.nrow(), nGenes = m.ncol();
  FPNode* fproot = new FPNode(-1);
  for (int i=0; i<nCells; i++) {
    FPNode* n = fproot;
    for (int j=0; j<nGenes; j++) {
      if (m(i, j)>0) {
	int childID = n->getChild(j)->id; //nodes correspond to genes
	//FPNode *c = n->getChild(j);
	if (childID==-1) 
	  n = n->addChild(j);
	else {
	  n->getChild(j)->inc();
	  n->getChild(j)->addItem(i); //items correspond to cells
	  n = n->getChild(j);
	}
      }
    }
  }
  //Now we can prune the tree to remove all nodes where the count is too low
  fproot->prune(minCells);
  fproot->print();
  //return fproot->asList();
  std::vector<Pattern*> potentialItemsets;
  potentialItemsets = fproot->generatePotentialItemsets((minGenes-1)*(minCells-1), minGenes, potentialItemsets, false);
  Rprintf("Found %d itemsets\n", potentialItemsets.size());
  //Convert to a list and return to R
  Rcpp::List itemsets = List::create();//_["cells"] = itemsets[0]->cellIDs, _["genes"] = itemsets[0]->geneIDs);
  for (int i=0; i<potentialItemsets.size(); i++)
    itemsets.push_back(potentialItemsets[i]->asList());
  return itemsets;
}

// [[Rcpp::export]]
Rcpp::List buildFPTreeT(IntegerMatrix m, int minGenes, int minCells) {
  //Same as before, but using the transposed matrix, reading it one row at a time where each row corresponds t
  int nCells = m.ncol(), nGenes = m.nrow();
  FPNode* fproot = new FPNode(-1);
  for (int i=0; i<nGenes; i++) {
    FPNode* n = fproot;
    for (int j=0; j<nCells; j++) {
      if (m(i, j)>0) {
	int childID = n->getChild(j)->id; //nodes correspond to cells
	if (childID==-1) 
	  n->addChild(j);
	n->getChild(j)->inc();
	n->getChild(j)->addItem(i); //items correspond to genes
	n = n->getChild(j);
      }
    }
  }
  //Now we can prune the tree to remove all nodes where the count is too low
  Rprintf("ID\tParentID\tcount\t#children\titems\n");
  //fproot->print();
  fproot->pruneItems(minGenes);
  fproot->print();
  //Now we can identify the potential itemsets of interest through a depth-first search
  std::vector<Pattern*> potentialItemsets;
  potentialItemsets = fproot->generatePotentialItemsets((minGenes-1)*(minCells-1), minCells, potentialItemsets, true);
  Rprintf("Found %d itemsets\n", potentialItemsets.size());
  //Convert to a list and return to R
  Rcpp::List itemsets = List::create();//_["cells"] = itemsets[0]->cellIDs, _["genes"] = itemsets[0]->geneIDs);
  for (int i=0; i<potentialItemsets.size(); i++)
    itemsets.push_back(potentialItemsets[i]->asList());
  return itemsets;
}




//NumericMatrix findTriplets(int minIntersection, List Hs, List Ls, NumericMatrix ms, NumericMatrix ls) {
//  std::forward_list<int[]> triplets;
//  for (int i=0; i<ms.nrow(); i++) {
//    NumericVector ids1 = eliasFanoDecodingCpp(Hs[i], Ls[i], ms(i,_), ls(i,_), 0, 0);
//    for (int j=i+1; j<ms.nrow(); j++) {
//      NumericVector ids2 = eliasFanoDecodingCpp(Hs[j], Ls[j], ms(j,_), ls(j,_), 0, 0);
//      //Calculate intersection btwn i and j
//      NumericVector ids12 = sortedIntersect(ids1, ids2);
//      for (int k=j+1; k<ms.nrow(); k++) {
//	NumericVector ids3 = eliasFanoDecodingCpp(Hs[k], Ls[k], ms(k,_), ls(k,_), 0, 0);
//	//Do the intersection
//	NumericVector ids = sortedIntersect(ids12, ids3);
//	if (ids.size()>minIntersection) {
//	  triplets.push_front(new int[]{i, j, k});
//	}
//      }
//    }
//  }
//  //Store in a matrix instead
//  return triplets;
//}

//Rcpp::List createIndex(NumericMatrix d, int u) {
  //Create Elias-Fano index for the genes in the matrix d. It is assumed that they all come from the same cell-type
  //NumericVector m(d.nrow());
  //for (int i=0; i<d.nrow(); i++) {
    //Find the indexes of the cells expressing the gene
    
    //if (m[i]>0) {
//int l = floor(log2(u/m[i]));
