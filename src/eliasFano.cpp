#include <Rcpp.h>
#include <iostream>
#include <bitset>
#include <algorithm>
#include <cmath>

#include "typedefs.h"
// the bits used for the encoding
#define BITS 32


using namespace Rcpp;

typedef std::pair<unsigned short, std::bitset<BITS> > BitSet32;

std::bitset<BITS> int2bin_core(const unsigned int id)
{
  // we have a fixed BITS bitset so we have to offset 
  // the leading zeros whenever we need to access the elements
  return std::bitset<BITS>(id);
}



BitSet32 int2bin_bounded(const unsigned int id, unsigned int min_bit_length)
{
  // Check if number is larger than the desired size in bits
  unsigned short bit_set_length = id == 0 ? 0 : __builtin_clz(id);
  bit_set_length = bit_set_length > min_bit_length ? bit_set_length : min_bit_length;
  return make_pair(bit_set_length, int2bin_core(id));
}


// These functions work with 1-based indexes on arrays, 
// __builtin_clz has undefined behavior with 0
BitSet32 int2bin(unsigned int id)
{
  return make_pair(__builtin_clz(id), int2bin_core(id));
}

// [[Rcpp::export]]
BoolVec quantizeCpp(IntegerVector ids, int n)
{
  BoolVec tf(n * ids.size(), false);
  auto tf_iter = tf.begin();
  // Quantize each of the values to n bits 
  for (size_t i = 0; i < ids.size(); i++)
  {
    auto b = int2bin_bounded( ids[i] << n, n);
    for (int i = 0; i < n; i++, ++tf_iter)
    {
  
      *tf_iter = b.second[n - 1 - i];
                       
    } 
  }

  return tf; 
}

// [[Rcpp::export]]
Rcpp::List intersect_cells(Rcpp::List a, Rcpp::List b)
{
  
  return a;
}


// [[Rcpp::export]]
Rcpp::List eliasFanoCodingCpp(const IntegerVector& expression_vector) {
  
  
  
  std::deque<int> sparse_index;
  
  int i = 0;
  
  for(auto const& expr: expression_vector)
  {
    i++; 
    if (expr > 0)
    {
      sparse_index.push_back(i);
    }
  }

  if(sparse_index.empty())
  {
    // std::cerr << "Empty vector" << std::endl;
    return List::create();
  }
  unsigned int l = int(log2(expression_vector.size() / (float)sparse_index.size()) + 0.5) + 1;
  std::vector<int> ids(sparse_index.begin(), sparse_index.end());
  

  int prev_indexH = 0;
  BoolVec L(l*ids.size(), false);
  std::vector<bool> H;

  auto l_iter = L.begin();
  for (auto const& expr : ids)
  {
    auto c = int2bin_bounded(expr, l);
    
    for( int i = 0; i < l; i++, ++l_iter)
    {
      *l_iter = c.second[i];
    }
    // L.insert(L.end(), c.begin(), c.begin() + l); //append the lower bits to L
    
    //Use a unary code for the high bits
    unsigned int upper_bits = (expr >> l);
    unsigned int m =  H.size() + upper_bits - prev_indexH + 1;
    prev_indexH = upper_bits;
    H.resize(m, false);
    H[m - 1] = true;
  }
  return List::create(_["H"] = H, _["L"] = L, _["l"] = l);
}




// [[Rcpp::export]]
NumericVector eliasFanoDecodingCpp(BoolVec H, BoolVec L, int l)
{
  NumericVector ids(L.size() / l);
  unsigned int cells_found = 0;
  unsigned int H_i = 0;
  auto prev_it = H.begin() - 1;
  int i = 0;
  for (auto true_it = std::find(H.begin(), H.end(), true); 
       true_it != H.end() && i < ids.size(); 
       true_it = std::find(true_it + 1, H.end(), true), ++i)
  {
    auto offset  = std::distance(prev_it, true_it);
    prev_it = true_it;
    H_i += offset - 1;
    int id = H_i << l;
    for (unsigned short k = 0; k < l; ++k)
    {
      id |= (L[(i * l) + k] << k); 
    }
    ids[i] = id;
  }
  return ids;
}
