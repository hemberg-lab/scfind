#pragma once


#include <map>
#include <set>
#include <Rcpp.h>


#include "const.h"
#include "typedef.h"
#include "EliasFano.h"


std::map< EliasFanoDB::CellTypeName, std::map<int, Transaction> >  transposeResultToCell(const Rcpp::List& genes_results);


std::set<Pattern> FPGrowthFrequentItemsetMining(const Rcpp::List& genes_results, const unsigned int min_support_cutoff);
std::set<Pattern> exhaustiveFrequentItemsetMining(const Rcpp::List& gene_results, const unsigned int min_support_cutoff);


// helper functions

std::string str_join( const std::vector<std::string>& elements, const char* const separator);


inline std::bitset<BITS> int2bin_core(const unsigned int id)
{
  // we have a fixed BITS bitset so we have to offset 
  // the leading zeros whenever we need to access the elements
  return std::bitset<BITS>(id);
}


inline BitSet32 int2bin_bounded(const unsigned int id, unsigned int min_bit_length)
{
  // Check if number is larger than the desired size in bits
  unsigned short bit_set_length = id == 0 ? 0 : __builtin_clz(id);
  bit_set_length = bit_set_length > min_bit_length ? bit_set_length : min_bit_length;
  return make_pair(bit_set_length, int2bin_core(id));
}


// These functions work with 1-based indexes on arrays, 
// __builtin_clz has undefined behavior with 0
inline BitSet32 int2bin(unsigned int id)
{
  return make_pair(__builtin_clz(id), int2bin_core(id));
}

inline double normalCDF(const double& x, const double& mu, const double& sigma)
{
  // this is an inline function for the cdm of the normal distribution
  // it depends on the cmath library where it contains the erfc function
  // it return a value ranging from zero to one 
  
  return 1 - (0.5 * erfc( (x - mu)/ (sigma * M_SQRT1_2) ));
   
}

// Accepts a vector, transforms and returns a quantization logical vector
// This function aims for space efficiency of the expression vector
Quantile lognormalcdf(const std::vector<int>& ids, const Rcpp::NumericVector& v, unsigned int bits, bool raw_counts = true);


int byteToBoolVector(const std::vector<char> buf, std::vector<bool>& bool_vec);

int getSizeBoolVector(const std::vector<bool>& v);


std::vector<double> decompressValues(const Quantile& q, const unsigned char& quantization_bits);



