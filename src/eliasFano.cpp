//#include <Rcpp.h>
#include <iostream>
#include <bitset>
#include <algorithm>
#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "scfind_types.h"
// the bits used for the encoding
#define BITS 32


using namespace Rcpp;

typedef std::pair<unsigned short, std::bitset<BITS> > BitSet32;
typedef std::vector<bool> BoolVec;


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
BitSet32 int2bin(unsigned int id)
{
  return make_pair(__builtin_clz(id), int2bin_core(id));
}





class EliasFanoDB
{
  
  typedef struct
  {
    BoolVec H;
    BoolVec L;
    unsigned char l;
  } EliasFano;
  
  typedef std::map<std::string, const EliasFano*> EliasFanoIndex;
  typedef std::map<std::string, EliasFanoIndex > CellTypeIndex;
  typedef std::deque<EliasFano> ExpressionMatrix;

 private:
  CellTypeIndex metadata;
  ExpressionMatrix ef_data;
  bool global_indices;
  int warnings;
  void insertToDB(int ef_index, const std::string& gene_name, const std::string& cell_type)
  {
    if(ef_index == -1)
    {
      // Something went wrong so do not do anything
      this->warnings++;
      return;
    }

    const EliasFano* ef = &(this->ef_data[ef_index]);
    
    if(metadata.find(gene_name) == metadata.end())
    {
      metadata[gene_name] = EliasFanoIndex();
    }

    metadata[gene_name][cell_type] = ef;
  }

  long eliasFanoCoding(const std::vector<int>& ids, int items) 
  {
    
    if(ids.empty())
    {
      return -1;
    }

    EliasFano ef;
    ef.l = int(log2(items / (float)ids.size()) + 0.5) + 1;
    unsigned char& l = ef.l;

    int prev_indexH = 0;
    ef.L.resize(l * ids.size(), false);

    BoolVec::iterator l_iter = ef.L.begin();
    for (auto expr = ids.begin(); expr != ids.end(); ++expr)
    {
      BitSet32 c = int2bin_bounded(*expr, l);
    
      for( int i = 0; i < l; i++, ++l_iter)
      {
        // TODO(Nikos) optimize this ? check if c.second[i] is false and THEN assign
        *l_iter = c.second[i];
      }
    
      //Use a unary code for the high bits
      unsigned int upper_bits = (*expr >> l);
      unsigned int m =  ef.H.size() + upper_bits - prev_indexH + 1;
      prev_indexH = upper_bits;

      ef.H.resize(m, false);
      ef.H[m - 1] = true;
    }
    ef_data.push_back(ef);
    // return the index of the ef_data in the deque
    return ef_data.size() - 1;
  }

  std::vector<int> eliasFanoDecodingCpp(const EliasFano& ef)
  {
    
    std::vector<int> ids(ef.L.size() / ef.l);

    // This step inflates the vector by a factor of 8
    std::vector<char> H(ef.H.size());
    H.insert(H.end(), ef.H.begin(), ef.H.end());
    
    unsigned int H_i = 0;
    // Warning: Very very dodgy I might want to replace this with a check in the loop
    auto prev_it = H.begin() - 1;
    int i = 0;
    for (auto true_it = std::find(H.begin(), H.end(), true); 
         true_it != H.end() && i < ids.size(); 
         true_it = std::find(true_it + 1, H.end(), true), ++i)
    {
      size_t offset  = std::distance(prev_it, true_it);
      prev_it = true_it;
      H_i += offset - 1;
      int id = H_i << ef.l;
      for (unsigned short k = 0; k < ef.l; ++k)
      {
        id |= (ef.L[(i * ef.l) + k] << k); 
      }
      ids[i] = id;
    }
    return ids;
  }



 public:
  EliasFanoDB(): global_indices(false), warnings(0)
  {
    
  }
  
  
  // long encodeSparseMatrix(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& cell_type, const arma::sp_mat& mat)
  // {
  //   int items = mat.n_cols;
  //   for(unsigned int gene_row = 0; gene_row < mat.n_rows; ++gene_row)
  //   {
  //     std::deque<int> sparse_index;

  //     arma::sp_mat::col_iterator end_mat = mat.end_col(gene_row);
  //     for(arma::sp_mat::col_iterator it = mat.begin_col(gene_row); it != end_mat; ++it )
  //     {
  //       sparse_index.push_back(it.col());
  //     }
  //     std::vector<int> ids(sparse_index.begin(), sparse_index.end());
  //     insertToDB(gene_names[gene_row], cell_type, eliasFanoCoding(ids, items));
  //   }
    
  //   std::cerr << "Total Warnings: " << warnings << std::endl;
  // }

  long encodeMatrix(const std::string& cell_type, const Rcpp::NumericMatrix& gene_matrix)
  {
    int items = gene_matrix.ncol();
    Rcpp::CharacterVector genes = Rcpp::rownames(gene_matrix);
    std::vector<std::string> gene_names(genes.size());
    for(Rcpp::CharacterVector::iterator gene_it = genes.begin(); gene_it != genes.end(); ++gene_it)
    {
      gene_names.push_back(Rcpp::as<std::string>(*gene_it));
    }
     
    for(unsigned int gene_row = 0; gene_row < gene_matrix.nrow(); ++gene_row)
    {
      const Rcpp::NumericVector& expression_vector = gene_matrix(gene_row, _);
      std::deque<int> sparse_index;
      int i = 0;
      for (Rcpp::NumericVector::const_iterator expr = expression_vector.begin(); expr != expression_vector.end(); ++expr)
      {
        i++; 
        if (*expr > 0)
        {
          sparse_index.push_back(i);
        }
      }
      std::vector<int> ids(sparse_index.begin(), sparse_index.end());
      insertToDB(eliasFanoCoding(ids, items), gene_names[gene_row], cell_type);
    }

    std::cerr << "Total Warnings: "<<warnings << std::endl;
  }

  Rcpp::List queryAllGenes(const Rcpp::CharacterVector& gene_names)
  {
    //std::map<std::string,std::vector<int>>
    //if(){}

  }

  Rcpp::List queryMarkerGenes(const Rcpp::CharacterVector& gene_names)
  {
    

  }


};


using namespace Rcpp;
RCPP_EXPOSED_CLASS(EliasFanoDB)
RCPP_MODULE(EliasFanoClass){
  class_<EliasFanoDB>("EliasFanoDB")
    .constructor("Initializes object")
    // .method("indexSparseMatrix", &EliasFanoDB::encodeSparseMatrix, "Encodes sparse Matrix")
    .method("indexMatrix", &EliasFanoDB::encodeMatrix, "Encodes regular Matrix");
    }




