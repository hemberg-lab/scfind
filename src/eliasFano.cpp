#include <Rcpp.h>
#include <iostream>
#include <bitset>
#include <algorithm>
#include <iterator>
#include <set>
#include <cmath>
// #include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include "scfind_types.h"


// the bits used for the encoding
#define BITS 32


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
inline BitSet32 int2bin(unsigned int id)
{
  return make_pair(__builtin_clz(id), int2bin_core(id));
}


class EliasFanoDB;
RCPP_EXPOSED_CLASS(EliasFanoDB)


class EliasFanoDB
{
  
  typedef struct
  {
    BoolVec H;
    BoolVec L;
    int l;
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
    int l = ef.l;

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

  std::vector<int> eliasFanoDecoding(const EliasFano& ef)
  {
    
    std::vector<int> ids(ef.L.size() / ef.l);

    // This step inflates the vector by a factor of 8
    std::vector<char> H;
    H.reserve(ef.H.size());
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

  // constructor
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
    std::vector<std::string> gene_names;
    gene_names.reserve(genes.size());
    for(Rcpp::CharacterVector::iterator gene_it = genes.begin(); gene_it != genes.end(); ++gene_it)
    {
      gene_names.push_back(Rcpp::as<std::string>(*gene_it));
    }
     
    for(unsigned int gene_row = 0; gene_row < gene_matrix.nrow(); ++gene_row)
    {
      const Rcpp::NumericVector& expression_vector = gene_matrix(gene_row, Rcpp::_);
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
      if ( i == 0)
      {
        warnings++;
        continue;
      }
      std::vector<int> ids(sparse_index.begin(), sparse_index.end());
      insertToDB(eliasFanoCoding(ids, expression_vector.size()), gene_names[gene_row], cell_type);
    }
    return 0;
    //std::cerr << "Total Warnings: "<<warnings << std::endl;
  }

  int total_genes()
  {
    for(auto & d : metadata)
    {
      std::cout << d.first << ",";
      
    }
    std::cout << std::endl;
    return 0;
  }

  Rcpp::List queryGenes(const Rcpp::CharacterVector& gene_names)
  {
    Rcpp::List t;
    for (Rcpp::CharacterVector::const_iterator it = gene_names.begin(); it != gene_names.end(); ++it)
    {
      
      
      std::string gene_name = Rcpp::as<std::string>(*it);
      std::cout << gene_name << std::endl;
      
      //t.add(Rcpp::wrap(gene_names[i]), Rcpp::List::create());
      Rcpp::List cell_types;
      
      if (metadata.find(gene_name) == metadata.end())
      {
        
        std::cout << "Gene " << gene_name << " not found in the index " << std::endl;
        continue;
      }

      auto gene_meta = metadata[gene_name];
      for (auto const& dat : gene_meta)
      {

        std::vector<int> ids = eliasFanoDecoding(*(dat.second));
        cell_types[dat.first] = Rcpp::wrap(ids);
      }
      t[gene_name] = cell_types;
      // t.names() = gene_names;
    }

    return t;
  }
  
  size_t dbMemoryFootprint()
  {
    size_t bytes = 0;
    for(auto & d : ef_data)
    {
      bytes += int((d.H.size() / 8) + 1);
      bytes += int((d.L.size() / 8) + 1);
      bytes += 4 + 8;
    }

    std::cout << "Raw elias Fano Index size " << bytes/(1024*1024) << "MB" << std::endl;

    for(auto& d : metadata)
    {
      bytes += d.first.size();
      bytes += d.second.size() * 4;
      for (auto& ct : d.second)
      {
        bytes += ct.first.size();
      }
    }
    return bytes;
  }


  // And query
  Rcpp::List findCellTypes(const Rcpp::CharacterVector& gene_names)
  {
    
    std::map<std::string, std::set<std::string> > cell_types;
    std::vector<std::string> genes;
    for (Rcpp::CharacterVector::const_iterator it = gene_names.begin(); it != gene_names.end(); ++it)
    {
      std::string gene_name = Rcpp::as<std::string>(*it);
      bool empty_set = false;
      // check if gene exists in the database
      auto db_it = metadata.find(gene_name);
      if ( db_it == metadata.end())
      {
        std::cout << gene_name << " is ignored, not found in the index"<< std::endl;
        continue;
      }

      // iterate cell type
      for (auto const& ct_it : db_it->second)
      {
        if(cell_types.find(ct_it.first) == cell_types.end())
        {
          cell_types[ct_it.first] = std::set<std::string>();
        }
        cell_types[ct_it.first].insert(gene_name);
      }
      genes.push_back(gene_name);
      
    }
    
    std::cout << "Proceeding with "<< cell_types.size() << " cell types " << std::endl;
//    for ( auto ct : cell_types)
//    {
//      std::cout << ct.first << ", " << ct.second.size() << std::endl;
//    }
    std::cout << std::endl;
    Rcpp::List t;
    
    for (auto const& ct : cell_types)
    {
      // TODO(fix) empty case?
      bool empty_set = false;
      if (ct.second.size() != genes.size())
      {
        continue;
      }
      std::vector<int> ef = eliasFanoDecoding(*(metadata[*(ct.second.begin())][ct.first]));
      std::set<int> int_cells(ef.begin(), ef.end());
      for (auto const& g : ct.second)
      {
        auto cells = eliasFanoDecoding(*(metadata[g][ct.first]));
        std::set<int> new_set;
        std::set_intersection(int_cells.begin(), int_cells.end(), cells.begin(), cells.end(), std::inserter(new_set, new_set.begin()));
        if(new_set.size() != 0)
        {
          int_cells = new_set;
        }
        else
        {
          empty_set = true;
          break;
        }
      } 
      if (!empty_set)
      {
        // std::vector<int> res(int_cells.begin(), int_c;
        t[ct.first] = Rcpp::wrap(int_cells);
      }
    }
    return t;
  }


  int dbSize()
  {
    std::cout << metadata.size() << "genes in the DB" << std::endl;
    return ef_data.size();
    
  }

  int sample(int index)
  {
    auto iter = metadata.begin();
    for(int i = 0; i < index; i++, ++iter);

    std::cout << "Gene: " << iter->first << std::endl;
    for(auto const& ct : iter->second)
    {
      std::cout << "Cell Type:" << ct.first << std::endl;
      auto v = eliasFanoDecoding(*(ct.second));
      for( auto& cell : v)
      {

        std::cout << cell << ", ";
      }
      std::cout << std::endl;
      
    }
    return 0;
    
  }

  std::vector<int> decode(int index)
  {
    if(index >= dbSize())
    {
      std::cerr << "Invalid index for database with size "<< dbSize() << std::endl;
      return std::vector<int>();
    }
    return eliasFanoDecoding(ef_data[index]);
  }
  Rcpp::List queryMarkerGenes(const Rcpp::CharacterVector& gene_names)
  {
    

  }

  int mergeDB(const EliasFanoDB& db)
  {    
    EliasFanoDB extdb(db);
 
    // Copy data on the running object
    // and update pointers
    // In order to maintain consistency
    for ( auto& gene: extdb.metadata)
    {
      for( auto& ct : gene.second)
      {
        ef_data.push_back(*ct.second);
        // Update with the new entry
        ct.second = &ef_data.back();
      }
    }

    for(auto const& gene : extdb.metadata)
    {
      // if gene does not exist in the index then insert the whole record from
      // the whole database
      if(metadata.find(gene.first) == metadata.end())
      {
        metadata[gene.first] = extdb.metadata[gene.first];
      }
      else
      {
        // Insert new cell types
        for(auto& ct: extdb.metadata[gene.first])
        {
          metadata[gene.first][ct.first] = ct.second;
        }
      }
    }
    return 0;
  }

};


RCPP_MODULE(EliasFanoDB)
{
  Rcpp::class_<EliasFanoDB>("EliasFanoDB")
    .constructor()
    .method("indexMatrix", &EliasFanoDB::encodeMatrix)
    .method("queryGenes", &EliasFanoDB::queryGenes)
    .method("dbSize", &EliasFanoDB::dbSize)
    .method("decode", &EliasFanoDB::decode)
    .method("mergeDB", &EliasFanoDB::mergeDB)
    .method("sample", &EliasFanoDB::sample)
    .method("genes", &EliasFanoDB::total_genes)
    .method("findCellTypes", &EliasFanoDB::findCellTypes)
    .method("dbMemoryFootprint", &EliasFanoDB::dbMemoryFootprint);
}




