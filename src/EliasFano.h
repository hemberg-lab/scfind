#pragma once
#include <Rcpp.h>
#include <iostream>

#include <algorithm>

#include <utility>
#include <map>
#include <vector>
#include <set>



#include "functions.h"

#define SERIALIZATION_VERSION 4


class CellType 
{
 public:
  std::string name;
  int total_cells;
  size_t operator()() const
  {
    return std::hash<std::string>{}(name);
  }
};



struct CellTypeCompare
{
  using is_transparent = std::true_type;
  bool operator()(const CellType& lhs, const CellType& rhs) const
  {
    return lhs.name < rhs.name;
  }
  
  bool operator()(const CellType& lhs, const std::string& name) const
  {
    return lhs.name < name;
  }
  
  bool operator()(const std::string& name, const CellType& rhs) const
  {
    return name < rhs.name;
  }
};


typedef float IDFtype;


typedef struct
{
  BoolVec H;
  BoolVec L;
  int l;
  IDFtype idf; // tfidf
  Quantile expr;
} EliasFano;

typedef struct
{
  int gene;
  int cell_type;
  int index;
} IndexRecord;


struct Cell_ID
{
  unsigned int num;
  const int cell_type;
  // Hashing
  size_t operator()() const
  {
    return std::hash<int>{}(cell_type) ^ std::hash<unsigned int>{}(num);
  }
  
  size_t operator==(const struct Cell_ID& obj) const 
  {
    return (num == obj.num) && (cell_type == obj.cell_type);
  }
};


typedef struct Cell_ID CellID;

namespace std
{
  template<>
  struct hash<CellID>
  {
    size_t operator()(const CellID& obj) const
    {
      // Return overloaded operator
      return obj();
    }
  };

  template<>
  struct hash<CellType>
  {
    size_t operator()(const CellType& obj) const
    {
      return obj();
    }
  };
  
}




class EliasFanoDB;
RCPP_EXPOSED_CLASS(EliasFanoDB)

typedef int EliasFanoID;
typedef int CellTypeID;
typedef std::string Gene;
typedef struct
{
  int reads;
  int feature;
}Cell;




class EliasFanoDB
{

 public:
  typedef std::unordered_map<CellTypeID, EliasFanoID> EliasFanoIndex;
  // gene -> cell type -> eliasFano
  typedef std::unordered_map<Gene, EliasFanoIndex > CellTypeIndex;
  typedef std::deque<EliasFano> ExpressionMatrix;
  
  // Store the gene metadata, gene support in cells at the index
  typedef std::map<std::string, unsigned int> GeneIndex;

 // private:
  CellTypeIndex metadata;
  ExpressionMatrix ef_data;
  std::map<CellType, int, CellTypeCompare> cell_types_id;
  std::deque<CellType> inverse_cell_type;
  std::map<CellType, std::vector<Cell> > cells;
  
  GeneIndex gene_counts;
  unsigned int total_cells;
  
  unsigned char quantization_bits;
  

  EliasFanoDB();
  
  
  
  bool global_indices;
  int warnings;
  
  void dumpGenes();

  void clearDB();


  int loadByteStream(const Rcpp::RawVector& stream);

  Rcpp::RawVector getByteStream();
  
  

  void insertToDB(int ef_index, const std::string& gene_name, const CellType& cell_type);
   

  long eliasFanoCoding(const std::vector<int>& ids, const Rcpp::NumericVector& values);
  

  std::vector<int> eliasFanoDecoding(const EliasFano& ef);


  // This is invoked on slices of the expression matrix of the dataset 
  long encodeMatrix(const std::string& cell_type_name, const Rcpp::NumericMatrix& gene_matrix);

  Rcpp::List total_genes();

  Rcpp::NumericVector getTotalCells();

  Rcpp::NumericVector getCellTypeSupport(Rcpp::CharacterVector& cell_types);
  

  Rcpp::List queryGenes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& datasets_active);
  
  size_t dataMemoryFootprint();

  size_t dbMemoryFootprint();


  // And query
  Rcpp::List findCellTypes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& datasets_active);


  // TODO(Nikos) this function can be optimized.. It uses the native quering mechanism
  // that casts the results into native R data structures
  Rcpp::DataFrame findMarkerGenes(const Rcpp::CharacterVector& gene_list, const Rcpp::CharacterVector datasets_active, unsigned int min_support_cutoff);

  std::map<std::string, std::vector<int> > intersect_cells(std::set<std::string> gene_set, Rcpp::List genes_results);

  int dbSize();

  int sample(int index);

  std::vector<int> decode(int index);
  
  int insertNewCellType(const CellType& cell_type);
  
  int mergeDB(const EliasFanoDB& db);


};
