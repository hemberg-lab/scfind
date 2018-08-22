#pragma once
#include <Rcpp.h>
#include <iostream>

#include <algorithm>

#include <utility>
#include <map>
#include <vector>
#include <set>
#include <unordered_map>


#include "functions.h"

#define SERIALIZATION_VERSION 5





// template<typename T>
// struct TypeCompare
// {
//   using is_transparent = std::true_type;
//   bool operator()(const T& lhs, const T& rhs) const
//   {
//     return lhs.name < rhs.name;
//   }
  
//   bool operator()(const T& lhs, const std::string& name) const
//   {
//     return lhs.name < name;
//   }
  
//   bool operator()(const std::string& name, const T& rhs) const
//   {
//     return name < rhs.name;
//   }
// };



typedef struct
{
  BoolVec H;
  BoolVec L;
  int l;
  float idf; // tfidf
  Quantile expr;
} EliasFano;

typedef struct
{
  int gene;
  int cell_type;
  int index;
} IndexRecord;


class EliasFanoDB;
RCPP_EXPOSED_CLASS(EliasFanoDB)

typedef int EliasFanoID;
typedef int CellTypeID;



class GeneMeta
{
public:
  int total_reads;
  GeneMeta();
  void merge(const GeneMeta& other);
};

class CellMeta
{
public:                                               
  int reads;
  int features;
  CellMeta();
};

class CellID
{
public:
  CellTypeID cell_type;
  int cell_id;
  CellID(CellTypeID, int);
  bool operator==(const CellID& obj) const
  {
    return (obj.cell_type == cell_type) && (obj.cell_id == cell_id);

  }
};

namespace std
{
  template<>
  struct hash<CellID>
  {
    inline size_t operator()(const CellID& cid) const
    {
      return hash<CellTypeID>()(cid.cell_type) ^ hash<int>()(cid.cell_id);
    }
  };

}


class QueryScore
{
public:
  friend class EliasFanoDB;
  int cells_in_query;
  int cell_types_in_query;
  double query_score;
  std::map<std::string, double> gene_scores;
  std::unordered_map<CellID , std::vector<double> > tfidf;
  

  
  QueryScore();
  void reset();
  void cell_type_relevance(const EliasFanoDB&, const Rcpp::List&, const std::set<std::string>&);
  void cell_tfidf(const EliasFanoDB&, const std::set<std::string>&);
  void estimateExpression(const Rcpp::List& gene_results, const EliasFanoDB& db);
};




typedef struct
{
  std::string name;
  int total_cells;
} CellType;


class EliasFanoDB
{
 public:
  
  typedef std::string GeneName;
  typedef std::unordered_map<CellTypeID, EliasFanoID> GeneContainer;
  typedef std::map<GeneName, GeneContainer> GeneExpressionDB;

  typedef std::map<GeneName, GeneMeta> GeneIndex;
  
  typedef std::unordered_map<CellID, CellMeta> CellIndex;
  
  typedef std::string CellTypeName;
  typedef std::unordered_map<CellTypeName, CellTypeID> CellTypeIndex;
  typedef std::deque<EliasFano> ExpressionMatrix;
  
 // private:
  
  GeneExpressionDB index;
  CellIndex cells;
  CellTypeIndex cell_types;
  
  std::deque<CellType> inverse_cell_type;
  
  GeneIndex genes;
  
  ExpressionMatrix ef_data;

  unsigned int total_cells;
  
  unsigned char quantization_bits;
  

  EliasFanoDB();
  
  bool global_indices;
  
  int warnings;
  
  void dumpGenes();

  void clearDB();
 
  const EliasFano& getEntry(const GeneName& gene_name, const CellTypeName& cell_type) const;
 
  int loadByteStream(const Rcpp::RawVector& stream);

  Rcpp::RawVector getByteStream();

  long eliasFanoCoding(const std::vector<int>& ids, const Rcpp::NumericVector& values);
  
  std::vector<int> eliasFanoDecoding(const EliasFano& ef);

  // This is invoked on slices of the expression matrix of the dataset 
  long encodeMatrix(const std::string& cell_type_name, const Rcpp::NumericMatrix& gene_matrix);

  Rcpp::List total_genes();

  int getTotalCells() const;

  Rcpp::NumericVector getCellTypeSupport(Rcpp::CharacterVector& cell_types);
  
  std::vector<double> getQuantizedExpressionLevels(const std::string& gene_name, const std::string& cell_type);
  
  Rcpp::List queryGenes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& datasets_active);
  
  size_t dataMemoryFootprint();

  size_t dbMemoryFootprint();


  // And query
  Rcpp::List findCellTypes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& datasets_active);


  // TODO(Nikos) this function can be optimized.. It uses the native quering mechanism
  // that casts the results into native R data structures
  Rcpp::DataFrame findMarkerGenes(const Rcpp::CharacterVector& gene_list, const Rcpp::CharacterVector datasets_active, unsigned int min_support_cutoff);

  std::map<std::string, std::vector<int> > intersect_cells(std::set<std::string> gene_set, Rcpp::List genes_results) const ;

  int dbSize();

  int sample(int index);

  std::vector<int> decode(int index);
  
  int insertNewCellType(const CellType& cell_type);
  
  int mergeDB(const EliasFanoDB& db);


};
