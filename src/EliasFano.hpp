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

#define SERIALIZATION_VERSION 6





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
  int getSize() const
  {
    return L.size() / l;
  }
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
  int getReads()const 
  {
    return reads;
  }
  
  int getFeatures() const
  {
    return features;
  }
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
  typedef struct 
  {
    double tfidf;
    size_t index;
    int cartesian_product_sets;
  } GeneScore;
  friend class EliasFanoDB;
  int cells_in_query;
  int cell_types_in_query;
  double query_score;
  std::map<std::string, GeneScore> genes;
  std::unordered_map<CellID , std::pair<std::vector<double>, int> > tfidf;
  

  
  QueryScore();
  void reset();
  void cell_type_relevance(const EliasFanoDB&, const Rcpp::List&, const std::set<std::string>&);
  void cell_tfidf(const EliasFanoDB&, const std::set<std::string>&);
  void estimateExpression(const Rcpp::List& gene_results, const EliasFanoDB& db, const Rcpp::CharacterVector& datasets);
};




class CellType
{
public:
  std::string name;
  int total_cells;
  int getTotalCells()const 
  {
    return total_cells;
  }
};

typedef struct
{
  int tp;
  int fp;
  int tn;
  int fn;
  float inv_precision() const
  {
    return  (tp + fp) / float(tp);
  }
  float inv_recall() const
  {
    return (fn + tp) /  float(tp);
  }

  float recall() const
  {
    return 1 / inv_recall();
  }

  float precision() const
  {
    return 1 / inv_precision();
  }

  float f1() const
  {
    return 2/(inv_precision() + inv_recall());
  }
} CellTypeMarker;


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
  
  // Get a vector that represents support for a set of genes with respect to a specific dataset
  Rcpp::IntegerVector totalCells(const Rcpp::CharacterVector&, const Rcpp::CharacterVector&);
  
  // 
  Rcpp::CharacterVector getGenesInDB();
  
  int getTotalCells(const Rcpp::CharacterVector&) const;

  int numberOfCellTypes(const Rcpp::CharacterVector&) const;

  int cellsInDB() const;
  
  CellTypeIndex getCellTypeIDs(const std::set<std::string>& datasets) const;
  
  Rcpp::NumericVector getCellTypeSupport(Rcpp::CharacterVector& cell_types);
  
  std::vector<double> getQuantizedExpressionLevels(const std::string& gene_name, const std::string& cell_type);
  
  Rcpp::List queryGenes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& datasets_active);
  
  size_t dataMemoryFootprint();

  size_t quantizationMemoryFootprint();
  
  size_t dbMemoryFootprint();
  
  

  // And query
  Rcpp::List findCellTypes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& datasets_active);


  // TODO(Nikos) this function can be optimized.. It uses the native quering mechanism
  // that casts the results into native R data structures
  Rcpp::DataFrame findMarkerGenes(const Rcpp::CharacterVector& gene_list, const Rcpp::CharacterVector datasets_active, unsigned int min_support_cutoff);


  Rcpp::DataFrame _findCellTypeMarkers(const Rcpp::CharacterVector& cell_types, 
                                       const Rcpp::CharacterVector& background, 
                                       const std::vector<GeneName>&);

  Rcpp::DataFrame findCellTypeMarkers(const Rcpp::CharacterVector& cell_types, 
                                      const Rcpp::CharacterVector& background);

  Rcpp::DataFrame evaluateCellTypeMarkers(const Rcpp::CharacterVector& cell_types, 
                                          const Rcpp::CharacterVector& gene_set, 
                                          const Rcpp::CharacterVector& background);
  

  
  
  std::map<GeneName, CellTypeMarker> _cellTypeScore(const std::string& cell_type, const std::vector<std::string>& universe, const std::vector <GeneName>&) const;
  const std::set<std::string> _getActiveCellTypes(std::vector<std::string> universe) const;
    
  const std::vector<CellTypeName> _getCellTypes() const;
  
  
  Rcpp::List getCellMeta(const std::string&, const int&) const;

  Rcpp::List getCellTypeMeta(const std::string&) const;

  std::map<std::string, std::vector<int> > intersect_cells(std::set<std::string> gene_set, Rcpp::List genes_results) const ;

  int dbSize();
  
  void dumpEFsize(int);
  int sample(int index);

  std::vector<int> decode(int index);
  
  int insertNewCellType(const CellType& cell_type);
  
  int mergeDB(const EliasFanoDB& db);
  

};
