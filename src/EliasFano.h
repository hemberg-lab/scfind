#pragma once
#include <iostream>
#include <Rcpp.h>
#include <algorithm>
#include <utility>
#include <map>
#include <vector>
#include <set>
#include <unordered_map>

#include "const.h"
#include "typedef.h"

//' @export EliasFanoDB 
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
  int warnings;
  unsigned int total_cells;
  unsigned char quantization_bits;
  

  EliasFanoDB();

  EliasFanoDB(SEXPREC*&);
  
  void dumpGenes() const;

  void clearDB();
 
  int setQuantizationBits(const unsigned int value);

  unsigned int getQuantizationBits() const;

  const EliasFano& getEntry(const GeneName& gene_name, const CellTypeName& cell_type) const;
 
  int loadByteStream(const Rcpp::RawVector& stream);

  Rcpp::RawVector getByteStream() const;

  long eliasFanoCoding(const std::vector<int>& ids, const Rcpp::NumericVector& values);
  
  std::vector<int> eliasFanoDecoding(const EliasFano& ef) const; 

  int queryZeroGeneSupport(const Rcpp::CharacterVector&) const;

  // This is invoked on slices of the expression matrix of the dataset 
  long encodeMatrix(const std::string& cell_type_name, const Rcpp::NumericMatrix& gene_matrix);


  Rcpp::List total_genes() const;
  
  // Get a vector that represents support for a set of genes with respect to a specific dataset
  Rcpp::IntegerVector totalCells(const Rcpp::CharacterVector&, const Rcpp::CharacterVector&) const;
  
  // 
  Rcpp::CharacterVector getGenesInDB() const;
  
  int getTotalCells(const Rcpp::CharacterVector&) const;


  Rcpp::List geneSupportInCellTypes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector&) const;
  
  const CellType& getCellType(const CellTypeName& name ) const;

  const Rcpp::NumericMatrix getCellTypeMatrix(const CellTypeName& cell_type) const;

  int numberOfCellTypes(const Rcpp::CharacterVector&) const;

  int cellsInDB() const;
  
  CellTypeIndex getCellTypeIDs(const std::set<std::string>& datasets) const;
  
  Rcpp::NumericVector getCellTypeSupport(Rcpp::CharacterVector& cell_types);
  
  Rcpp::List queryGenes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& datasets_active) const;
  
  size_t dataMemoryFootprint() const;

  size_t quantizationMemoryFootprint() const;
  
  size_t dbMemoryFootprint() const;

  // And query
  Rcpp::List findCellTypes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& datasets_active) const;
  Rcpp::List _findCellTypes(const std::vector<std::string>& gene_names, const std::vector<CellTypeName>& cell_types_bg) const;


  // TODO(Nikos) this function can be optimized.. It uses the native quering mechanism
  // that casts the results into native R data structures
  Rcpp::DataFrame findMarkerGenes(const Rcpp::CharacterVector& gene_list, const Rcpp::CharacterVector datasets_active, bool exhaustive = false, const int user_cutoff = -1) const;
  

  Rcpp::DataFrame _findCellTypeMarkers(const Rcpp::CharacterVector& cell_types, 
                                       const Rcpp::CharacterVector& background, 
                                       const std::vector<GeneName>&,
                                       int mode = ALL) const;


  
  Rcpp::DataFrame findCellTypeMarkers(const Rcpp::CharacterVector& cell_types, 
                                      const Rcpp::CharacterVector& background) const;

  Rcpp::DataFrame evaluateCellTypeMarkers(const Rcpp::CharacterVector& cell_types, 
                                          const Rcpp::CharacterVector& gene_set, 
                                          const Rcpp::CharacterVector& background);

  Rcpp::DataFrame evaluateCellTypeMarkersAND(const Rcpp::CharacterVector& cell_types, 
                                          const Rcpp::CharacterVector& gene_set, 
                                          const Rcpp::CharacterVector& background);
  

  
  
  std::map<GeneName, CellTypeMarker> _cellTypeScore(const std::string& cell_type, const std::vector<std::string>& universe, const std::vector <GeneName>&, int mode = ALL) const;
  
  const std::set<std::string> _getValidCellTypes(std::vector<std::string> universe) const;


  const std::vector<CellTypeName> _getCellTypes() const;
  const std::vector<CellTypeName> _getCellTypes(const std::vector<std::string>& datasets) const;
  
  const std::vector<CellTypeName> getCellTypes() const;
  
  Rcpp::List getCellMeta(const std::string&, const int&) const;

  Rcpp::List getCellTypeMeta(const std::string&) const;

  int dbSize() const;
  
  void dumpEFsize(int) const;
  int sample(int index) const;

  std::vector<int> decode(int index) const;
  
  int insertNewCellType(const CellType& cell_type);
  
  int mergeDB(const EliasFanoDB& db);
  

};


