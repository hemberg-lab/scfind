#pragma once

#include <Rcpp.h>
#include <vector>
#include <set>

#include "typedef.h"

class EliasFanoDB;

class QueryScore
{
public:
  typedef struct 
  {
    double tfidf;
    size_t index;
    int cartesian_product_sets;
    int support_in_datasets;
  } GeneScore;

  friend class EliasFanoDB;
  std::map<std::string, GeneScore> genes;
  std::unordered_map<CellID , std::pair<std::vector<double>, int> > tfidf;
  QueryScore();
  void cell_type_relevance(const EliasFanoDB&, const Rcpp::List&, const std::set<std::string>&);
  float cell_tfidf(const EliasFanoDB&, const std::set<std::string>&);
  void estimateExpression(const Rcpp::List& gene_results, const EliasFanoDB& db, const Rcpp::CharacterVector& datasets);
  unsigned int geneSetCutoffHeuristic(const float = 0.5);
  int calculate_cell_types(const std::set<std::string>&gene_set);
};
