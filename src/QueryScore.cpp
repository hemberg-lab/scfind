#include "QueryScore.h"
#include "typedef.h"
#include "utils.h"
#include "EliasFano.h"



QueryScore::QueryScore() : cells_in_query(0), query_score(0)
{}

void QueryScore::cell_type_relevance(const EliasFanoDB& db, const Rcpp::List& genes_results, const std::set<std::string>& gene_set)
{
  // for a specific gene set get the intersection of the cells as a map
  std::map<std::string, std::vector<int> > ct_map = db.intersect_cells(gene_set, genes_results);  
  
  // Remove all empty records
  for (auto it = ct_map.begin(); it != ct_map.end();)
  {
    it = it->second.empty() ? ct_map.erase(it) : ++it;
  }

  for (auto const& _ct : ct_map)
  {
    // const std::string& ct = _ct.first;
    auto cct =  db.cell_types.find(_ct.first);
    
    double tfidf = 0;
    for (auto const&  g : gene_set)
    {
      auto dbg = db.index.find(g);
      auto gct = dbg->second.find(cct->second);
      const auto& ef = db.ef_data[gct->second];
      tfidf += ef.idf;
    }
    query_score += tfidf * log(_ct.second.size());
    cells_in_query += _ct.second.size();
  }
  // Get the mean ct size
  query_score /= ct_map.size();
  cell_types_in_query = ct_map.size();
}

void QueryScore::reset()
{
  this->query_score = 0;
  this->cells_in_query = 0;
}

void QueryScore::estimateExpression(const Rcpp::List& gene_results, const EliasFanoDB& db, const Rcpp::CharacterVector& datasets, bool console_message = false)
{

  Rcpp::Rcout << "calculating tfidf for the reduced expression matrix... " << std::endl;


  // Store temporarily the strings so we can insert those in the map
  const auto& tmp_strings = Rcpp::as<std::vector<std::string>>(gene_results.names());
  const auto tmpl_cont = std::vector<double>(tmp_strings.size(), 0);
  int total_cells_in_universe = db.getTotalCells(datasets);

  
  Rcpp::IntegerVector gene_support = db.totalCells(gene_results.names(), datasets);

  std::vector<int> gs = Rcpp::as<std::vector<int>>(gene_support);
  
  std::vector<std::string> gs_names = Rcpp::as<std::vector<std::string>>(gene_support.names());
  std::vector<std::pair<std::string, int>> gs_pairs;
  gs_pairs.reserve(gs_names.size());
  std::transform(
    gs.begin(), 
    gs.end(), 
    gs_names.begin(),
    std::back_inserter(gs_pairs),
    [](const int& support, const std::string& gene_name){
      return std::make_pair(gene_name, support);
    });
  
  // Iterate through the genes to calculate the tfidf matrix
  // and do the cutoff estimation at once (for performance reasons)
  
  // for the cutoff estimation each gene is assigned a score .
  // That way we can estimate the distribution of the input gene list 
  // and do more accurate cutoff estimations

  // Build the reduced expression matrix
    
  bool estimate_cutoff = tmp_strings.size() > 7 ? true : false;
    
  for (size_t gene_row = 0; gene_row < tmp_strings.size(); ++gene_row)
  {

    const std::string& gene = tmp_strings[gene_row];
    // Rcpp::Rcerr << "Gene: " << gene << std::endl;
    float gene_idf = total_cells_in_universe / ((float)db.genes.at(gene).total_reads);
    // get the current score of the gene
    GeneScore g = {0, gene_row, 0};
    double& gene_score = this->genes.insert(std::make_pair(gene, g)).first->second.tfidf;
    
    const Rcpp::List& cts = gene_results[gene];
    const auto ct_names = Rcpp::as<std::vector<std::string>>(cts.names());
    for(auto const& cell_type : ct_names)
    {
      const Rcpp::IntegerVector& expr_indices  = cts[cell_type];
      const auto ctid_it = db.cell_types.find(cell_type);
      CellTypeID ct_id = ctid_it->second;
      std::vector<double> expr_values = decompressValues(db.getEntry(gene, cell_type).expr, db.quantization_bits);
      if (expr_values.size() != (unsigned int)expr_indices.size())
      {
        Rcpp::Rcerr << "Corrupted DB!" << std::endl;
      }
      int expr_index = 0;
      for (auto const& cell_id : expr_indices)
      {
        CellID cell(ct_id, cell_id);
        
        // insert if it does not exist
        auto ins_res = tfidf.insert(std::make_pair(cell, std::make_pair(tmpl_cont, 0)));
        
        // assign the decompressed expression vector to the data structure
        std::vector<double>& tfidf_vec = ins_res.first->second.first;
        
        // increase gene_support
        auto& gene_support = ins_res.first->second.second;
        gene_support++;
        
        // tfidf calculation ( the expression value , the total reads of that cell and the gene transcript abundance)
        tfidf_vec[gene_row] = (expr_values[expr_index++] / db.cells.at(cell).reads) * gene_idf;

        console_message == true ? Rcpp::Rcout << "gene " << gene << " cell type " << cell_type << " " << db.cells.at(cell).reads << std::endl : Rcpp::Rcout << ""; 

        gene_score += tfidf_vec[gene_row];
        
      }
    }
  }
    
  if (estimate_cutoff)
  {
    // iterate through cells for (gene cutoff)
    std::vector<int> genes_subset(this->genes.size(), 0);
    for (auto const& c : this->tfidf)
    {
      size_t i = 0;
      for (auto const& v : c.second.first)
      {
        genes_subset[i++] += v > 0 ? c.second.second - 1 : 0;
      }
    }

 
    for (auto& v : this->genes)
    {
      v.second.cartesian_product_sets = genes_subset[v.second.index];
      v.second.cartesian_product_sets /= this->genes.size();
      const auto& current_gene_name = v.first;
      int union_sum = std::accumulate(
        gs_pairs.begin(), 
        gs_pairs.end(), 
        0, 
        [&current_gene_name](const int& sum, const std::pair<std::string,int>& gs_pair){
          if(gs_pair.first == current_gene_name)
          {
            return sum;
          }
          return sum + gs_pair.second;
        });
      float mean_overlap = float(union_sum) / tfidf.size(); // normalize by cell size
      // how much the genes contribute to the overlap
      mean_overlap /= log(this->genes.size());
      v.second.cartesian_product_sets *= mean_overlap;
    }
    int i = 0;
    for (auto const& v : this->genes)
    {  
      Rcpp::Rcerr << "Cutoff proposed for gene " << v.first << ": " << v.second.cartesian_product_sets <<" with support " << gs_pairs[i++].second << std::endl;
    }
  }
  else
  {
    for (auto& v : this->genes)
    {
      v.second.cartesian_product_sets = 1;
    }

  }
}



void QueryScore::cell_tfidf(const EliasFanoDB& db, const std::set<std::string>& gene_set)
{
 
  this->query_score = 0;
  float min = genes[*(gene_set.begin())].tfidf;
  for(auto const& g : gene_set)
  {
    float tfidf = genes[g].tfidf;
    min = tfidf < min ? tfidf : min;
    this->query_score += tfidf;

  }
  this->query_score *= min;

}



int QueryScore::calculate_cell_types(const std::set<std::string>&gene_set)
{
  return 0;
}


