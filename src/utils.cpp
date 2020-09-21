
#include <vector>
#include "utils.h"
#include "typedef.h"
#include "fp_growth.h"

std::string str_join( const std::vector<std::string>& elements, const char* const separator)
{
  switch (elements.size())
  {
    case 0:
      return "";
    case 1:
      return *(elements.begin());
    default:
      std::ostringstream os; 
      std::copy(elements.cbegin(), elements.cend() - 1, std::ostream_iterator<std::string>(os, separator));
      os << *elements.crbegin();
      return os.str();
  }
}



// Accepts a vector, transforms and returns a quantization logical vector
// This function aims for space efficiency of the expression vector
Quantile lognormalcdf(const std::vector<int>& ids, const Rcpp::NumericVector& v, unsigned int bits, bool raw_counts)
{
  
  std::function<double(const double&)> expr_tran = raw_counts ? [](const double& x) {return log(x + 1);}: [](const double& x){return x;};  

  Quantile expr;
  expr.mu = std::accumulate(ids.begin(),ids.end(), 0.0, [&v, &expr_tran](const double& mean, const int& index){
                                                          return  mean + expr_tran(v[index - 1]);
                                                        }) / ids.size();
  expr.sigma = sqrt(
    std::accumulate(
      ids.begin(), 
      ids.end(), 
      0.0, 
      [&v, &expr, &expr_tran](const double& variance, const int& index){
        return pow(expr.mu - expr_tran(v[index - 1]), 2);
      }) / ids.size());
  // initialize vector with zeros
  expr.quantile.resize(ids.size() * bits, 0);
  //Rcpp::Rcerr << "Mean,std" << expr.mu << "," << expr.sigma << std::endl;
  //Rcpp::Rcerr << "ids size " << ids.size() << " v size " << v.size() << std::endl;
  int expr_quantile_i = 0;
  for (auto const& s : ids)
  {
    unsigned int t = round(normalCDF(expr_tran(v[s]), expr.mu, expr.sigma) * (1 << bits));
    std::bitset<BITS> q = int2bin_core(t);
    for (unsigned int i = 0; i < bits; ++i)
    {
      expr.quantile[expr_quantile_i++] = q[i];
    }
  }
  return expr;
}

float inverf(float x)
{
  float tt1, tt2, lnx, sgn;
  sgn = (x < 0) ? -1.0f : 1.0f;
   
  x = (1 - x)*(1 + x);        // x = 1 - x*x;uo
  lnx = logf(x);

  tt1 = 2 / (PI * 0.147) + 0.5f * lnx;
  tt2 = 1 / 0.147 * lnx;
   
  return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}


double lognormalinv(const double& p, const double& mu, const double& sigma)
{
  return exp((inverf(2*p - 1) * sqrt(2) * sigma) + mu);
}



std::vector<double> decompressValues(const Quantile& q, const unsigned char& quantization_bits)
{
  int vector_size = q.quantile.size() / quantization_bits;
  std::vector<double> result(vector_size,0);

  if(quantization_bits > 16)
  {
    Rcpp::Rcerr << "Too much depth in the quantization bits!" << std::endl;
  }
  std::vector<double> bins((1 << quantization_bits));
  double bins_size = bins.size();
  // min value
  
  for(int i = 0; i < bins_size; ++i)
  {
    double cdf = (i + 0.5) / bins_size;
    bins[i] = lognormalinv(cdf, q.mu, q.sigma);
    // Rcpp::Rcout << "Quantile:" << bins[i] << " for " << cdf << " " << q.mu << " " << q.sigma << std::endl;
    
  }
 
  for (size_t i = 0; i < result.size(); ++i)
  {
    int quantile = 0;
    for (size_t j = 0; j < quantization_bits; ++j)
    {
      quantile |= (1 << q.quantile[(quantization_bits * i) + j]);
    }
    result[i] = bins[quantile];
  }
  return result;
}




int byteToBoolVector(const std::vector<char> buf, std::vector<bool>& bool_vec)
{
  bool_vec.resize((buf.size() << 3), 0);
  int c = 0;
  for (auto const& b : buf)
  {
    for ( int i = 0; i < 8; i++)
    {
      bool_vec[c++] = ((b >> i) & 1);
    }
  }

  return 0;
    
}


int getSizeBoolVector(const std::vector<bool>& v)
{
  int size = v.size() / 8;
  if (v.size() % 8 != 0)
  {
    ++size;
  }
  return size;
}



std::map<EliasFanoDB::CellTypeName, std::map<int, Transaction> >  transposeResultToCell(const Rcpp::List& genes_results)
{
  std::map<EliasFanoDB::CellTypeName, std::map<int, Transaction> > cells;
  
  const std::vector<std::string> gene_names = Rcpp::as<std::vector<std::string>>(genes_results.names());
  // Start inversing the list to a cell level
  for (auto const& gene_name : gene_names)
  {
    // Gene hits contains the cell type hierarchy
    const Rcpp::List& gene_hits = genes_results[gene_name];
    const Rcpp::CharacterVector& cell_type_names = gene_hits.names();
    for (auto const& _ct : cell_type_names)
    {
      const auto ct = Rcpp::as<std::string>(_ct);
        
      if (cells.find(ct) == cells.end())
      {
        cells[ct] = std::map<int, Transaction>();
      }

      std::vector<unsigned int> ids  = Rcpp::as<std::vector<unsigned int> >(gene_hits[ct]);
      auto& cell_index = cells[ct];
      // For all the hits
      for (auto const& id : ids)
      {
        // search for the id in the cell type space!
        if (cell_index.find(id) == cell_index.end())
        {
          // insert new cell
          cell_index[id] = Transaction();
        }
        // Add gene hit in the cell
        cell_index[id].push_back(gene_name);
      }
    }
  }

  return cells;
}

// Wrapper function
std::set<Pattern> FPGrowthFrequentItemsetMining(const Rcpp::List& genes_results, const unsigned int min_support_cutoff)
{
  const std::map<EliasFanoDB::CellTypeName, std::map<int, Transaction> > result_by_celltype = transposeResultToCell(genes_results);
 
  // Find out how many cells are in the query
  const unsigned int cells_present = std::accumulate(result_by_celltype.begin(),
                                                     result_by_celltype.end(),
                                                     0 , 
                                                     [](const int& sum, const std::pair<EliasFanoDB::CellTypeName, std::map<int, Transaction> >& celltype){
                                                       return sum + celltype.second.size();
                                                     });
  
  Rcpp::Rcerr << "Query Results Transposed: found " << cells_present << " sets" << std::endl;
  
  // Collect all transactions for fp-growth
  std::vector<Transaction> transactions;
  transactions.reserve(cells_present);
  for (auto const & ct : result_by_celltype)
  {
    // Iterate Cells of cell type
    for (auto const & cl : ct.second)
    {
      // Maybe sort? Should be sorted
      // std::sort(cl.second.begin(), cl.second.end());
      if (cl.second.size() != 1)
      {
        transactions.push_back(std::vector<Item>(cl.second.begin(), cl.second.end()));
      }
    }
  }

  Rcpp::Rcerr << transactions.size() << " transactions" << std::endl;

  const FPTree fptree{transactions, min_support_cutoff};
  return fptree_growth(fptree);
}



class CellIDs
{
 public:
  EliasFanoDB::GeneName gene_name;
  std::deque<CellID> cell_ids;
  CellIDs(const EliasFanoDB::GeneName& gene_name) : gene_name(gene_name)
  {

  }
};

typedef std::vector< CellIDs > CellVectors;

int findAllGeneSets(const CellVectors& query_results, std::set<Pattern>& gene_sets, const unsigned int min_support_cutoff)
{

  int gene_number = query_results.size();
  unsigned long gene_limit = 1 << (gene_number + 1) ;
  
  if(gene_number > 32)
  {
    return 1;
  }
  
  Rcpp::Rcerr << "Starting Exhaustive search with " << gene_number << " genes" << std::endl;
 
  
  // mask is a set of genes that each bit set states the gene presence
  for (unsigned long mask = 1; mask < gene_limit; ++mask)
  {
    // auto& vector it->second;
    std::deque<CellID> current_cell_vector;
    bool set = false;
    bool fail = false;
    Pattern current_pattern;
    // https://stackoverflow.com/a/7774362/4864088
    for(unsigned char i = 0; i < gene_number; ++i)
    {
      if (char(mask >> i) & (0x01))
      {
        if (!set)
        {
          // on first gene set the set;
          current_cell_vector = query_results[i].cell_ids;
          current_pattern.first.insert(query_results[i].gene_name);
          current_pattern.second = current_cell_vector.size();
          set = true;
        }
        else
        {
          const auto& current_gene_cells = query_results[i].cell_ids;
          std::deque<CellID> intersection;
          std::set_intersection(
            current_gene_cells.begin(),
            current_gene_cells.end(), 
            current_cell_vector.begin(),
            current_cell_vector.end(),
            std::back_inserter(intersection)
                                );
          if (intersection.size() < min_support_cutoff)
          {
            // no reason to go further
            fail = true;
            break;
          }
          else
          {
            current_pattern.first.insert(query_results[i].gene_name);
            current_pattern.second = intersection.size();
            current_cell_vector = std::move(intersection);
          }
        }
      }
    }
    if (not fail)
    {
      gene_sets.insert(current_pattern);  
    }
    
  }
  return 0;

}


std::set<Pattern> exhaustiveFrequentItemsetMining(const Rcpp::List& genes_results, const unsigned int min_support_cutoff)
{
  
  std::set<Pattern> results;

  // Boiler Plate code
  const std::vector<std::string> gene_names = Rcpp::as<std::vector<std::string>>(genes_results.names());  
  std::unordered_map<EliasFanoDB::CellTypeName,CellTypeID> celltype_ids;
  CellVectors gene_cell_vector;
  
  // Start inversing the list to a cell level
  for (auto const& gene_name : gene_names)
  {
    
    CellIDs cell_vector(gene_name);
    auto& cells = cell_vector.cell_ids;

    // Gene hits contains the cell type hierarchy
    const Rcpp::List& gene_hits = genes_results[gene_name];
    const Rcpp::CharacterVector& cell_types_in_gene = gene_hits.names();
    for (auto const& celltype : cell_types_in_gene)
    {
      // Generate Cell Type ID 
      auto celltype_name = Rcpp::as<std::string>(celltype);
      CellTypeID ct_id = celltype_ids.insert(
        std::make_pair(
          celltype_name,
          celltype_ids.size()
                       )).first->second;
      
      std::vector<unsigned int> ids  = Rcpp::as<std::vector<unsigned int> >(gene_hits[celltype_name]);
      for (auto const& id : ids)
      {
        cells.push_back(CellID(ct_id, id));
      }
    }

    // (Optimization) if genes have not a minimum support then remove them
    if (not (cell_vector.cell_ids.size() < min_support_cutoff))
    {
      gene_cell_vector.push_back(cell_vector);
    }
    
  }
  for(auto& cell: gene_cell_vector)
  {
    auto& cell_vector = cell.cell_ids;
    std::sort(cell_vector.begin(), cell_vector.end());
  }

  findAllGeneSets(gene_cell_vector, results, min_support_cutoff);

  return results;
}


