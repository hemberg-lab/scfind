#include <cmath>
#include <iterator>
#include <numeric>
#include <functional>
#include <exception>
#include <stdexcept>


#include "EliasFano.h"
#include "Serialization.h"
#include "QueryScore.h"
#include "typedef.h"
#include "utils.h"


CellMeta::CellMeta() : reads(0), features(0)
{}


GeneMeta::GeneMeta() : total_reads(0)
{}

void GeneMeta::merge(const GeneMeta& other) 
{
  this->total_reads += other.total_reads;
}

CellID::CellID(CellTypeID ct, int cid) : cell_type(ct), cell_id(cid)
{}

void EliasFanoDB::clearDB()
{
  // Clear the database
  index.clear();
  cell_types.clear();
  inverse_cell_type.clear();
}

int EliasFanoDB::setQuantizationBits(unsigned int qbvalue)
{
 
  if (ef_data.empty() and qbvalue < 32)
  {
    quantization_bits = qbvalue;
  }
  else
  {
    Rcpp::Rcerr << "Quantized bits not set, DB not empty or qbvalue to high!" << std::endl;
    return 1;
  }

  if (qbvalue > 10)
  {
    Rcpp::Rcerr << "Setting to high value may be a performance hog in retrieving cell expression" << std::endl;
    
  }
  return 0;
}

unsigned int EliasFanoDB::getQuantizationBits() const
{
  return this->quantization_bits;
}


int EliasFanoDB::loadByteStream(const Rcpp::RawVector& stream)
{
  clearDB();
  SerializationDB ser;
  ser.loadByteStream(stream);
  ser.deserializeDB(*this);
  return 0;
}

Rcpp::RawVector EliasFanoDB::getByteStream()
{
  SerializationDB ser;
  return ser.getByteStream(*this);
  
} 

long EliasFanoDB::eliasFanoCoding(const std::vector<int>& ids, const Rcpp::NumericVector& values) 
{
    
  if(ids.empty())
  {
    return -1;
  }
  int items = values.size();

    
  EliasFano ef;
  ef.l = int(log2(items / (float)ids.size()) + 0.5) + 1;
  ef.idf = log2(items / (float)ids.size());
  int l = ef.l;

  int prev_indexH = 0;
  ef.L.resize(l * ids.size(), false);

  BoolVec::iterator l_iter = ef.L.begin();
  ef.expr = lognormalcdf(ids, values, this->quantization_bits);
  
    
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
  // Rcpp::Rcerr <<"New index" << ef_data.size() - 1 << std::endl;
  // return the index of the ef_data in the deque
  return ef_data.size() - 1;
}


std::vector<int> EliasFanoDB::eliasFanoDecoding(const EliasFano& ef) const 
{
    
  
  // This step inflates the vector by a factor of 8
  std::vector<char> H;
  std::vector<int> ids(ef.L.size() / ef.l);
  H.reserve(ef.H.size());
  H.insert(H.end(), ef.H.begin(), ef.H.end());
    

  unsigned int H_i = 0;
  // Warning: Very very dodgy I might want to replace this with a check in the loop
  auto prev_it = H.begin() - 1;
  size_t i = 0;
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



const CellType& EliasFanoDB::getCellType(const CellTypeName& name ) const
{
  
  auto id = this->cell_types.at(name);
  return  this->inverse_cell_type.at(id);
}


// TODO(Nikos) this needs refactoring
const Rcpp::NumericMatrix EliasFanoDB::getCellTypeMatrix(const CellTypeName& cell_type) const
{
  const CellType ct = getCellType(cell_type);
  const CellTypeID ct_id = this->cell_types.at(cell_type);
  std::vector<GeneName>  feature_names;
  
  // Feature number mill be the feature names size
  for (auto const& record : index)
  {
    auto rec_it = record.second.find(ct_id);
    if (rec_it != record.second.end())
    {
      feature_names.push_back(record.first);
    }
  }

  
  // Initialize matrix
  Rcpp::NumericMatrix mat(feature_names.size(), ct.total_cells);
  
  int qb = quantization_bits;
  // for the sparse expression  vector matrix get the indices and deconvolute the quantized values
  for (int row = 0; row < mat.nrow(); ++row)
  {
    // for each gene in the database extract the values
    const auto& rec = getEntry(feature_names[row], cell_type);
    const auto indices_val = eliasFanoDecoding(rec);
    const auto exp_val = decompressValues(rec.expr, qb);

    if (indices_val.size() != exp_val.size())
    {
      Rcpp::Rcerr << "not equal number of genes" << std::endl;
      continue;
    }
    Rcpp::NumericVector na_vec(ct.total_cells);
    auto exp_it = exp_val.begin();
    
    if (exp_val.size() != indices_val.size())
    {
      Rcpp::Rcerr << "Sparse vector representation mismatch" << std::endl;
      Rcpp::Rcerr << feature_names[row] << std::endl;
      continue;
    }
    
    for (auto const& index : indices_val)
    {
      na_vec[index - 1] = (*exp_it);
      ++exp_it;
    }
    

    mat(row, Rcpp::_) = na_vec;
  }
  

  Rcpp::rownames(mat) = Rcpp::wrap(feature_names);

  return mat;

}



const EliasFano& EliasFanoDB::getEntry(const GeneName& gene_name, const CellTypeName& cell_type) const 
{
  try
  {
    return this->ef_data.at(this->index.at(gene_name).at(this->cell_types.at(cell_type)));                       
  } 
  catch(const std::out_of_range& e)
  {
    Rcpp::Rcerr << e.what() << std::endl;
    auto g_it = index.find(gene_name);
    if (g_it == index.end())
    {
      Rcpp::Rcerr << gene_name << "Gene not found" << std::endl;
    }
    auto ct_it = this->cell_types.find(cell_type);

    if (ct_it == this->cell_types.end())
    {
      Rcpp::Rcerr << "Cell type"<<cell_type<<" not found in the database" << std::endl;
    }
    else
    {
      auto ef_it = g_it->second.find(ct_it->second);
      
      if(ef_it == g_it->second.end())
      {
        Rcpp::Rcerr << "Cell type "<< cell_type<<" not found for gene " << gene_name << std::endl;
        
      }
    }
    throw std::invalid_argument("Unable to retrieve entry from database");
  }
}



// constructor
EliasFanoDB::EliasFanoDB(): 
  warnings(0), 
  total_cells(0),
  quantization_bits(2)
{
    
}

EliasFanoDB::EliasFanoDB(SEXPREC*& obj)
{}



int EliasFanoDB::queryZeroGeneSupport(const Rcpp::CharacterVector& datasets) const
{
  int zs = 0;
  for(auto const& g : this->index)
  {
    auto cell_support = this->totalCells(Rcpp::wrap(g.first), datasets);
    // std::vector<int> cell_support =  Rcpp::as<std::vector<int>>();
    if (cell_support[0] == 0){
      zs++;
      Rcpp::Rcerr << "Gene " << g.first << " found no support with " <<  g.second.size() << " cell types"<< std::endl;
    }
    
  }
  return zs;
}


// This is invoked on slices of the expression matrix of the dataset 
long EliasFanoDB::encodeMatrix(const std::string& cell_type_name, const Rcpp::NumericMatrix& gene_matrix)
{
  Rcpp::CharacterVector cell_type_genes = Rcpp::rownames(gene_matrix);

  CellType cell_type;
  cell_type.name = cell_type_name;
  cell_type.total_cells = gene_matrix.ncol();
  
  int cell_type_id = insertNewCellType(cell_type);

  auto gene_names =  Rcpp::as<std::vector<std::string>>(cell_type_genes);
  
  // Increase the cell number present in the index
  this->total_cells += gene_matrix.ncol();


  // Store the metadata for the cell
  std::vector<CellMeta> current_cells(gene_matrix.ncol());
  


  // This can be scaled!
  for (int gene_row = 0; gene_row < gene_matrix.nrow(); ++gene_row)
  {
    
    // Collect indices of sparse vector for this gene
    const Rcpp::NumericVector& expression_vector = gene_matrix(gene_row, Rcpp::_);
    std::deque<int> sparse_index;
    int i = 0;
    for (Rcpp::NumericVector::const_iterator expr = expression_vector.begin(); expr != expression_vector.end(); ++expr)
    {
      i++; // 1 based indexing!
      
      if (*expr > 0)
      {
        // Not Thread safe as well
        current_cells[i - 1].reads += *expr;
        current_cells[i - 1].features++;
        // End of critical area
        
        sparse_index.push_back(i);
      }
    }
    if (sparse_index.empty())
    {
      continue;
    }
    
    // Insert the gene gene metadata database
    auto gene_it = this->genes.insert(std::make_pair(gene_names[gene_row], GeneMeta())).first;
    // Insert the gene in the index
    auto db_entry = this->index.insert(std::make_pair(gene_names[gene_row], GeneContainer())).first;
    
    // Cast it as a vector
    std::vector<int> ids(sparse_index.begin(), sparse_index.end());

    // Is total reads or non zero reads ?
    gene_it->second.total_reads += ids.size();
    auto ef_index = eliasFanoCoding(ids, expression_vector);
    // Insert the entry to index
    if (ef_index != -1)
    {
      db_entry->second.insert(std::make_pair(cell_type_id, ef_index));
    }
  }
  
  int i = 0; // 1 based indexing
  for (auto const& cell : current_cells)
  {
    if (cell.reads == 0)
    {
      Rcpp::Rcerr << "Vector of zeros detected for cell " << cell_type_name << " " << i << std::endl;
    }
    this->cells.insert({CellID(cell_type_id, ++i), cell});
  }
  
  return 0;
  //Rcpp::Rcerr << "Total Warnings: "<<warnings << std::endl;
}



const std::vector<EliasFanoDB::CellTypeName> EliasFanoDB::getCellTypes() const
{
  return this->_getCellTypes();
}

const std::vector<EliasFanoDB::CellTypeName> EliasFanoDB::_getCellTypes(const std::vector<std::string>& datasets) const
{
  auto cts = this->_getCellTypes();
  std::vector<EliasFanoDB::CellTypeName> results;
  results.reserve(cts.size());
  std::copy_if(cts.begin(), 
               cts.end(), 
               std::back_inserter(results), 
               [&datasets](const CellTypeName& ct){
                 auto it = ct.find_first_of(".");
                 if (it != std::string::npos)
                 {
                   return std::find(datasets.begin(), datasets.end(),ct.substr(0, it)) != datasets.end();
                 }
                 return false;
               });
  
  return results;

}
Rcpp::List EliasFanoDB::geneSupportInCellTypes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& datasets_active) const
{
  auto cell_types = this->_getCellTypes(Rcpp::as<std::vector<std::string>>(datasets_active));
  auto genes = Rcpp::as<std::vector<EliasFanoDB::GeneName>>(gene_names);
  Rcpp::List results;

  for (auto const& g : genes)
  {
    Rcpp::IntegerVector gene_results;
    for (auto const& ct : cell_types)
    {
      // Querying cell types
      int size = 0;
      try{
   
        const auto r = this->ef_data.at(this->index.at(g).at(this->cell_types.at(ct)));
        size = r.getSize();
      }catch(const std::out_of_range& e)
      {
        continue;
      }

      gene_results.push_back(size, ct);
      
    }
    results[g] = gene_results;
  }
  return results;
}

Rcpp::List EliasFanoDB::total_genes()
{
  Rcpp::List t; 
  for(auto & d : index)
  {
    t.push_back(Rcpp::wrap(d.first));
  }
  return t;
}

Rcpp::IntegerVector EliasFanoDB::totalCells(const Rcpp::CharacterVector& genes, 
                                            const Rcpp::CharacterVector& datasets_active) const
{
  Rcpp::IntegerVector t(genes.size(), 0);
  std::vector<std::string> datasets = Rcpp::as<std::vector<std::string>>(datasets_active);
  t.names() = genes;
  std::vector<std::string> str = Rcpp::as<std::vector<std::string>> (genes);
  
  // Building the inverse index for index cell_types
  std::unordered_map<CellTypeID, CellTypeName> inv_ct;
  for (auto const& ct : this->cell_types)
  {
    inv_ct[ct.second] = ct.first;
  }

  int i = 0;
  for (auto const& g : str)
  {
    auto git = this->index.find(g);
    if (git != this->index.end())
    {
      for (auto const& ct : git->second)
      {
        const std::string& ct_name = inv_ct[ct.first];
        std::string ct_dataset = ct_name.substr(0, ct_name.find("."));
        auto find_dataset = std::find(datasets.begin(), datasets.end(), ct_dataset);
        // check if the cells are in active datasets
        if (find_dataset == datasets.end())
        {
          continue;
        }
        t[i] += this->ef_data[ct.second].getSize();
      }
    }
  
    i++;
  }
  return t;
  
}

EliasFanoDB::CellTypeIndex EliasFanoDB::getCellTypeIDs(const std::set<std::string>& datasets) const 
{
  CellTypeIndex cts;
  for(auto const& ct : this->inverse_cell_type)
  {
    auto index = ct.name.find_last_of(".");
    if(index != std::string::npos)
    {
      auto dataset = ct.name.substr(0, index);
      const auto it = datasets.find(dataset);
      if (it != datasets.end())
      {
        cts[ct.name] = this->cell_types.find(ct.name)->second;
      }
    }

  }
  return cts;
}

int EliasFanoDB::cellsInDB() const
{
  return this->total_cells;
}

int EliasFanoDB::getTotalCells(const Rcpp::CharacterVector& datasets) const
{
  std::vector<std::string> act = Rcpp::as<std::vector<std::string> >(datasets);
  std::set<std::string> act_set(act.begin(), act.end());
  CellTypeIndex active_cell_types = getCellTypeIDs(act_set);
  int total_number_of_cells = 0;
  for (auto const& ct : active_cell_types)
  {
    total_number_of_cells += this->inverse_cell_type[ct.second].total_cells;
  }
  return total_number_of_cells;
}

int EliasFanoDB::numberOfCellTypes(const Rcpp::CharacterVector& datasets) const
{
  std::vector<std::string> act = Rcpp::as<std::vector<std::string> >(datasets);
  std::set<std::string> act_set(act.begin(), act.end());
  CellTypeIndex active_cell_types = getCellTypeIDs(act_set);
  
  return active_cell_types.size();
}


Rcpp::NumericVector EliasFanoDB::getCellTypeSupport(Rcpp::CharacterVector& cell_types)
{
  std::vector<std::string> cts = Rcpp::as<std::vector<std::string>>(cell_types);
  std::vector<int> ct_support;
  ct_support.reserve(cts.size());
  for (auto const& ct : cts)
  {
    // TODO(Nikos) fix this, otherwise we will get a nice segfault no error control
    auto cit = this->cell_types.find(ct);
    if(cit != this->cell_types.end())
    {
      ct_support.push_back(this->inverse_cell_type[cit->second].total_cells);
    }
    else
    {
      ct_support.push_back(0);
    }
  }
    
    
  return Rcpp::wrap(ct_support);
}
  

Rcpp::List EliasFanoDB::queryGenes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& datasets_active)
{
  Rcpp::List t;
  for (Rcpp::CharacterVector::const_iterator it = gene_names.begin(); it != gene_names.end(); ++it)
  {
      
    std::string gene_name = Rcpp::as<std::string>(*it);
    //t.add(Rcpp::wrap(gene_names[i]), Rcpp::List::create());
    Rcpp::List cell_types;
      
    if (index.find(gene_name) == index.end())
    {
        
      Rcpp::Rcout << "Gene " << gene_name << " not found in the index " << std::endl;
      continue;
    }
    std::vector<std::string> datasets = Rcpp::as<std::vector<std::string>>(datasets_active);
    const auto& gene_meta = index[gene_name];
    for (auto const& dat : gene_meta)
    {
      CellType current_cell_type = this->inverse_cell_type[dat.first];
        
      std::string dataset = current_cell_type.name.substr(0, current_cell_type.name.find("."));
      auto ct_find = std::find(datasets.begin(), datasets.end(), dataset);
        
      if (ct_find == datasets.end())
      {
        continue;
      }
      std::vector<int> ids = eliasFanoDecoding(ef_data[dat.second]);
      cell_types[current_cell_type.name] = Rcpp::wrap(ids);
    }
    t[gene_name] = cell_types;
      
  }

  return t;
}
  
size_t EliasFanoDB::dataMemoryFootprint()
{
  size_t bytes = 0;
  for(auto & d : ef_data)
  {
    bytes += int((d.H.size() / 8) + 1);
    bytes += int((d.L.size() / 8) + 1);
    bytes += int((d.expr.quantile.size() / 8) + 12);
      
  }
  bytes += ef_data.size() * 32; // overhead of l idf and deque struct
  return bytes;
}

size_t EliasFanoDB::quantizationMemoryFootprint()
{
  size_t bytes = 0;
  for (auto & d : ef_data)
  {
    bytes += int((d.expr.quantile.size() / 8) + 12);
  }
  bytes += ef_data.size() * 32; // overhead of l idf and deque struct
  return bytes;
}


size_t EliasFanoDB::dbMemoryFootprint()
{
  size_t bytes = dataMemoryFootprint();

  Rcpp::Rcout << "Raw elias Fano Index size " << bytes / (1024 * 1024) << "MB" << std::endl;
    

  // GeneIndex genes GeneExpressionDB
  for (auto const& d : index)
  {
    // One for each
    bytes += d.first.size() * 2;
    bytes += d.second.size() * 4;
  }

  bytes += index.size() * (sizeof(GeneMeta) + 8);
  
  // CellIndex cells
  bytes += cells.size() * (sizeof(CellID) + sizeof(CellMeta) + 4);
  
  // CellTypeIndex cell_types std::deque<cellType> inverse_cell_type
  for (auto const& c : cell_types)
  {
    bytes += (c.first.size()*2) + 4 + sizeof(CellType);
  }
  
  bytes += 16;
  
  return bytes;
}


Rcpp::CharacterVector EliasFanoDB::getGenesInDB()
{
  std::vector<std::string> gene_names;
  gene_names.reserve(this->genes.size());
  for(auto const& g : this->genes)
  {
    gene_names.push_back(g.first);
  }
  return Rcpp::wrap(gene_names);
}


Rcpp::List EliasFanoDB::findCellTypes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& datasets_active) const 
{
  std::vector<std::string> datasets = Rcpp::as<std::vector<std::string>>(datasets_active);
  std::vector<CellTypeName> cell_types_bg;
  for (auto const& ct : this->cell_types)
  {
    cell_types_bg.push_back(ct.first);
  }
  cell_types_bg.erase(std::remove_if(
    cell_types_bg.begin(), 
    cell_types_bg.end(), 
    [&datasets](const CellTypeName& ct_name){
      std::string ct_dataset = ct_name.substr(0, ct_name.find("."));
      return std::find(datasets.begin(), datasets.end(), ct_dataset) == datasets.end();
    }), cell_types_bg.end());
  
  return _findCellTypes(Rcpp::as<std::vector<std::string>>(gene_names), cell_types_bg);
}


// And query
Rcpp::List EliasFanoDB::_findCellTypes(const std::vector<GeneName>& gene_names, const std::vector<EliasFanoDB::CellTypeName>& cell_types_bg) const
{

  // Store the results here
  Rcpp::List t;
  std::vector<GeneName> genes(gene_names);

  // Remove genes not found in index
  genes.erase(std::remove_if(genes.begin(), genes.end(),[&](const GeneName& g){
                                                          auto is_missing = (index.find(g) == index.end());
                                                          if (is_missing)
                                                            Rcpp::Rcerr << g << " is ignored, not found in the index"<< std::endl;
                                                          return is_missing;
                                                        }), genes.end());
  
  // Get Cell types that have all the genes present
  std::vector<CellTypeName> cts = cell_types_bg;
  std::vector<const GeneContainer*> gene_set;
  for(auto const& g : genes)
  {
    gene_set.push_back(&(this->index.at(g)));
  }
  cts.erase(std::remove_if(cts.begin(), cts.end(),[&](const CellTypeName& ct){
                                                    CellTypeID cid = this->cell_types.at(ct);
                                                    for( auto const& g : gene_set)
                                                      if (g->find(cid) == g->end())
                                                        return true;
                                                    return false;
                                                  }), cts.end());

  

  // intersect cells
  for (auto const& ct : cts)
  {
    auto last_intersection = eliasFanoDecoding(getEntry(*(genes.begin()),ct));
    std::vector<int> curr_intersection;
    curr_intersection.reserve(last_intersection.size());
    for(std::size_t i = 1; i < genes.size(); ++i)
    {
      std::vector<int> cells = eliasFanoDecoding(getEntry(genes.at(i), ct));
      std::set_intersection(
        cells.begin(), 
        cells.end(), 
        last_intersection.begin(), 
        last_intersection.end(), 
        std::back_inserter(curr_intersection));
      std::swap(last_intersection,curr_intersection);
      curr_intersection.clear();
      curr_intersection.reserve(last_intersection.size());
      if(last_intersection.empty())
      {
        break;
      }
    }
    if(not last_intersection.empty())
      t[ct] = Rcpp::wrap(last_intersection);
      
  }
  return t;
}




// that casts the results into native R data structures
Rcpp::DataFrame EliasFanoDB::findMarkerGenes(const Rcpp::CharacterVector& gene_list, const Rcpp::CharacterVector datasets_active, const int user_cutoff = -1)
{
    
  std::vector<std::string> query;
  std::vector<double> query_scores;
  std::vector<double> query_tfidf;
  std::vector<int> query_cell_type_cardinality;
  std::vector<int> query_cell_cardinality;
  std::vector<int> query_gene_cardinality;
  

  
    
  // Perform an OR query on the database as a first step
  const Rcpp::List genes_results = queryGenes(gene_list, datasets_active);
  const std::map<CellTypeName, std::map<int, Transaction> > result_by_celltype = transposeResultToCell(genes_results);
  
  // Find out how many cells are in the query
  const unsigned int cells_present = std::accumulate(result_by_celltype.begin(),
                                               result_by_celltype.end(),
                                               0 , 
                                               [](const int& sum, const std::pair<CellTypeName, std::map<int, Transaction> >& celltype){
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
  

  QueryScore qs;
  qs.estimateExpression(genes_results, *this, datasets_active);

  unsigned int min_support_cutoff = 0;
  if (user_cutoff < 0)
  {
    min_support_cutoff = qs.geneSetCutoffHeuristic();
  }
  else
  {
    min_support_cutoff = user_cutoff;
  }

  // Run fp-growth algorithm
  Rcpp::Rcerr << "Running fp-growth tree with " << min_support_cutoff << " cutoff "<< std::endl;
  const FPTree fptree{transactions, min_support_cutoff};
  const std::set<Pattern> patterns = fptree_growth(fptree);
  Rcpp::Rcerr << "Found " << patterns.size() << " geneset patterns " << std::endl;
  
  // Iterate through the calculated frequent patterns
  for (auto const& item : patterns)
  {
    
    Rcpp::List gene_query;
    const auto& gene_set = item.first;
    
    // We do not care for queries with cardinality less than 2
    if (gene_set.size() < 2)
    {
      continue;
    }
    
    std::string view_string = str_join(std::vector<Item>(gene_set.begin(), gene_set.end()), ",");

    
    // cell_type_relevance
    query_cell_cardinality.push_back(item.second);
    query_tfidf.push_back(qs.cell_tfidf(*this, gene_set));
    
    // other fields
    query_gene_cardinality.push_back(gene_set.size());
    query.push_back(view_string);
  }

  // Dump the list
  return Rcpp::DataFrame::create(Rcpp::Named("Genes") = Rcpp::wrap(query_gene_cardinality),
                                 Rcpp::Named("Query") = Rcpp::wrap(query),
                                 Rcpp::Named("tfidf") = Rcpp::wrap(query_tfidf),
                                 Rcpp::Named("Cells") = Rcpp::wrap(query_cell_cardinality));
}





const std::set<std::string> EliasFanoDB::_getValidCellTypes(std::vector<std::string> universe) const
{
  std::set<std::string> active_cell_types;
  std::vector<std::string> db_cell_types(this->_getCellTypes());
  std::sort(universe.begin(), universe.end());
  std::sort(db_cell_types.begin(), db_cell_types.end());
  std::set_intersection(
    db_cell_types.begin(), 
    db_cell_types.end(), 
    universe.begin(), 
    universe.end(), 
    std::inserter(active_cell_types, active_cell_types.begin()));
  if (universe.size() != active_cell_types.size())
  {
    std::vector<std::string> cts_not_found;
    std::set_difference(
      universe.begin(), 
      universe.end(), 
      active_cell_types.begin(), 
      active_cell_types.end(), 
      std::back_inserter(cts_not_found));
    for(auto const& ct : cts_not_found)
    {
      Rcpp::Rcerr << "Ignoring cell type "<< ct <<" Not found in DB" << std::endl;
    }
  }

  return active_cell_types;

}


Rcpp::DataFrame EliasFanoDB::findCellTypeMarkers(const Rcpp::CharacterVector& cell_types, const Rcpp::CharacterVector& background) const
{
  std::vector<GeneName> gene_set;
  gene_set.reserve(this->genes.size());
  for (auto const& g : this->genes)
  {
    gene_set.push_back(g.first);
  }
  return _findCellTypeMarkers(cell_types, background, gene_set);
}




Rcpp::DataFrame EliasFanoDB::_findCellTypeMarkers(const Rcpp::CharacterVector& cell_types, const Rcpp::CharacterVector& background, const std::vector<EliasFanoDB::GeneName>& gene_set, int mode) const
{
  std::vector<std::string> 
    bk_cts(Rcpp::as< std::vector<std::string> >(background)),
    cts(Rcpp::as< std::vector<std::string> >(cell_types)),
    genes, df_cell_type;
  
  
  std::vector<int> tp, fp, tn, fn;
  std::vector<float> precision, recall, f1;
  for (auto const& ct : cts)
  {
    auto marker_genes = this->_cellTypeScore(ct, bk_cts, gene_set, mode);
    if (marker_genes.empty())
    {
      Rcpp::Rcerr << "Marker genes could not be found for cell type " << ct << std::endl;
      continue;
    }
    
    for (auto& t : marker_genes)
    {
      const CellTypeMarker& ctm = t.second;
      genes.push_back(t.first);
      df_cell_type.push_back(ct);
      tp.push_back(ctm.tp);
      fp.push_back(ctm.fp);
      tn.push_back(ctm.tn);
      fn.push_back(ctm.fn);
      precision.push_back(ctm.precision());
      recall.push_back(ctm.recall());
      f1.push_back(ctm.f1());
    }
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("cellType") = Rcpp::wrap(df_cell_type),
    Rcpp::Named("genes") = Rcpp::wrap(genes),
    Rcpp::Named("tp") = Rcpp::wrap(tp),
    // Rcpp::Named("tn") = Rcpp::wrap(tn),
    Rcpp::Named("fp") = Rcpp::wrap(fp),
    Rcpp::Named("fn") = Rcpp::wrap(fn),
    Rcpp::Named("precision") = Rcpp::wrap(precision),
    Rcpp::Named("recall") = Rcpp::wrap(recall),
    Rcpp::Named("f1") = Rcpp::wrap(f1)
                                 );
  
}


Rcpp::DataFrame EliasFanoDB::evaluateCellTypeMarkersAND(const Rcpp::CharacterVector& cell_types, 
                                                     const Rcpp::CharacterVector& gene_set,  
                                                     const Rcpp::CharacterVector& background)
{
  return _findCellTypeMarkers(cell_types, background, Rcpp::as<std::vector<GeneName>>(gene_set), AND);
}


Rcpp::DataFrame EliasFanoDB::evaluateCellTypeMarkers(const Rcpp::CharacterVector& cell_types, 
                                                     const Rcpp::CharacterVector& gene_set,  
                                                     const Rcpp::CharacterVector& background)
{
  return _findCellTypeMarkers(cell_types, background, Rcpp::as<std::vector<GeneName>>(gene_set), ALL);
}



std::map<EliasFanoDB::GeneName, CellTypeMarker> EliasFanoDB::_cellTypeScore(const std::string& cell_type, const std::vector<std::string>& universe, const std::vector<EliasFanoDB::GeneName>& gene_names, int mode) const
{
  auto ct_it = this->cell_types.find(cell_type);
  if ( ct_it == this->cell_types.end())
  {
    Rcpp::Rcerr <<"Cell type "<< cell_type << " not found. exiting..." << std::endl;
    return std::map<EliasFanoDB::GeneName, CellTypeMarker>();
  }

  const CellTypeID cell_type_id = ct_it->second;
  int total_cells_in_ct = this->inverse_cell_type[cell_type_id].total_cells;
  const auto active_cell_types = this->_getValidCellTypes(universe);

  // Calculate background universe total cells
  const std::deque<CellType>& all_cts = this->inverse_cell_type;
  const CellTypeIndex& cts_index = this->cell_types;
  const int act_total_cells = std::accumulate(active_cell_types.begin(), 
                                              active_cell_types.end(),
                                              0,
                                              [&all_cts, &cts_index](const int& sum, const CellTypeName& name){
                                                const auto ct_id = cts_index.find(name);
                                                return sum + all_cts[ct_id->second].total_cells;
                                              });
  
  std::map<GeneName, CellTypeMarker> scores;
  if (mode == ALL)
  {
    CellTypeMarker ctm_template = {0, 0, 0, 0};
  
    for (auto const& gene_name : gene_names)
    {  
      const auto index_it = this->index.find(gene_name);
      if (index_it == this->index.end())
      {
        Rcpp::Rcerr << "Gene " << gene_name << " not found in the database, Ignoring... " << std::endl;
        continue;
      }

      const auto& gene_entry = *index_it;
      // Make sure the cell type is in then batch
      auto ctm = gene_entry.second.find(cell_type_id);
      if (ctm == gene_entry.second.end())
      {
        continue;
      }

      auto dit = scores.insert(std::make_pair(gene_entry.first, ctm_template));
      CellTypeMarker& gene_ctm_score = dit.first->second;
      // this->cellTypeMarkerGeneMetrics(gene_ctm_score);
      const EliasFano& ex_vec = this->ef_data[ctm->second];
      // Rcpp::Rcout << ctm->second << " size: " << this->ef_data[ctm->second]. << std::endl;
      int cells_in_ct = ex_vec.getSize();
      gene_ctm_score.tp = cells_in_ct;
      gene_ctm_score.fn = total_cells_in_ct - gene_ctm_score.tp;
      for (auto const& ct : gene_entry.second)
      {
        auto bct_it = active_cell_types.find(all_cts[ct.first].name);
        // if we are not interested in the cell type continue
        if ( bct_it == active_cell_types.end())
        {
          continue;
        }
      
        int bkg_cell_number = this->ef_data[ct.second].getSize();
        gene_ctm_score.fp += bkg_cell_number;
      }

      // total cells in the universe - total cells expressing gene - total cells not in cell type
      gene_ctm_score.tn = act_total_cells - gene_ctm_score.fp - total_cells_in_ct;

      // subtract the cells in the cell_type of interest
      gene_ctm_score.fp -= cells_in_ct;
    }
  }
  else if (mode == AND)
  {
    CellTypeMarker ctm_template = {0, 0, 0, 0};
    std::ostringstream imploded;
    std::copy(gene_names.begin(), gene_names.end(),
           std::ostream_iterator<std::string>(imploded, ","));
    auto dit = scores.insert(std::make_pair(imploded.str(), ctm_template));
    CellTypeMarker& gene_ctm_score = dit.first->second;
    
    
    Rcpp::List result = this->_findCellTypes(gene_names, std::vector<CellTypeName>(active_cell_types.begin(), active_cell_types.end()));
    std::vector<CellTypeName> cell_types = Rcpp::as< std::vector<CellTypeName> >(result.names());
    
    
    for (auto const& ct : cell_types)
    {
      auto cells_positive = Rcpp::as<Rcpp::IntegerVector>(result[ct]).size();
      auto cells_negative = this->inverse_cell_type.at(this->cell_types.at(ct)).total_cells - cells_positive;
      if ( ct == cell_type)
      {
        gene_ctm_score.tp = cells_positive;
        gene_ctm_score.fn += cells_negative;
      }
      else
      {
        gene_ctm_score.fp += cells_positive;
        gene_ctm_score.tn = cells_negative;
      }
    }
    

    
  }


  return scores;
}


const std::vector<EliasFanoDB::CellTypeName> EliasFanoDB::_getCellTypes() const
{
  std::vector<CellTypeName> cts;
  cts.reserve(this->cell_types.size());
  for (auto const& ct : this->inverse_cell_type)
  {
    cts.push_back(ct.name);
  }
  return cts;
}

// Returns a map of cell types and cell id's that contain all genes in the gene_set
std::map<std::string, std::vector<int>> EliasFanoDB::intersect_cells(std::set<std::string> gene_set, Rcpp::List genes_results) const
{
  
  std::map<std::string, std::vector<int> > ct_map;
  
  // First gene init set for intersections
  auto git = gene_set.begin();
     
  const Rcpp::List& first_gene = genes_results[*git];
  const std::vector<std::string> initial_cell_types = Rcpp::as<std::vector<std::string> >(first_gene.names());
  for (auto const& ct : initial_cell_types)
  {
    ct_map[ct] = Rcpp::as< std::vector<int> >(first_gene[ct]);
  }

  for (++git; git != gene_set.end(); ++git)
  {
    const Rcpp::List& g_res = genes_results[*git];
    std::vector<std::string> cell_types = Rcpp::as<std::vector<std::string> >(g_res.names());
    std::vector<std::string> cell_types_in_ct_map;
    cell_types_in_ct_map.reserve(ct_map.size());
    for (auto const& ct : ct_map)
    {
      cell_types_in_ct_map.push_back(ct.first);
    }
        
    std::sort(cell_types_in_ct_map.begin(), cell_types_in_ct_map.end());
    std::sort(cell_types.begin(), cell_types.end());
        
        
    std::vector<std::string> ct_intersection;
    std::set_intersection(cell_types_in_ct_map.begin(),
                          cell_types_in_ct_map.end(),
                          cell_types.begin(),
                          cell_types.end(),
                          std::back_inserter(ct_intersection));

    // Update ct_map with the values that have been thrown out cause of the intersection
    for(auto m_it = ct_map.begin(); m_it != ct_map.end();)
    {
      auto f_it = std::find(ct_intersection.begin(),ct_intersection.end(), m_it->first);
      if (f_it == ct_intersection.end())
      {
        m_it = ct_map.erase(m_it);
      }
      else
      {
        ++m_it;
      }
    }
        
  }

  return ct_map;

}

int EliasFanoDB::dbSize()
{
  Rcpp::Rcout << index.size() << "genes in the DB" << std::endl;
  return ef_data.size();
    
}

std::vector<int> EliasFanoDB::decode(int index)
{
  if(index >= dbSize())
  {
    Rcpp::Rcerr << "Invalid index for database with size "<< dbSize() << std::endl;
    return std::vector<int>();
  }
  return eliasFanoDecoding(ef_data[index]);
}

  
int EliasFanoDB::insertNewCellType(const CellType& cell_type)
{ 
  auto ct_it = this->cell_types.insert(std::make_pair(cell_type.name, this->cell_types.size()));
  
  if (not ct_it.second)
  {

    Rcpp::Rcerr << "This should not happen!! Duplicate Cell Type: " << cell_type.name << std::endl;
  }
  else
  {
    this->inverse_cell_type.push_back(cell_type);
  }
  
  return ct_it.first->second;
}

int EliasFanoDB::mergeDB(const EliasFanoDB& db)
{    
  EliasFanoDB extdb(db);
  if (extdb.getQuantizationBits() != this->getQuantizationBits())
  {
    Rcpp::Rcerr << "Can not perform merging.. Quantization bits are not equal in the two databases. Please fix" << std::endl;
    return 1;
  }
  // the DB will grow by this amount of cells
  this->total_cells += extdb.total_cells;

  // Insert new cell types in the database
  for (auto const& ct : extdb.inverse_cell_type)
  {
    insertNewCellType(ct);
  }
  
  // Iterate through the data model
  for ( auto& gene: extdb.index)
  {
    // Update cell counts for the individual gene
    auto gene_it = this->genes.insert({gene.first, extdb.genes[gene.first]});
    
    if(not gene_it.second)
    {
      gene_it.first->second.merge(extdb.genes[gene.first]);
    }
    
    for (auto& ct : gene.second)
    {
      int new_id = ef_data.size();
      // Push the new elias fano index in the database
      ef_data.push_back(extdb.ef_data[ct.second]);
      // Update with the new entry
      int cell_type_id = this->cell_types[extdb.inverse_cell_type[ct.first].name];
      index[gene.first][cell_type_id] = new_id;
    }
  }

  for (auto const& cell : extdb.cells)
  {

    CellID clone(cell.first);
    int old_cell_type_id = clone.cell_type;
    // Trace the new entry in the database
    int new_cell_type_id = this->cell_types[extdb.inverse_cell_type[old_cell_type_id].name];
    clone.cell_type = new_cell_type_id;
    this->cells.insert({clone, cell.second});
  }


  return 0;
}


// Rcpp will take care of the wrapping
Rcpp::List EliasFanoDB::getCellMeta(const std::string& ct, const int& num) const
{
  const auto ct_it = this->cell_types.find(ct);
  const CellID cid(ct_it->second, num);
  const auto cmeta_it = this->cells.find(cid);
  const CellMeta& cmeta = cmeta_it->second;
  return Rcpp::List::create(
    Rcpp::Named("total_reads") = Rcpp::wrap(cmeta.getReads()),
    Rcpp::Named("total_features") = Rcpp::wrap(cmeta.getFeatures()));
                            
                            
}
    
Rcpp::List EliasFanoDB::getCellTypeMeta(const std::string& ct_name) const 
{
  const auto ct_it = this->cell_types.find(ct_name);
  const CellType& ctmeta = this->inverse_cell_type[ct_it->second];
  return Rcpp::List::create(Rcpp::Named("total_cells") = ctmeta.getTotalCells());
}


RCPP_MODULE(EliasFanoDB)
{
  Rcpp::class_<EliasFanoDB>("EliasFanoDB")
    .constructor()
    .method("setQB", &EliasFanoDB::setQuantizationBits)
    .method("indexMatrix", &EliasFanoDB::encodeMatrix)
    .method("queryGenes", &EliasFanoDB::queryGenes)
    .method("zgs", &EliasFanoDB::queryZeroGeneSupport)
    .method("decode", &EliasFanoDB::decode)
    .method("mergeDB", &EliasFanoDB::mergeDB)
    .method("findCellTypes", &EliasFanoDB::findCellTypes)
    .method("efMemoryFootprint", &EliasFanoDB::dataMemoryFootprint)
    .method("dbMemoryFootprint", &EliasFanoDB::dbMemoryFootprint)
    .method("quantFootprint", &EliasFanoDB::quantizationMemoryFootprint)
    .method("findMarkerGenes", &EliasFanoDB::findMarkerGenes)
    .method("numberOfCellTypes", &EliasFanoDB::numberOfCellTypes)
    .method("getByteStream", &EliasFanoDB::getByteStream)
    .method("loadByteStream", &EliasFanoDB::loadByteStream)
    .method("getTotalCells", &EliasFanoDB::getTotalCells)
    .method("genes", &EliasFanoDB::getGenesInDB)
    .method("genesSupport", &EliasFanoDB::totalCells)
    .method("geneSupportInCellTypes", &EliasFanoDB::geneSupportInCellTypes)
    .method("cellTypeMarkers", &EliasFanoDB::findCellTypeMarkers)
    .method("getCellTypes", &EliasFanoDB::getCellTypes)
    .method("getCellMeta", &EliasFanoDB::getCellMeta)
    .method("getCellTypeExpression", &EliasFanoDB::getCellTypeMatrix)
    .method("getCellTypeMeta", &EliasFanoDB::getCellTypeMeta)
    .method("evaluateCellTypeMarkersAND", &EliasFanoDB::evaluateCellTypeMarkersAND)
    .method("evaluateCellTypeMarkers", &EliasFanoDB::evaluateCellTypeMarkers)
    .method("getCellTypeSupport", &EliasFanoDB::getCellTypeSupport);
}


