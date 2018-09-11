#include <cmath>
#include <iterator>
#include <numeric>
#include <functional>
#include <exception>
#include <stdexcept>

// [[Rcpp::depends(RcppArmadillo)]]

#include "scfind_types.h"
#include "fp_growth.hpp"
#include "functions.h"
#include "EliasFano.hpp"
#include "Serialization.hpp"




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
 
QueryScore::QueryScore() : query_score(0), cells_in_query(0)
{}

void QueryScore::cell_type_relevance(const EliasFanoDB& db, const Rcpp::List& genes_results, const std::set<std::string>& gene_set)
{

  std::map<std::string, std::vector<int> > ct_map = db.intersect_cells(gene_set, genes_results);
  
  for (auto const& _ct : ct_map)
  {
    const auto& ct = _ct.first;
    for (auto& gene : gene_set)
    {
      
      
      auto current_cells = std::move(ct_map[ct]);
      auto current_cells_results = std::vector<int>();
      
      const Rcpp::List& g_res = genes_results[gene];
      std::vector<int> cells_in_gct = Rcpp::as<std::vector<int> >(g_res[ct]);
        
      // We do not need sorting the index is already sorted
      // std::sort(current_cells.begin(), current_cells.end());
      // std::sort(cells_in_gct.begin(), cell_in_gct.end
      std::set_intersection(
        current_cells.begin(),
        current_cells.end(),
        cells_in_gct.begin(),
        cells_in_gct.end(),
        std::back_inserter(current_cells_results));
      
      ct_map[ct] = std::move(current_cells_results);
      if (current_cells_results.empty())
      {
        break;
      }
    }
  }
    
  
    // Remove all empty records
  for (auto it = ct_map.begin(); it != ct_map.end();)
  {
    it = it->second.empty() ? ct_map.erase(it) : ++it;
  }
  for (auto const& _ct : ct_map)
  {
    const std::string& ct = _ct.first;
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

void QueryScore::estimateExpression(const Rcpp::List& gene_results, const EliasFanoDB& db)
{
  std::cout << "calculating tfidf for the reduced expression matrix... " << std::endl;
  // Store temporarily the strings so we can insert those in the map
  // TODO(Nikos) check if genes are unique in the set.. Possibly this can be done in the R side ?
  
  const auto& tmp_strings = Rcpp::as<std::vector<std::string>>(gene_results.names());
  const auto tmpl_cont = std::vector<double>(tmp_strings.size(), 0);
  int gene_row = 0;
  for (size_t gene_row = 0; gene_row < tmp_strings.size(); ++gene_row)
  {

    const std::string& gene = tmp_strings[gene_row];
    // std::cerr << "Gene: " << gene << std::endl;
    
    float gene_idf = db.getTotalCells() / ((float)db.genes.at(gene).total_reads);
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
      if (expr_values.size() != expr_indices.size())
      {
        std::cerr << "Corrupted DB!" << std::endl;
      }
      int expr_index = 0;
      for (auto const& cell_id : expr_indices)
      {
        CellID cell(ct_id, cell_id);
        
        auto ins_res = tfidf.insert(std::make_pair(cell, std::make_pair(tmpl_cont, 0)));
        auto& tfidf_vec = ins_res.first->second.first;
        auto& gene_support = ins_res.first->second.second;
        gene_support++;
        
        tfidf_vec[gene_row] = (expr_values[expr_index++] / db.cells.at(cell).reads) * gene_idf;
        gene_score += tfidf_vec[gene_row];
        
      }
    }
  }
  // iterate through cells
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
    std::cerr << v.first << " " << genes_subset[v.second.index] << std::endl;
  }
  
  std::cout << "Done!" << std::endl;
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
  this->cells_in_query *= min;
  
}



void EliasFanoDB::clearDB()
{
  // Clear the database
  index.clear();
  cell_types.clear();
  inverse_cell_type.clear();
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
  // std::cerr <<"New index" << ef_data.size() - 1 << std::endl;
  // return the index of the ef_data in the deque
  return ef_data.size() - 1;
}


std::vector<int> EliasFanoDB::eliasFanoDecoding(const EliasFano& ef)
{
    
  
  // This step inflates the vector by a factor of 8
  std::vector<char> H;
  std::vector<int> ids(ef.L.size() / ef.l);
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



// TODO(Nikos) returning the first entry of the database is dodgy, fix that
const EliasFano& EliasFanoDB::getEntry(const GeneName& gene_name, const CellTypeName& cell_type) const 
{
  try
  {
    return this->ef_data.at(this->index.at(gene_name).at(this->cell_types.at(cell_type)));                       
  } 
  catch(const std::out_of_range& e)
  {
    std::cerr << e.what() << std::endl;
    auto g_it = index.find(gene_name);
    if (g_it == index.end())
    {
      std::cerr << gene_name << "Gene not found" << std::endl;
    }
    auto ct_it = this->cell_types.find(cell_type);
    
    if (ct_it == this->cell_types.end())
    {
      std::cerr << "Cell type"<<cell_type<<" not found in the database" << std::endl;
    }
    else
    {
      auto ef_it = g_it->second.find(ct_it->second);
      
      if(ef_it == g_it->second.end())
      {
        std::cerr << "Cell type "<< cell_type<<" not found for gene " << gene_name << std::endl;
        
      }
    }
    throw std::invalid_argument("Unable to retrieve entry from database");
  }
}



// constructor
EliasFanoDB::EliasFanoDB(): 
  global_indices(false), 
  warnings(0), 
  total_cells(0), 
  quantization_bits(2)
{
    
}


// This is invoked on slices of the expression matrix of the dataset 
long EliasFanoDB::encodeMatrix(const std::string& cell_type_name, const Rcpp::NumericMatrix& gene_matrix)
{
  int items = gene_matrix.ncol();
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
  for (unsigned int gene_row = 0; gene_row < gene_matrix.nrow(); ++gene_row)
  {
    // Insert the gene gene metadata database
    auto gene_it = this->genes.insert(std::make_pair(gene_names[gene_row], GeneMeta())).first;
    // Insert the gene in the index
    auto db_entry = this->index.insert(std::make_pair(gene_names[gene_row], GeneContainer())).first;

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
      std::cerr << "Vector of zeros detected for cell " << cell_type_name << " " << i << std::endl;
    }
    this->cells.insert({CellID(cell_type_id, ++i), cell});
  }
  
  return 0;
  //std::cerr << "Total Warnings: "<<warnings << std::endl;
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

Rcpp::IntegerVector EliasFanoDB::totalCells(const Rcpp::CharacterVector& genes)
{
  Rcpp::IntegerVector t(genes.size(), 0);
  t.names() = genes;
  std::vector<std::string> str = Rcpp::as<std::vector<std::string>> (genes);
  int i = 0;
  for (auto const& g : str)
  {
    auto git = this->index.find(g);
    if(git != this->index.end())
    {
      for (auto const& ct : git->second)
      {
        t[i] += this->ef_data[ct.second].getSize();
      }
    }
  
    i++;
  }
  return t;
  
}

int EliasFanoDB::getTotalCells() const
{
  return this->total_cells;
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
  
std::vector<double> EliasFanoDB::getQuantizedExpressionLevels(const std::string& gene_name, const std::string& cell_type)
{
  const EliasFano& ef = this->ef_data[this->index[gene_name][this->cell_types[cell_type]]];
  return decompressValues(ef.expr, this->quantization_bits);
  
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
        
      std::cout << "Gene " << gene_name << " not found in the index " << std::endl;
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

  std::cout << "Raw elias Fano Index size " << bytes / (1024 * 1024) << "MB" << std::endl;
    

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


// And query
Rcpp::List EliasFanoDB::findCellTypes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& datasets_active)
{
    
  std::unordered_map<CellTypeID, std::set<std::string> > cell_types;
  std::vector<std::string> genes;

  std::vector<std::string> datasets = Rcpp::as<std::vector<std::string>>(datasets_active);

  // Fast pruning if there is not an entry we do not need to consider
  for (Rcpp::CharacterVector::const_iterator it = gene_names.begin(); it != gene_names.end(); ++it)
  {
    std::string gene_name = Rcpp::as<std::string>(*it);
      
    // check if gene exists in the database
    auto db_it = index.find(gene_name);
    if (db_it == index.end())
    {
      std::cerr << gene_name << " is ignored, not found in the index"<< std::endl;
      continue;
    }

    // iterate cell type
    for (auto const& ct_it : db_it->second)
    {
      CellType ct = this->inverse_cell_type[ct_it.first];
        
      // Remove cells if not in the selected datasets
      std::string ct_dataset = ct.name.substr(0, ct.name.find("."));
      auto find_dataset = std::find(datasets.begin(), datasets.end(), ct_dataset);
      // check if the cells are in active datasets
      if (find_dataset == datasets.end())
      {
        continue;
      }

      if (cell_types.find(ct_it.first) == cell_types.end())
      {
        cell_types[ct_it.first] = std::set<std::string>();
      }
      cell_types[ct_it.first].insert(gene_name);
    }
    genes.push_back(gene_name);
  }
    
  // Store the results here
  Rcpp::List t;

  for (auto const& ct : cell_types)
  {
    bool empty_set = false;
    bool initial_set = true;
    if (ct.second.size() != genes.size())
    {
      continue;
    }
      
    auto g_it = ct.second.begin();
    std::vector<int> ef = eliasFanoDecoding(this->ef_data[index[*g_it][ct.first]]);

    for (++g_it; g_it != ct.second.end();++g_it)
    {
      auto cells = eliasFanoDecoding(this->ef_data[index[*g_it][ct.first]]);
        
      std::vector<int> intersected_cells;
      std::set_intersection(ef.begin(), 
                            ef.end(), 
                            cells.begin(), 
                            cells.end(), 
                            std::back_inserter(intersected_cells));
      ef = intersected_cells;
      if (ef.empty())
      {
        break;
      }
    }
      
    if(!ef.empty())
    {
      t[this->inverse_cell_type[ct.first].name] = Rcpp::wrap(ef);
    }
  }
  return t;
}





// TODO(Nikos) this function can be optimized.. It uses the native quering mechanism
// that casts the results into native R data structures
Rcpp::DataFrame EliasFanoDB::findMarkerGenes(const Rcpp::CharacterVector& gene_list, const Rcpp::CharacterVector datasets_active,unsigned int min_support_cutoff = 5)
{
    
  std::vector<std::string> query;
  std::vector<double> query_scores;
  std::vector<double> query_tfidf;
  std::vector<int> query_cell_type_cardinality;
  std::vector<int> query_cell_cardinality;
  std::vector<int> query_gene_cardinality;

  std::map<CellTypeName, std::map<int, Transaction> > cells;

  int cells_present = 0;
    
  // Perform an OR query on the database as a first step
  const Rcpp::List genes_results = queryGenes(gene_list, datasets_active);

  // Start inversing the list to a cell level
  const std::vector<std::string> gene_names = Rcpp::as<std::vector<std::string>>(genes_results.names());
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
          cells_present++;
        }
        // Add gene hit in the cell
        cell_index[id].push_back(gene_name);
      }
    }
  }

  std::cerr << "Query Done: found " << cells_present << " rules" << std::endl;
    
  // Collect all transactions for fp-growth
  std::vector<Transaction> transactions;
  transactions.reserve(cells_present);
  for (auto & ct : cells)
  {
    for (auto & cl : ct.second)
    {
      // Maybe sort?        
      // std::sort(cl.second.begin(), cl.second.end());
      if (cl.second.size() != 1)
      {
        transactions.push_back(std::vector<std::string>(cl.second.begin(), cl.second.end()));
      }
    }
  }
    
  std::cerr << transactions.size() << " transactions" << std::endl;
  // Run fp-growth algorithm

  QueryScore qs;
  qs.estimateExpression(genes_results, *this);
  
  std::cerr << "Running fp-growth tree with " << min_support_cutoff << " cutoff"<< std::endl;
  const FPTree fptree{transactions, min_support_cutoff};
  const std::set<Pattern> patterns = fptree_growth(fptree);
  std::cerr << "Found " << patterns.size() << " geneset patterns " << std::endl;
  
  // Iterate through the calculated frequent patterns
  
  
  for (auto const& item : patterns)
  {
    
    
    Rcpp::List gene_query;
    const auto& gene_set = item.first;
    int fp_support = item.second;
    
    // We do not care for queries with cardinality less than 2
    if (gene_set.size() < 2)
    {
      continue;
    }
    
    std::string view_string = str_join(std::vector<Item>(gene_set.begin(), gene_set.end()), ",");

    
    // cell_type_relevance
    qs.reset();
    qs.cell_type_relevance(*this, genes_results, gene_set);
    query_scores.push_back(qs.query_score);
    query_cell_cardinality.push_back(qs.cells_in_query);
    query_cell_type_cardinality.push_back(qs.cell_types_in_query);
    
    // query tfidf
    qs.reset();
    qs.cell_tfidf(*this, gene_set);
    query_tfidf.push_back(qs.query_score);
    
    // other fields
    query_gene_cardinality.push_back(gene_set.size());
    query.push_back(view_string);
  }

  std::vector<int> query_rank(query_scores.size());    
  std::iota(query_rank.begin(), 
            query_rank.end(), 
            1);
  
  std::sort(query_rank.begin(), 
            query_rank.end(), 
            [&query_scores](const int& i1, const int& i2){
              return query_scores[i1] < query_scores[i2];
            });

  // Dump the list
  return Rcpp::DataFrame::create(Rcpp::Named("Genes") = Rcpp::wrap(query_gene_cardinality),
                                 Rcpp::Named("Query") = Rcpp::wrap(query),
                                 Rcpp::Named("Rank") = Rcpp::wrap(query_rank),
                                 Rcpp::Named("tfidf") = Rcpp::wrap(query_tfidf),
                                 Rcpp::Named("Cells") = Rcpp::wrap(query_cell_cardinality),
                                 Rcpp::Named("Cell Types") = Rcpp::wrap(query_cell_type_cardinality));
}





const std::set<std::string> EliasFanoDB::_getActiveCellTypes(std::vector<std::string> universe) const
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
      std::cerr << "Ignoring cell type "<< ct <<". Not found in DB" << std::endl;
    }
  }

  return active_cell_types;

}


Rcpp::DataFrame EliasFanoDB::findCellTypeMarkers(const Rcpp::CharacterVector& cell_types, const Rcpp::CharacterVector& background)
{
  std::vector<GeneName> gene_set;
  gene_set.reserve(this->genes.size());
  for (auto const& g : this->genes)
  {
    gene_set.push_back(g.first);
  }
  return _findCellTypeMarkers(cell_types, background, gene_set);
}


Rcpp::DataFrame EliasFanoDB::_findCellTypeMarkers(const Rcpp::CharacterVector& cell_types, const Rcpp::CharacterVector& background, const std::vector<EliasFanoDB::GeneName>& gene_set)
{
  std::vector<std::string> 
    bk_cts(Rcpp::as< std::vector<std::string> >(background)),
    cts(Rcpp::as< std::vector<std::string> >(cell_types)),
    genes, df_cell_type;
  
  
  std::vector<int> tp, fp, tn, fn;
  std::vector<float> precision, recall, f1;
  for (auto const& ct : cts)
  {
    auto marker_genes = this->_cellTypeScore(ct, bk_cts, gene_set);
    if (marker_genes.empty())
    {
      std::cerr << "Marker genes could not be found for cell type " << ct << std::endl;
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
    Rcpp::Named("tn") = Rcpp::wrap(tn),
    Rcpp::Named("fp") = Rcpp::wrap(fp),
    Rcpp::Named("fn") = Rcpp::wrap(fn),
    Rcpp::Named("precision") = Rcpp::wrap(precision),
    Rcpp::Named("recall") = Rcpp::wrap(recall),
    Rcpp::Named("f1") = Rcpp::wrap(f1)
                         );
  
}


Rcpp::DataFrame EliasFanoDB::evaluateCellTypeMarkers(const Rcpp::CharacterVector& cell_types, 
                                                     const Rcpp::CharacterVector& gene_set,  
                                                     const Rcpp::CharacterVector& background)
{
  return _findCellTypeMarkers(cell_types, background, Rcpp::as<std::vector<GeneName>>(gene_set));
}



std::map<EliasFanoDB::GeneName, CellTypeMarker> EliasFanoDB::_cellTypeScore(const std::string& cell_type, const std::vector<std::string>& universe, const std::vector<EliasFanoDB::GeneName>& gene_names) const
{
  auto ct_it = this->cell_types.find(cell_type);
  if ( ct_it == this->cell_types.end())
  {
    std::cerr <<"Cell type "<< cell_type << " not found. exiting..." << std::endl;
    return std::map<EliasFanoDB::GeneName, CellTypeMarker>();
  }

  const CellTypeID cell_type_id = ct_it->second;
  int total_cells_in_ct = this->inverse_cell_type[cell_type_id].total_cells;
  const auto active_cell_types = this->_getActiveCellTypes(universe);

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
  CellTypeMarker ctm_template = {0, 0, 0, 0};
  
  for (auto const& gene_name : gene_names)
  {  
    const auto index_it = this->index.find(gene_name);
    if (index_it == this->index.end())
    {
      std::cerr << "Gene " << gene_name << " not found in the database, Ignoring... " << std::endl;
      continue;
    }

    const auto& gene_entry = *index_it;
    // Make sure the cell type is in the batch
    auto ctm = gene_entry.second.find(cell_type_id);
    if (ctm == gene_entry.second.end())
    {
      continue;
    }

    auto dit = scores.insert(std::make_pair(gene_entry.first, ctm_template));
    CellTypeMarker& gene_ctm_score = dit.first->second;
    // this->cellTypeMarkerGeneMetrics(gene_ctm_score);
    const EliasFano& ex_vec = this->ef_data[ctm->second];
    // std::cout << ctm->second << " size: " << this->ef_data[ctm->second]. << std::endl;
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
      // if this is the current cell type we are assesing continue as well
      // if ( bct_it->first == cell_type_id) // not needed if we subtract the cells_in_ct later
      // {
      //   continue;
      // }
      
      int bkg_cell_number = this->ef_data[ct.second].getSize();
      gene_ctm_score.fp += bkg_cell_number;
    }

    // total cells in the universe - total cells expressing gene - total cells not in cell type
    gene_ctm_score.tn = act_total_cells - gene_ctm_score.fp - total_cells_in_ct;
    // subtract the cells in the cell_type of interest
    gene_ctm_score.fp -= cells_in_ct;



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
  std::cout << index.size() << "genes in the DB" << std::endl;
  return ef_data.size();
    
}

std::vector<int> EliasFanoDB::decode(int index)
{
  if(index >= dbSize())
  {
    std::cerr << "Invalid index for database with size "<< dbSize() << std::endl;
    return std::vector<int>();
  }
  return eliasFanoDecoding(ef_data[index]);
}
  
int EliasFanoDB::insertNewCellType(const CellType& cell_type)
{ 
  auto ct_it = this->cell_types.insert(std::make_pair(cell_type.name, this->cell_types.size()));
  
  if (not ct_it.second)
  {

    std::cerr << "This should not happen!! Duplicate Cell Type: " << cell_type.name << std::endl;
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
    auto entry_db = this->index.insert({gene.first, GeneContainer()});
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



RCPP_MODULE(EliasFanoDB)
{
  Rcpp::class_<EliasFanoDB>("EliasFanoDB")
    .constructor()
    .method("indexMatrix", &EliasFanoDB::encodeMatrix)
    .method("queryGenes", &EliasFanoDB::queryGenes)
    .method("decode", &EliasFanoDB::decode)
    .method("mergeDB", &EliasFanoDB::mergeDB)  
    .method("findCellTypes", &EliasFanoDB::findCellTypes)
    .method("efMemoryFootprint", &EliasFanoDB::dataMemoryFootprint)
    .method("dbMemoryFootprint", &EliasFanoDB::dbMemoryFootprint)
    .method("findMarkerGenes", &EliasFanoDB::findMarkerGenes)
    .method("getByteStream", &EliasFanoDB::getByteStream)
    .method("loadByteStream", &EliasFanoDB::loadByteStream)
    .method("getCellsInDB", &EliasFanoDB::getTotalCells)
    .method("genes", &EliasFanoDB::getGenesInDB)
    .method("genesSupport", &EliasFanoDB::totalCells)
    .method("cellTypeMarkers", &EliasFanoDB::findCellTypeMarkers)
    .method("getCellTypes", &EliasFanoDB::_getCellTypes)
    .method("evaluateCellTypeMarkers", &EliasFanoDB::evaluateCellTypeMarkers)
    .method("getCellTypeSupport", &EliasFanoDB::getCellTypeSupport);
  
}


