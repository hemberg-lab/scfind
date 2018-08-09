
#include <cmath>
#include <iterator>
#include <numeric>
#include <functional>
// [[Rcpp::depends(RcppArmadillo)]]

#include "scfind_types.h"
#include "fp_growth.hpp"
#include "functions.h"
#include "EliasFano.h"
#include "Serialization.hpp"

void EliasFanoDB::dumpGenes()
{
    
  std::cout << "Total Genes:" << gene_counts.size() << std::endl;
  for (auto const c : this->cell_types_id)
  {
    std::cout << c.first.name << " ";
  }
  std::cout << "Total Cell types" << std::endl;
    
  for (auto const& m : metadata)
  {
    std::cout << "Gene:" << m.first << std::endl;
    long total = 0;
    for (auto const& t : m.second)
    {
      auto ct = inverse_cell_type[t.first];
      // str.size();
      total += ct.name.size() + 4;
    }
    std::cout << total <<std:: endl;
  }
} 

void EliasFanoDB::clearDB()
{
  // Clear the database
  metadata.clear();
  cell_types_id.clear();
  inverse_cell_type.clear();
  gene_counts.clear();
  // serialized_bytestream.clear();
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


void EliasFanoDB::insertToDB(int ef_index, const std::string& gene_name, const CellType& cell_type)
{
  if (ef_index == -1)
  {
    // Something went wrong so do not do anything
    this->warnings++;
    return;
  }

  EliasFano* ef = &(this->ef_data[ef_index]);
    
  if (metadata.find(gene_name) == metadata.end())
  {
    metadata[gene_name] = EliasFanoIndex();
  }
    
  auto celltype_record = this->cell_types_id.find(cell_type);
  if (celltype_record == this->cell_types_id.end())
  {
    int id = this->cell_types_id.size();
    this->cell_types_id[cell_type] = id;
    this->inverse_cell_type.push_back(cell_type);
    metadata[gene_name][id] = ef_index;
  }else
  {
    metadata[gene_name][this->cell_types_id[cell_type]] = ef_index;
  }
    
    
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
  Rcpp::CharacterVector genes = Rcpp::rownames(gene_matrix);

  CellType cell_type;
  cell_type.name = cell_type_name;
  cell_type.total_cells = gene_matrix.ncol();

  std::vector<std::string> gene_names;
  gene_names.reserve(genes.size());


  // Increase the cell number present in the index
  this->total_cells += gene_matrix.ncol();
    
  for(Rcpp::CharacterVector::iterator gene_it = genes.begin(); gene_it != genes.end(); ++gene_it)
  {
    std::string gene_name = Rcpp::as<std::string>(*gene_it);
    gene_names.push_back(gene_name);
    // If this is the first time this gene occurs initialize a counter
    if (this->gene_counts.find(gene_name) == this->gene_counts.end())
    {
      this->gene_counts[gene_name] = 0;
    }
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
    this->gene_counts[gene_names[gene_row]] += ids.size();
      
    insertToDB(eliasFanoCoding(ids, expression_vector), gene_names[gene_row], cell_type);
  }
  return 0;
  //std::cerr << "Total Warnings: "<<warnings << std::endl;
}

Rcpp::List EliasFanoDB::total_genes()
{
  Rcpp::List t;
  for(auto & d : metadata)
  {
    t.push_back(Rcpp::wrap(d.first));
  }
  return t;
}

Rcpp::NumericVector EliasFanoDB::getTotalCells()
{
  return Rcpp::wrap(this->total_cells);
}

Rcpp::NumericVector EliasFanoDB::getCellTypeSupport(Rcpp::CharacterVector& cell_types)
{
  std::vector<std::string> cts = Rcpp::as<std::vector<std::string>>(cell_types);
  std::vector<int> ct_support;
  ct_support.reserve(cts.size());
  for (auto const& ct : cts)
  {
    // TODO(Nikos) fix this, otherwise we will get a nice segfault no error control
    auto cit = this->cell_types_id.find(ct);
    if(cit != this->cell_types_id.end())
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
      
    if (metadata.find(gene_name) == metadata.end())
    {
        
      std::cout << "Gene " << gene_name << " not found in the index " << std::endl;
      continue;
    }
    std::vector<std::string> datasets = Rcpp::as<std::vector<std::string>>(datasets_active);
    const auto& gene_meta = metadata[gene_name];
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
    bytes += int((d.expr.quantile.size() / 8) + 1);
      
  }
  bytes += ef_data.size() * 32; // overhead of l idf and deque struct
  return bytes;
}

size_t EliasFanoDB::dbMemoryFootprint()
{
  size_t bytes = dataMemoryFootprint();

  std::cout << "Raw elias Fano Index size " << bytes / (1024 * 1024) << "MB" << std::endl;
    
  for(auto& d : metadata)
  {
    bytes += d.first.size();
    bytes += d.second.size() * 8;
  }
  return bytes;
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
    auto db_it = metadata.find(gene_name);
    if (db_it == metadata.end())
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
    std::vector<int> ef = eliasFanoDecoding(this->ef_data[metadata[*g_it][ct.first]]);

    for (++g_it; g_it != ct.second.end();++g_it)
    {
      auto cells = eliasFanoDecoding(this->ef_data[metadata[*g_it][ct.first]]);
        
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
  // Store the results in this list
  Rcpp::List results;
    
  std::vector<std::string> query;
  std::vector<double> query_scores;
  std::vector<int> query_cell_type_cardinality;
  std::vector<int> query_cell_cardinality;
  std::vector<int> query_gene_cardinality;

  std::map<std::string, std::map<int, Transaction> > cells;

  int cells_present = 0;
    
  // Perform an OR query on the database as a first step
  const Rcpp::List genes_results = queryGenes(gene_list, datasets_active);

  // Start inversing the list to a cell level
  const Rcpp::CharacterVector& gene_names = genes_results.names();
  for (auto const& gene : gene_names)
  {
    const auto gene_name = Rcpp::as<std::string>(gene);
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
  const FPTree fptree{transactions, min_support_cutoff};
  const std::set<Pattern> patterns = fptree_growth(fptree);

  // Patterns contain genesets
  std::vector<std::pair<std::string, double> > tfidf;
    
  // Iterate through the calculated frequent patterns
  for (auto const& item : patterns)
  {
    Rcpp::List gene_query;
    const auto& gene_set = item.first;
    int fp_support = item.second;
      
    if (gene_set.size() == 1)
    {
      continue;
    }
    
    // TODO(Nikos) optimize wrapping 
    std::map<std::string, std::vector<int> > ct_map = intersect_cells(gene_set, genes_results);
      
    if (ct_map.empty())
    {
      continue;
    }
    // std::cout << " Common cell types accross genes" << ct_map.size() << std::endl;
    for (auto const& _ct : ct_map)
    {
      auto ct = _ct.first;
      for(auto& gene : gene_set){
        auto& current_cells = ct_map[ct];
        const Rcpp::List& g_res = genes_results[gene];
        std::vector<int> cells;
        std::vector<int> cells_in_gct = Rcpp::as<std::vector<int>>(g_res[ct]);
        
        // We do not need sorting the index is already sorted
        // std::sort(current_cells.begin(), current_cells.end());
        // std::sort(cells_in_gct.begin(), cell_in_gct.end
        std::set_intersection(
          current_cells.begin(),
          current_cells.end(),
          cells_in_gct.begin(),
          cells_in_gct.end(),
          std::back_inserter(cells));
          
        // This invalidates the current cells reference, careful there
        ct_map[ct] = cells;
        
        if (cells.empty())
        {
          ct_map.erase(ct_map.find(ct));
          break;
        }
      }
    }
    
    
    // Remove all empty records
    for (auto it = ct_map.begin(); it != ct_map.end();)
    {
      it = it->second.empty() ? ct_map.erase(it) : ++it;
    }
    
    // Continue to the next pattern
    if (ct_map.empty())
    {
      continue;
    }

    // std::cout << "FP-Growth support " << fp_support << ", Discovered support: " << total_support << std::endl;
     
    // Give query a score
    std::string view_string = str_join(std::vector<Item>(gene_set.begin(), gene_set.end()), ",");
    double query_score = 0;
    int cells_in_query = 0;
    std::map<std::string, double> ct_tfidf;
    for ( auto const& _ct : ct_map)
    {
      const std::string& ct = _ct.first;
      auto cct =  this->cell_types_id.find(_ct.first);
      double tfidf = 0;
      for (auto const&  g : gene_set)
      {

        tfidf += this->ef_data[this->metadata[g][cct->second]].idf;
      }
      query_score += tfidf * log(_ct.second.size());
      cells_in_query += _ct.second.size();
    }
    // Get the mean ct size
    query_score /= ct_map.size();
      
      
    query_cell_cardinality.push_back(cells_in_query);
    query_gene_cardinality.push_back(gene_set.size());
    query.push_back(view_string);
    query_scores.push_back(query_score);
    query_cell_type_cardinality.push_back(ct_map.size());
      
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
                                 Rcpp::Named("Cells") = Rcpp::wrap(query_cell_cardinality),
                                 Rcpp::Named("Cell Types") = Rcpp::wrap(query_cell_type_cardinality));
}


// Returns a map of cell types and cell id's that contain all genes in the gene_set
std::map<std::string, std::vector<int>> EliasFanoDB::intersect_cells(std::set<std::string> gene_set, Rcpp::List genes_results)
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
      
  // std::cout << *git << ", cell types " << ct_map.size() << std::endl;

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
  std::cout << metadata.size() << "genes in the DB" << std::endl;
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
  // Fetch the last id
  int id = this->inverse_cell_type.size();

  if ( this->cell_types_id.find(cell_type) != this->cell_types_id.end())
  {

    std::cerr << "This should not happen!! Duplicate Cell Type: " << cell_type.name << std::endl;
    id = this->cell_types_id[cell_type];
  }
  else
  {
    id = this->inverse_cell_type.size();
    this->inverse_cell_type.push_back(cell_type);
    this->cell_types_id[cell_type] = id;
  }
  return id;
    
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
  for ( auto& gene: extdb.metadata)
  {
    // Update cell counts for the individual gene
    if (this->gene_counts.find(gene.first) == this->gene_counts.end())
    {
      this->gene_counts[gene.first] = 0;
    }
      
    this->gene_counts[gene.first] += extdb.gene_counts[gene.first];
      
      
    // if gene does not exist yet initialize entry in metadata
    if(metadata.find(gene.first) == metadata.end())
    {
      metadata[gene.first] = EliasFanoIndex();
    }
      
    for( auto& ct : gene.second)
    {
      int new_id = ef_data.size();
      // Push the new elias fano index in the database
      ef_data.push_back(extdb.ef_data[ct.second]);
      // Update with the new entry
      int cell_type_id = this->cell_types_id[extdb.inverse_cell_type[ct.first]];
      metadata[gene.first][cell_type_id] = new_id;
    }
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
    .method("dumpGenes", &EliasFanoDB::dumpGenes)
    .method("getByteStream", &EliasFanoDB::getByteStream)
    .method("loadByteStream", &EliasFanoDB::loadByteStream)
    .method("getCellsInDB", &EliasFanoDB::getTotalCells)
    .method("getCellTypeSupport", &EliasFanoDB::getCellTypeSupport);
  
}


