#' merge two elias fano indices
#'
#' @param efdb.root the root index
#' @param efdb the index to be merged
#'
#' @importFrom hash values keys
#' 
mergeIndices <- function(efdb.root, efdb)
{
    genes.to.merge <- intersect(keys(efdb.root), keys(efdb))
    new.genes <-  setdiff(keys(efdb), keys(efdb.root))
    
    ## print(paste(length(new.genes),length(genes.to.merge)))
    for(gene in new.genes)
    {
        efdb.root[[gene]] <- efdb[[gene]]
        
    }
    for ( gene in genes.to.merge)
    {
        ## maybe we want to merge??
        efdb.root[[gene]][keys(efdb[[gene]])] <- values(efdb[[gene]], simplify = F)
    }
    return(efdb.root)
    
}


contigency.table <- function(query.results)
{
    data <- as.data.frame(lapply(values(query.results), length))
    names(data) <-  keys(query.results)
    return(data)
    
}

caseCorrect <- function(object, gene.list)
{
    gene.list <- gene.list[gene.list != ""]

    if(length(gene.list) != 0)
    {
        gene.corr <- object@index$genes()[match(tolower(gene.list), tolower(object@index$genes()), nomatch = 0)]
        
        if(length(setdiff(tolower(gene.list), tolower(gene.corr))) != 0) message(paste0("Ignored ", toString(setdiff(gene.list, gene.list[match(tolower(gene.corr), tolower(gene.list), nomatch=0)])), ". Not found in the index"))

        return(unique(gene.corr))
    }
    else
    {
        return(c())
    }
}



#' @importFrom stats aggregate p.adjust phyper setNames
phyper.test <- function(object, result, datasets)
{
    df <- query.result.as.dataframe(result)
    datasets <- select.datasets(object, datasets)
    
    cell.types.df <- aggregate(cell_id ~ cell_type, df, FUN = length)
    colnames(cell.types.df)[colnames(cell.types.df) == 'cell_id'] <- 'cell_hits'
    cell.types.df$total_cells<- object@index$getCellTypeSupport(cell.types.df$cell_type)
    query.hits <- nrow(df)
    
    cell.types.df$pval <- p.adjust(1 - phyper(cell.types.df$cell_hits, # total observed successes ( query.hits for cell type)
                                 cell.types.df$total_cells, # total successes ( cell type size )
                                 sum(cell.types.df$total_cells) - cell.types.df$total_cells, # total failures( total cells excluding cell type)
                                 query.hits # sample size 
                                 ), n = length(cellTypeNames(object, datasets)))


    return(cell.types.df)

}

query.result.as.dataframe <- function(query.result)
{
    if (is.data.frame(query.result))
    {
        return(query.result)
    }
    if (length(query.result) == 0)
    {
        return(data.frame(cell_type = c() , cell_id = c()))
    }
    else
    {
        result <- setNames(unlist(query.result, use.names=F), rep(names(query.result), lengths(query.result)))
        return(data.frame(cell_type = names(result), cell_id = result))

    }            
    
}

select.datasets <- function(object, datasets)
{
    
    if (missing(datasets))
    {
        ## Select all available datasets
        datasets <- object@datasets
    }
    else
    {
        ## datasets should not be a superset of the data
        if(length(setdiff(datasets, object@datasets)) != 0)
        {
            stop(paste("Dataset", setdiff(datasets,object@datasets), "does not exist in the database"))
        }
    }
    return(datasets)

}



scfind.get.genes.in.db <- function(object){
    
    return(object@index$genes())

}

pair.id <- function(cell.list = list()){
    if(length(cell.list) == 0) 
    {
        return(c())
    } 
    else
    {
        pair.vec <- stack(cell.list)
        return (paste0(pair.vec$ind, "#",pair.vec$values))
    }
    
}

find.signature <- function(object, cell.type, max.genes=1000, min.cells=10, max.pval=0) {
    # Use this method to find a gene signature for a cell-type. 
    # We do this by ranking genes by recall and then adding genes to the query until we exceed a target p-value threshold or until a minimum number of cells is returned from the query
    df <- cellTypeMarkers(object, cell.type, top.k=max.genes, sort.field="recall", message=F)
    genes <- as.character(df$genes)
    genes.list <- c()
    thres = max(c(min.cells, object@index$getCellTypeMeta(cell.type)$total_cells))
    for (j in 1:dim(df)[1]) {
        res <- hyperQueryCellTypes(object, c(genes.list, genes[j]))
        if (dim(res)[1]>0) {
            ind <- which(res[,1]==cell.type)
            if (length(ind)==0) {
                break
            }
            else {
                if (res[ind,4]>max.pval | res[ind,2]<thres) {
                    break
                }
            }
        }
        genes.list <- c(genes.list, genes[j])
    }
    return( genes.list )
}

## Free Text Search


tokenize <- function(object, query, strict = F, any_id = c("MESH", "CHEBI", "OMIM"))
{
    # Allow user to choose non-operator query treat as AND / OR search, since AND search can hardly find cell type
    mode <- if(strict == T) c('AND#', 'AND_gene', 'AND_word', 'AND_variant') else c('OR#', 'OR_gene', 'OR_word', 'OR_variant')
    
    queries <- unlist(strsplit(query, ",|\\;|\\. "))
    
    queries <- queries[queries != ""]
    queries <- paste0("#", gsub("\\s", "#", trimws(queries, which = "both")))
    
    # Allowing derivatives of the operators
    queries <- gsub('#yes#|#but#also#|#also#|#along#with#|#as#well#as#|#including#|#includes#|#include#|#plus#|#with#|#keep#|#keeping#|add#|#adding#|#both#|#contribute#to#|#connect#with#', '#AND#', queries, ignore.case = T) 
    queries <- gsub('#either#|#link#to#|#associate#with#|#relate#to#', '#OR#', queries, ignore.case = T) 
    queries <- gsub('#no#|#not#either#|#neither#|#nor#|#exclude#|#excluding#|#reject#|#rejecting#|#omit#|#omitting#|#lack#|#lacking#|#minus#|#without#|#ignore#|#ignoring#|#leave#out#|#cancel#|#discard#|#drop#|#take#out#|#except#|#excepting#|#but#not#', '#NOT#', queries, ignore.case = T)
    queries <- gsub('#or#not#', '#ORNOT#', queries, ignore.case = T)
    
    queries <- strsplit(paste(gsub("#", " ", queries), collapse = " # "), " ")[[1]]
    
    all_ids <- paste0(any_id, ":") ## Tweak here if more ids are included
    
    tokens <- list()
    
    # if only gene query without operators
    if(length(suppressMessages(caseCorrect(object, queries[-grep("#", queries)]))) == length(queries) ) 
    { 
        tokens[[mode[2]]] <- suppressMessages(caseCorrect(object, queries))
        return(tokens)
    }
    
    
    conds <- c("AND#", "OR#", "NOT#", "ORNOT#")
    ops <- conds[match(toupper(gsub("#.*", "", queries)), gsub("#", "", conds), nomatch=0)]
    ops.inds <- which(toupper(gsub("#.*", "", queries)) %in%  gsub("#", "", ops)) 
    
    
    if(length(ops) == 0)
    {
        # Process complex query without operators
        genes <- suppressMessages(caseCorrect(object, queries))
        snps <- grep("^rs+\\d+$", sub('*.\\s', '', queries), value = T, ignore.case = T)
        # words <- gsub("#", "\\ ", setdiff(tolower(queries), tolower(c(genes, snps))))
        words <- gsub("\"|'|`", "", setdiff(tolower(queries), tolower(c(genes, snps))))
        tokens[[mode[2]]] <- if(length(genes) != 0) genes else NULL
        tokens[[mode[4]]] <- if(length(snps) != 0) snps else NULL
        tokens[[mode[3]]] <- if(length(words) != 0) words else NULL
        # return(tokens)
        
    }
    else
    {	
        # Process complex query with operators
        
        for(i in length(ops) : if(any(match(ops[1], c(mode[1]), nomatch = F)))1 else 0)
        {
            tk <- NULL
            genes <- NULL
            snps <- NULL
            words <- NULL
            operator  <- if(i > 0) ops[i] else mode[1]
            
            if( i == 1 && any(match(ops[1], c(mode[1]), nomatch = F)))
            {
                # take both side of operators when AND/OR are used 
                tk <- c(queries, "#")
                genes <- suppressMessages(caseCorrect(object, sub(ops[1], "", tk[!tk %in% c(letters, LETTERS)], ignore.case = T)) )
                snps <- grep("^rs+\\d+$", setdiff(tolower(sub(ops[1], "", tk, ignore.case = T)), tolower(genes)), value = T, ignore.case = T)
                
            }
            else
            {
                if (length(queries) != 0) 
                {
                    tk <- if(i == 0) c(queries, "#") else c(queries[ops.inds[i]: length(queries)], "#")
                }
                
                genes <- if(length(tk) != 0) suppressMessages(caseCorrect(object, 
                                                                          sub(operator, "", tk[!tk %in% c(letters, LETTERS)], ignore.case = T))) else NULL
                snps <- if(length(tk) != 0) grep("^rs+\\d+$", 
                                                 setdiff(tolower(sub(operator, "", tk, ignore.case = T)), tolower(genes)), value = T, ignore.case = T) else NULL
                
            }
            
            words <- if(length(tk) != 0) gsub("\"|'|`", "", tk[!tolower(tk) %in% tolower(c(sub("#", "",operator),genes, snps, ""))] ) else NULL
            # queries <- setdiff(queries, tk) 
            queries <- if(i != 0) queries[-(ops.inds[i]:length(queries))] # PACMAN
            
            if (i == 0 && length(queries) != 0) 
            {
                # treat all non-operator query according to mode operator (i.e. AND or OR)
                tokens[[mode[2]]] <- c(tokens[[mode[2]]], if(length(genes) != 0) genes else NULL)
                tokens[[mode[4]]] <- c(tokens[[mode[4]]], if(length(snps) != 0) snps else NULL)
                tokens[[mode[3]]] <- c(tokens[[mode[3]]], if(length(words) != 0) words else NULL)
            }
            else
            {
                if(!is.null(genes) || !is.null(snps) || !is.null(words))
                {
                    tokens[[sub("#", "_gene", operator)]] <- c(tokens[[sub("#", "_gene", operator)]], if(length(genes) != 0) genes else NULL)
                    tokens[[sub("#", "_variant", operator)]] <- c(tokens[[sub("#", "_variant", operator)]], if(length(snps) != 0) snps else NULL)
                    tokens[[sub("#", "_word", operator)]] <- c(tokens[[sub("#", "_word", operator)]], if(length(words) != 0) words else NULL)
                }
            }
            
        }
    }
    
    if(any(grepl(paste(all_ids, collapse = "|"), tokens, ignore.case = T)))
    {
        tokens <- as.matrix(stack(tokens))
        ids_inds <- grep(paste(all_ids, collapse= "|"), tokens[,1], ignore.case=T)
        tokens[ids_inds,2] <- sub("_word", "_mesh", tokens[ids_inds,2])
        tokens <- unstack(data.frame(tokens))
    }
    
    if(class(tokens) != "list")
    {
        return(setNames(strsplit(as.character(tokens$res), ",", fixed = T), rownames(tokens)))
    }
    
    return(tokens)
}

cos.sim.Word2Vec <- function(model, tokens, candidates)
{
    similarities <- list()
    dists <- NULL
    
    julia_assign("m", model)
    for(i in tokens)
    {
        julia_assign("w1", i)
        for(j in candidates)
        {
            julia_assign("w2", j)
            # set value as -1, deg180 when no word is found in w2v model
            dists <- c(dists, 
                       tryCatch({ julia_eval("similarity(m, w1, w2)") },
                                error = function(err) { -1 },
                                finally = { -1 } )
            )			
        }
        similarities[[i]] <- dists
        dists <- NULL
    }
    similarities[['similarity']] <- Reduce(`+`, similarities)
    similarities[['candidates']] <- candidates
    return(similarities)
}


do.u.mean <- function(candidates, view.all = F)
{
    i <- length(candidates)
    if(view.all == F)
    {
        i <- 1
        pos <- c("", "y", "yes", "ok", "pos", "positive")
        neg <- c("n", "no", "neg", "negative")
        cancel <- c("c", "cancel", "stop", "break")
        while(i < length(candidates))
        {
            
            result <- readline(prompt=paste0("Do you mean '", candidates[i], "'? (y/n/c): "))
            if(any(pos %in% tolower(result)))
            {
                message(paste("User chose", candidates[i]))
                return(candidates[i])
            }
            if(any(cancel %in% tolower(result))) 
            {
                message("Warning: This word will be neglected.")
                return(c())
            }
            if(any(neg %in% tolower(result)))
            {
                
                i = i + 1
            }
        }
        message("No match found, this word will be neglected.")
        return(c())
    }
    else
    {
        message("Do You Mean?")
        message(paste0("[", 1:i, "] ", candidates, "\n"))
        result <- as.numeric(readline(prompt=paste0("Pleaes input a number: ")))
        if(!is.na(result) && result < length(candidates)) return(candidates[result]) else return(candidates[i])
    }
    
}


dedup.gene <- function (gene.list)
{
    genes <- c(gene.list$and, gene.list$not, gene.list$or, gene.list$ornot)
    dup.genes <- genes[duplicated(genes)]
    
    if(length(dup.genes) == 0) 
    {
        return(gene.list)
    }
    else
    {
        for(i in dup.genes)
        {
            if(length(grep(i, gene.list$genes)) > 1) message(paste0('Keeping ', grep(i, gene.list$genes, value = T)[1], ' in "', toString(grep(i, gene.list$genes, value = T)), '".'))
            gene.list[['genes']] <- setdiff(gene.list$genes, grep(i, gene.list$genes, value = T)[-1])
            return(gene.list)
        }
    }
    
}


weigh.gene <- function(page, inds, greedy)
{
    if(greedy == 1) return(strsplit(toString(page[inds,2]), ",|\\s")[[1]])
    genes <- strsplit(toString(page[inds,2]), ",|\\s")[[1]]
    freq <- as.numeric(strsplit(toString(page[inds,3]), ",|\\s")[[1]])
    cutoff <- if(greedy != 0) length(unique(freq)) * greedy else 1
    gene.list <- genes[freq %in% tail(unique(freq), cutoff)]
    
    if(length(gene.list) < length(genes)) message(paste0("Returned ", length(gene.list), " out of ", length(genes), " genes (argument greedy = ", greedy, ")"))
    
    return(gene.list)
}


token2phrases <- function( dictionary, token, pattern, spell.tolerate = F) ## need to consider ignore.case more in this function
{
    inds <- NULL
    candidates <- NULL
    bestmatch <- NULL
    sort.token <- paste(sort(strsplit(token, " |,|-")[[1]]), collapse=" ")
    all.token <- paste(token, sort.token, sep = "|")
    
    for ( i in pattern )
    {
        page <- match( i , names(dictionary))
        if(spell.tolerate == F)
        {
            tmp <- if(class(dictionary[[page]]) == "list") grep(all.token, dictionary[[page]][[1]][,1], ignore.case = T) else grep(all.token, dictionary[[page]][,1], ignore.case = T)
        }
        else 
        {
            tmp <- if(class(dictionary[[page]]) == "list") agrep(token, dictionary[[page]][[1]][,1], ignore.case = T) else agrep(token, dictionary[[page]][,1], ignore.case = T)
        }
        inds <- if(length(tmp) != 0) c(inds, paste0(i, "_", tmp)) else inds
        candidates <- if(class(dictionary[[page]]) == "list") c(candidates, as.character(dictionary[[page]][[1]][tmp,1])) else c(candidates, as.character(dictionary[[page]][tmp,1]))
    }
    
    df <- if(length(inds) != 0 ) data.frame("dictionary" = sub("\\$_.*|_.*", "", inds), "id" = as.numeric(sub(".*_", "", inds)), "phrase" = candidates, stringsAsFactors = F) else data.frame()
    
    if(any(df$phrase %in% c(token,sort.token))) 
    {
        df <- df[which(df$phrase %in% c(token,sort.token)), ] 
        df$phrase <- unname(sapply(df$phrase, function(x) {paste(sort(trimws(strsplit(df$phrase[1], ',')[[1]])), collapse=',')} ))
        if(any(duplicated(df$phrase)))
        {
            df <- head(df[order(match(df$dictionary, pattern)),], 1)
        }
        message(paste0( "Found '", token, "' in ", as.character(df$dictionary), " dictionary.") )
    }
    
    return(df)
}


id2genes <- function(dictionary, bestmatch, greedy)
{
    if(ncol(bestmatch) == 0) return(NULL)
    
    if(nrow(bestmatch) > 1) 
    {
        warning("argument 'bestmatch' has nrow > 1 and only the first row will be used")
        cat('\n')
        bestmatch <- bestmatch[1,]
    }
    
    page <- match(bestmatch$dictionary, names(dictionary) )
    
    # MESHID require 2 dictionaries
    if( class(dictionary[[page]]) == "list" )
    {
        mesh_id <- sub("\\s.*", "", as.character(dictionary[[page]][[1]][bestmatch$id,2]))
        mesh_id <- paste0(if( !any(grepl(":", mesh_id)) ) ":", mesh_id)
        
        # replace by row number in meshID2genename dictionary
        bestmatch$id <- grep(mesh_id, dictionary[[page]][[2]][,1])[1]
        if(is.na(bestmatch$id)) 
        {
            message(paste0("No relevant gene for '", bestmatch$phrase, "'.") )
            return(NULL)
        }
        else 
        {
            return(weigh.gene(page = dictionary[[page]][[2]], inds = bestmatch$id, greedy = greedy))
        }
        
    }
    
    # pairing row ID to gene names and return
    return(weigh.gene(page = dictionary[[page]], inds = bestmatch$id, greedy = greedy))
    
}


# scfindQ2ReadWord2Genes
read_dictionaries <-  function(paths) 
{	
    dictionary <- list()
    if(any(grepl(".rds$", paths, ignore.case = T)))
    {
        for( i in paths[grep(".rds$", paths, ignore.case = T)])
        {
            message(paste0("Reading '", sub(".*\\/", "", i ), "' file with size : ", floor(file.info(i)$size/1024/1024), " MB"))
            dictionary <- c(dictionary, readRDS(i))
        }
    }
    
    for(i in paths[grep(".rds$", paths, invert = T, ignore.case = T)])
    {
        name <-  gsub(".*/|\\..*", "", i)
        name <- gsub("_.*","" , name)
        eval(parse(text=paste0("dictionary$",name,"= read.delim2('",i,"', header=F)")))
    }
    
    return(dictionary)
}

# scfindQ2ReadWordVectors
read_w2v <- function(path)
{
    julia_install_package_if_needed("Word2Vec")
    julia_library("Word2Vec")
    model <- normalizePath(path) %>J% wordvectors(kind=binary)
    return(model)
}

# scfindQ2Wordcloud


gene2wordcloud <- function(dictionary, gene.list, word.window, not.genes = F, max.words = 100, return = F, page, any_id = c("MESH", "CHEBI", "OMIM"))
{
    word.window <- if(missing(word.window)) length(gene.list) * 10 else word.window# needs optimization
    result <- data.frame("phrases" = c(), "weight" = c())
    phrases <- c()
    weight <- NULL
    
    gene.list <- if(not.genes == F) grep("^-\\*|^\\*-|^-", gene.list, value = T, invert = T)
    gene.list <- gsub("^\\*|^-|^\\*-|^-\\*", "", gene.list)
    
    page <- if( missing (page) ) names(dictionary) else intersect(names(dictionary), page)
    all_ids <- paste0(any_id, ":") ## Tweak here if more ids are included
    
    message(paste("Calculating weight for", length(gene.list), "genes..."))
    for( i in grep("genename$", page, value = T, ignore.case=T) )
    {
        df <- NULL
        df.mesh <- NULL
        
        # exclude all phrases that do not include genes in gene.list
        if(class(dictionary[[i]]) != "list")
        {
            trim <- dictionary[[i]][-(grep(paste0(paste0(gene.list, ",|,", gene.list) , collapse="|"), dictionary[[i]][,2], invert=T, ignore.case=T)),]
        } else
        {
            trim <- dictionary[[i]][[2]][-(grep(paste0(paste0(gene.list, ",|,", gene.list) , collapse="|"), dictionary[[i]][[2]][,2], invert=T, ignore.case=T)),]
            
        }
        
        # create lists as the concept of hash table
        p2g <- setNames(strsplit(as.character(trim[,2]), ",", fixed=T), as.character(trim[,1]))
        p2w <- setNames(strsplit(as.character(trim[,3]), ",", fixed=T), as.character(trim[,1]))
        
        for( j in gene.list )
        {
            message(paste("Calculating weight for", j, "..."))
            df <- c(df, pair.id(p2w)[grep(paste0("#", j), pair.id(p2g), ignore.case = T)])
        }
        if(!is.null(df)){
            df <- data.frame("phrases" = gsub("#.*", "", df), "weight" = as.numeric(gsub(".*#", "", df)))
            df <- aggregate(df$weight, by=list(df$phrases), FUN=sum)
            
            mesh <- grep(paste(all_ids, collapse = "|"), df[,1], ignore.case = T)
            if(length(mesh) != 0)
            {
                
                best.mesh <- sub(".*:", "", df[tail(order(df[mesh,2]), 1),1])
                
                
                df.mesh <- dictionary[[i]][[1]][grep(best.mesh, dictionary[[i]][[1]][,2]),]
                
                df.mesh <- cbind(stack(setNames(strsplit(as.character(df.mesh[,2]), ",", fixed=T), as.character(df.mesh[,1]))),
                                 stack(setNames(strsplit(as.character(df.mesh[,3]), ",", fixed=T), as.character(df.mesh[,1]))) )
                
                df.mesh <- df.mesh[grep(best.mesh, df.mesh[,1], ignore.case=T),]
                df.mesh <- tail(df.mesh[order(as.numeric(df.mesh[,3])),c(2,3)], word.window)
                colnames(df.mesh) <- colnames(df)
                df <- df.mesh
            }
            
            phrases <- c(phrases, as.character(tail(df[order(df$x),1], word.window)))
            weight <- c(weight, tail(df[order(df$x),2], word.window))
        }
        
    }
    if(length(phrases) != 0)
    {
        result <- data.frame("phrases" = phrases, "weight" = as.numeric(weight) )
    }
    else
    {
        return(list())
    }
    
    if(return == T || max.words == 0) 
    {
        return(result)
    }
    else
    {
        wordcloud(phrases, rescale(as.numeric(weight), to = c(1, 100)), min.freq=1, max.words=max.words, random.order=F, rot.per=.35, colors=brewer.pal(8, "Dark2"))
    }
}

# scfindQ2Genes
query2genes <- function(object, dictionary, query, strict = F, automatch = T, greedy = 0.6, priority, spell.tolerate = T )
{
    model <- if(any(grepl("model", names(dictionary)))) dictionary[['model']] else warning('No word2vec model is provided in your dictionaries.')
    priority <- if(missing(priority)) NULL else priority
    query <- tokenize(object, query, strict)
    gene.list <- NULL 
    raw.genes <- NULL
    genes <- NULL
    phrase <- NULL
    weight <- NULL
    conds <- c("AND_", "OR_", "NOT_", "ORNOT_")
    symbl <- c("", "*", "-", "*-")
    
    results <- list()
    # if only gene query without operators
    if(identical(names(query),"AND_gene")) return(results <- ("genes" = query$AND_gene))
    
    
    
    if(any(grepl("_gene", names(query))))
    {
        gene.list <- c(
            if(!is.null(query$AND_gene)) query$AND_gene,
            if(!is.null(query$OR_gene)) paste0("*", query$OR_gene), 
            if(!is.null(query$NOT_gene)) paste0("-", query$NOT_gene),
            if(!is.null(query$ORNOT_gene)) paste0("*-", query$ORNOT_gene)
        )
        
        results <- list( "and" = if(!is.null(query$AND_gene)) query$AND_gene else c(),
                         "or" = if(!is.null(query$OR_gene)) query$OR_gene else c(),
                         "not" = if(!is.null(query$NOT_gene)) query$NOT_gene else c(),
                         "ornot" = if(!is.null(query$ORNOT_gene))  query$ORNOT_gene else c())
        # if only gene query
        if(all(grepl("_gene", names(query)))) 
        {
            results[["genes"]] <- gene.list	
            results <- dedup.gene(results)
            return(results[!sapply(results, is.null)])
        }
    }
    
    
    if(any(grepl("variant", names(dictionary), ignore.case = T)) &&  any(grepl("variant", names(query), ignore.case = T)))
    {
        page <- grep("variant", names(dictionary), ignore.case = T, value = T)
        
        for(i in grep("_variant", names(query), value = T))
        {
            raw.genes <- suppressMessages(caseCorrect(object, weigh.gene(page = dictionary[[page]], inds = which(dictionary[[page]][,1] %in% query[[i]]), greedy = greedy)))
            results[[tolower(sub("_variant", "", i))]] <- c(results[[tolower(gsub("_variant", "", i))]] , raw.genes)
            raw.genes <- if(!grepl('AND_variant', i) && length(raw.genes) != 0) paste0(symbl[grep(sub('variant', '', i), conds)], raw.genes) else raw.genes
            gene.list <- c(gene.list, if(length(raw.genes) != 0) raw.genes else NULL)
        }
        
        if(!any(grepl('_word|_mesh', names(query)))) 
        {
            results[["genes"]] <- unique(gene.list)
            results <- dedup.gene(results)
            return(results[!sapply(results, is.null)])
        }
    }
    
    if(any(grepl("meshID2genename", names(dictionary), ignore.case = T)) &&  any(grepl("mesh", names(query), ignore.case = T)))
    {
        
        for(page in grep("meshID2genename", names(dictionary), ignore.case = T))
        {
            for(i in grep("_mesh", names(query), value = T))
            {
                raw.genes <- suppressMessages(caseCorrect(object, weigh.gene(page = dictionary[[page]][[2]], inds = which(dictionary[[page]][[2]][,1] %in% toupper(query[[i]])), greedy = greedy)))
                
                results[[tolower(sub("_mesh", "", i))]] <- c(results[[tolower(gsub("_mesh", "", i))]] , raw.genes)
                raw.genes <- if(!grepl('AND_mesh', i) && length(raw.genes) != 0) paste0(symbl[grep(sub('mesh', '', i), conds)], raw.genes) else raw.genes
                gene.list <- c(gene.list, if(length(raw.genes) != 0) raw.genes else NULL)
            }
        }
        
        if(!any(grepl('_word', names(query)))) 
        {
            results[["genes"]] <- unique(gene.list)
            results <- dedup.gene(results)
            return(results[!sapply(results, is.null)])
        }
    }
    
    
    df <- stack(query)
    df <- df[grep("_word", df[,2]),]
    df <- subset(df, df$values != "")
    
    # clean stop words
    if(any(grepl("stop", names(dictionary), ignore.case = T)))
    { 
        stop.inds <- which(df[,1] %in% tolower(dictionary[[grep("stop", names(dictionary), ignore.case = T)]][,1]))
        df <- if(length(stop.inds) != 0) df[-stop.inds,] else df
    }
    
    pattern = c(priority, setdiff(dictionary$priority, priority))
    
    if(nrow(df) != 0 && any(grepl(paste(pattern, collapse = "|"), names(dictionary), ignore.case = T)))
    {
        bestmatch <- data.frame("dictionary" = c(), "id" = c(), "phrase" = c(), stringsAsFactors = F)
        
        
        for(i in unique(grep("_word", df[,2],value = T)))
        {
            
            tokens <- subset(df, df$ind == i)[,1]
            tokens <- strsplit(paste(tokens, collapse=" "), " # |#| #|# ")[[1]]
            tokens <- tokens[tokens != ""]
            if(length(tokens) != 0)
            {
                for ( j in tokens )
                {
                    res <- token2phrases(dictionary = dictionary, token = j, pattern = pattern, spell.tolerate = spell.tolerate )
                    
                    if(nrow(res) == 1)
                    {
                        # handling best match
                        raw.genes <- c(raw.genes, id2genes(dictionary = dictionary, bestmatch = res, greedy = greedy))
                        message(paste("Found", toString(raw.genes), "for", j))
                    }
                    else 
                    {
                        if(nrow(res) == 0)
                        {
                            # handling no phrase found by sorting words in alphabetical order and search against sorted database
                            tmp.tk <- strsplit(j, " ")[[1]]
                            
                            if(length(tmp.tk) == 1)
                            {
                                # if not match any with the single word
                                raw.genes <- NULL
                                message(paste0("No gene found for '", j,"'" ))
                            }
                            else
                            {
                                # narrow down dictionary
                                cands <- token2phrases(dictionary = dictionary, token = gsub("\\s","|",j), pattern = pattern)
                                # sort in alphabetical order for comparison
                                cands$phrase <-  unname(sapply(as.character(cands$phrase), function(x)paste(sort(unlist(strsplit(x, " "))), collapse=" ")))
                                
                                for( k in sort(tmp.tk))
                                {
                                    cands.tk <- grep(k, cands$phrase, value = T, ignore.case = T)
                                    
                                    ### compare every cands.tk to bestmatch, if same, skip, else find bestmatch 
                                    if(length(cands.tk) == 0) 
                                    {
                                        message(paste0("No gene found for '", k,"'" ))
                                    }
                                    else 
                                    {
                                        if(!any(unname(sapply(bestmatch$phrase, function(y)paste(sort(unlist(strsplit(y, " "))), collapse=" "))) %in% cands.tk))
                                        {
                                            cands.tk <- stack(setNames(strsplit(cands.tk, " |-|,"), cands.tk))  # for words has hivens and commas
                                            cands.sim <- cos.sim.Word2Vec(model, strsplit(j, " ")[[1]], unique(cands.tk[,1]))
                                            cands.tk$similarity <- cands.sim$similarity[match(cands.tk[,1], cands.sim$candidates)]
                                            cands.tk <- aggregate(cands.tk$similarity, by=list(cands.tk$ind), FUN=mean)
                                            cands.tk <- rev(as.character(cands.tk[order(cands.tk[,2]),1]))
                                            
                                            bestmatch <- rbind(bestmatch, cands[which(cands$phrase %in% cands.tk[1]),])
                                        }
                                    }
                                    
                                }
                                message(paste0("Found '", paste(bestmatch$phrase, collapse= "','"), "' for '", j, "'."))
                            }
                        }
                        else
                        {
                            # handling very non-specific query
                            if(length(j) == 1 && spell.tolerate == T)
                            {
                                # correct token for better match
                                tk.in.common <- tail(sort(table(strsplit(paste(res$phrase, collapse = " "), " |-|,")[[1]])),1)
                                j <- if(agrepl(j, names(tk.in.common)) && tk.in.common[[1]] == nrow(res)) names(tk.in.common) else j
                            }
                            
                            cands.tk <- setNames(strsplit(res$phrase, " |-|,"), res$phrase) # for words has hivens and commas
                            
                            if(length(cands.tk) > 1)
                            {
                                cands.tk <- stack(cands.tk)
                                cands.sim <- cos.sim.Word2Vec(model, strsplit(j, " ")[[1]], unique(cands.tk[,1]))
                                cands.tk$similarity <- cands.sim$similarity[match(cands.tk[,1], cands.sim$candidates)]
                                cands.tk <- aggregate(cands.tk$similarity, by=list(cands.tk$ind), FUN=mean)
                                cands.tk <- rev(as.character(cands.tk[order(cands.tk[,2]),1]))
                                
                                if(automatch == T) 
                                {	
                                    bestmatch <- res[which(res$phrase %in% head(cands.tk[1]))[1],]
                                    message(paste0("Found '", paste(bestmatch$phrase, collapse= "','"), "' for ", j, ".")) 
                                }
                                else 
                                {
                                    bestmatch <- res[which(res$phrase %in% do.u.mean(cands.tk, F))[1],]
                                }
                            }
                            else
                            {
                                bestmatch <- res
                            }
                            
                        }
                        
                        for( k in 1: nrow(bestmatch) ) 
                        {
                            raw.genes <- c(raw.genes, id2genes(dictionary = dictionary, bestmatch = bestmatch[k,], greedy = greedy))
                        }
                        
                    }
                }
                
                genes <- suppressMessages(caseCorrect(object, raw.genes))
                results[[tolower(sub("_word", "", i))]] <- c(results[[tolower(sub("_word", "", i))]], genes)
                genes <- if(length(genes) != 0) paste0(symbl[match(gsub("word", "", i), conds)], genes) else NULL
                gene.list <- if(!is.null(genes)) c(gene.list, genes) else gene.list
            }
            
        }
        results[["genes"]] <- unique(gene.list)
        results <- dedup.gene(results)
        return(results[!sapply(results, is.null)])
    }
    else
    {
        results[["genes"]] <- if(length(gene.list) != 0) gene.list else NULL
        results <- dedup.gene(results)
        return(results[!sapply(results, is.null)])
    }
    
}

# scfindQ2CellTypes

query2CellTypes <- function(object, dictionary, query, datasets, optimize = T, abstract = T, strict = F, greedy = 0.6, priority, spell.tolerate = T )
{
    
    # Run hyperQueryCellTypes according to the gene list generated from free form search
    datasets <- if(missing(datasets)) object@datasets else select.datasets(object, datasets)
    optimized.query <- NULL
    result <- list()
    gene.list <- query2genes(object = object, dictionary = dictionary, query = query, strict = strict, greedy = greedy, priority = priority, spell.tolerate = spell.tolerate )
    
    
    if(optimize == T && any(grepl("and|or", names(gene.list))))
    {
        for(i in grep("^and$|^or$", names(gene.list)))
        {
            if(length(gene.list[[i]]) > 1) 
            {
                optimized.query <- suppressMessages(markerGenes(object, gene.list[[i]])) 
                optimized.query <- toString(optimized.query$Query[tail(order(optimized.query$tfidf),1)])
                optimized.query <- strsplit(optimized.query, ",")[[1]]
            }
            else
            {
                optimized.query <- gene.list[[i]]
            }
            
            gene.list [['genes']] <- setdiff(gene.list[['genes']], paste0( if(names(gene.list[i]) == "or") "*" else "", setdiff(gene.list[[i]], optimized.query)))
            gene.list[[i]] <- if(length(optimized.query) != 0) optimized.query else gene.list[[i]]
        }
        result[['query_optimized']] <- gene.list
    }
    else
    {
        result[['query']] <- gene.list
    }
    
    
    cell_hits <- suppressMessages( findCellTypes(object, gene.list$genes, datasets) )
    
    if(length(cell_hits) != 0) 
    {
        q2CT <- phyper.test(object, cell_hits, datasets)
    } 
    else 
    {
        message("No Cell Is Found!")
        return(result)
    }
    message(paste0("Found ", nrow(q2CT), " celltypes for your search '", sub(',.*', '', query), "...'."))
    if(abstract == F)
    {
        q2CT <- q2CT[order(q2CT[,"pval"]),]
        result[['cell_hits']] <- cell_hits
        result[['celltypes']] <- q2CT
        return(result)
    }
    else
    {
        significant.celltype <- q2CT[which(q2CT$pval <= min(q2CT$pval, 0.001)),]
        if(nrow(significant.celltype) != 0) 
        {
            result[['celltypes_significant']] <- significant.celltype
            cell_hits <- cell_hits[as.character(significant.celltype$cell_type)]
            result[['cell_hits_significant']] <- cell_hits
            message(paste0("Found ", nrow(significant.celltype), " celltypes has pval <= ", min(q2CT$pval, 0.001), "."))
        }
        else
        {
            message("No Significant CellType Is Found!")
        }
        
    }
    return(result)
}

# Prepare for the drag and drop building index feature
hackCellTypeIndex.SCESet <- function(matrix, annotation, qb = 2)
{
    
    if (grepl(dataset.name,'.'))
    {
        stop("The dataset name should not contain any dots")
    }
    
    
    cell.types.all <- as.factor(annotation$cell_type1)
    cell.types <- levels(cell.types.all)
    new.cell.types <- hash(keys = cell.types, values = paste0(dataset.name, '.', cell.types))
    genenames <- unique(rownames(matrix))
    
    if (length(cell.types) > 0)
    {
        non.zero.cell.types <- c()
        index <- hash()

        exprs <- matrix
        
        ef <- new(EliasFanoDB)
        
        qb.set <- ef$setQB(qb)
        if (qb.set == 1)
        {
            stop("Setting the quantization bits failed")
        }
        for (cell.type in cell.types) {
            inds.cell <- which(cell.type == cell.types.all)
            if(length(inds.cell) < 2)
            {
                ## print(paste('Skipping', cell.type))
                next
            }
            non.zero.cell.types <- c(non.zero.cell.types, cell.type)
            message(paste("\tIndexing", cell.type, "as", new.cell.types[[cell.type]], " with ", length(inds.cell), " cells."))
            cell.type.exp <- exprs[,inds.cell]
            if(is.matrix(exprs))
            {
                ef$indexMatrix(new.cell.types[[cell.type]], cell.type.exp)
            }
            else
            {
                ef$indexMatrix(new.cell.types[[cell.type]], as.matrix(cell.type.exp))
            }
        }
    }
    
    
    index <- new("SCFind", index = ef, datasets = dataset.name)
    
    
    return(index)
}

