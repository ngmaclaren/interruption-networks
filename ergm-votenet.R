library(statnet)

loc <- './data/votenets/'
ext <- '.gml'
graph_files <- list.files(loc)
graph_files <- graph_files[grep(ext, graph_files)]
graphs <- list()
for (file in graph_files) {
    file_loc <- paste(loc, file, sep = "")
    g <- igraph::read_graph(file_loc, format = gsub("\\.", "", ext))
    g <- intergraph::asNetwork(g)
    gID <- gsub(ext, "", file)
    set.vertex.attribute(g, 'gID', gID)
    graphs[[gID]] <- g
}

sm_list <- list()
node_attr_list <- list()
i <- 1
for (g in graphs) {
    tm <- as.sociomatrix(g)
    rownames(tm) <- g %v% 'label'
    colnames(tm) <- g %v% 'label'
    tdf <- data.frame(node = g %v% 'label',
                      gID = g %v% 'gID',
                      gender = g %v% 'gender')
    sm_list[[i]] <- tm
    node_attr_list[[i]] <- tdf
    i <- i + 1
}
    
starter <- sm_list[[1]]
for (j in 2:length(sm_list)) {
    nextone <- merge(starter, sm_list[[j]], by = 'row.names', all = TRUE)
    nextone <- as.matrix(nextone[-1])
    rownames(nextone) <- colnames(nextone)
    starter <- nextone
}

sm <- starter
sm[is.na(sm)] <- 0

node_attr <- do.call('rbind', node_attr_list)
    
mmg <- as.network(sm) # for mega-matrix graph
set.vertex.attribute(mmg, 'gID', node_attr$gID)
set.vertex.attribute(mmg, 'gender', node_attr$gender)

mmgm <- ergm(mmg ~
                 edges +
                 mutual +
                 nodeifactor('gender') +
                 nodeofactor('gender'),
             constraints = ~ blockdiag('gID')
             )
summary(mmgm)
