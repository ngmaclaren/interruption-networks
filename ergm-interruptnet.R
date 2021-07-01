library(ergm)
library(ergm.count)
set.seed(12345)

### Uncomment one
loc <- './data/networks-iss/'
## loc <- './data/networks-nss/'
## loc <- './data/networks-both/'

ext <- '.gml'
graph_files <- list.files(loc, pattern = paste0(ext, "$"))

graphs <- list()
for (file in graph_files) {
    file_loc <- paste(loc, file, sep = "")
    g <- igraph::read_graph(file_loc, format = gsub("\\.", "", ext))
    g <- intergraph::asNetwork(g)
    gID <- gsub(ext, "", file)
    network::set.vertex.attribute(g, 'gID', gID)
    graphs[[gID]] <- g
}

sm_list <- list()
node_attr_list <- list()
i <- 1
for (g in graphs) {
    tm <- as.sociomatrix(g, attrname = 'weight')
    rownames(tm) <- g %v% 'label'
    colnames(tm) <- g %v% 'label'
    tdf <- data.frame(node = g %v% 'label',
                      gID = g %v% 'gID',
                      gender = g %v% 'gender',
                      tst = g %v% 'tst',
                      esl = g %v% 'esl',
                      isop = g %v% 'isop',
                      intel = g %v% 'intel',
                      gameknowl = g %v% 'gameknowl',
                      inst = g %v% 'inst',
                      sim = g %v% 'sim',
                      consc = g %v% 'consc',
                      agree = g %v% 'agree',
                      neur = g %v% 'neur',
                      open = g %v% 'open',
                      extra = g %v% 'extra')
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
node_attr$tst <- node_attr$tst/1000

smc <- sm
rownames(smc) <- seq(1:nrow(smc))
colnames(smc) <- seq(1:ncol(smc))
mmg <- as.network(smc, directed = TRUE, matrix.type = 'a', ignore.eval = FALSE, names.eval = 'weight')
set.vertex.attribute(mmg, 'label', node_attr$node)
set.vertex.attribute(mmg, 'gID', node_attr$gID)
set.vertex.attribute(mmg, 'gender', node_attr$gender)
set.vertex.attribute(mmg, 'tst', node_attr$tst)
set.vertex.attribute(mmg, 'esl', node_attr$esl)
set.vertex.attribute(mmg, 'isop', node_attr$isop)
set.vertex.attribute(mmg, 'intel', node_attr$intel)
set.vertex.attribute(mmg, 'gameknowl', node_attr$gameknowl)
set.vertex.attribute(mmg, 'inst', node_attr$inst)
set.vertex.attribute(mmg, 'sim', node_attr$sim)
set.vertex.attribute(mmg, 'consc', node_attr$consc)
set.vertex.attribute(mmg, 'agree', node_attr$agree)
set.vertex.attribute(mmg, 'neur', node_attr$neur)
set.vertex.attribute(mmg, 'open', node_attr$open)
set.vertex.attribute(mmg, 'extra', node_attr$extra)

bigmodel_notst <- ergm(
    mmg ~ sum + mutual + transitiveweights +
        nodeicov('intel') +
        nodeicov('gameknowl') +
        nodeicov('consc') + nodeicov('agree') + nodeicov('neur') +
        nodeicov('open') + nodeicov('extra') +
        nodeifactor('gender') + nodeifactor('esl') + nodeifactor('isop') +
        nodeocov('intel') +
        nodeocov('gameknowl') +
        nodeocov('consc') + nodeocov('agree') + nodeocov('neur') +
        nodeocov('open') + nodeocov('extra') +
        nodeofactor('gender') + nodeofactor('esl') + nodeofactor('isop') +
        nodematch('gender') + nodematch('esl')
   ,
    coef = -1, reference = ~ DiscUnif(0, 10),
    response = 'weight', constraints = ~ blockdiag('gID'))
summary(bigmodel_notst)

bigmodel_triads <- ergm(
    mmg ~ sum + mutual + transitiveweights +
        nodeicov('tst') +
        nodeicov('intel') +
        nodeicov('gameknowl') +
        nodeicov('consc') + nodeicov('agree') + nodeicov('neur') +
        nodeicov('open') + nodeicov('extra') +
        nodeifactor('gender') + nodeifactor('esl') + nodeifactor('isop') +
        nodeocov('tst') +
        nodeocov('intel') +
        nodeocov('gameknowl') +
        nodeocov('consc') + nodeocov('agree') + nodeocov('neur') +
        nodeocov('open') + nodeocov('extra') +
        nodeofactor('gender') + nodeofactor('esl') + nodeofactor('isop') +
        nodematch('gender') + nodematch('esl')
   ,
    coef = -1, reference = ~ DiscUnif(0, 10),
    response = 'weight', constraints = ~ blockdiag('gID'))
summary(bigmodel_triads)
