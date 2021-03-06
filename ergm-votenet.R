library(statnet)
set.seed(12345)

loc <- './data/votenets/'
ext <- '.gml'
graph_files <- list.files(loc, pattern = paste0(ext, "$"))
graphs <- list()
for (file in graph_files) {
    file_loc <- paste(loc, file, sep = "")
    g <- igraph::read_graph(file_loc, format = gsub("\\.", "", ext))
    g <- intergraph::asNetwork(g)
    gID <- gsub(ext, "", file)
    set.vertex.attribute(g, 'gID', gID)
    graphs[[gID]] <- g
}

                                        # By-hand correction
##set.vertex.attribute(graphs[["XSP"]], "age", NA, 6) # "XSP07"

### Isolated nodes appear to have an effect on ERGM results. To recover significance for the `mutual` coefficient, remove isolates:
## for (i in 1:length(graphs)) {
##     isolates <- get.vertex.attribute(graphs[[i]], "vertex.names")[which(has.edges(graphs[[i]]) == FALSE)]
##     delete.vertices(graphs[[i]], isolates)
## }

sm_list <- list()
node_attr_list <- list()
i <- 1
for (g in graphs) {
    tm <- as.sociomatrix(g)#, attrname = 'weight')
    rownames(tm) <- g %v% 'label'
    colnames(tm) <- g %v% 'label'
    tdf <- data.frame(node = g %v% 'label',
                      gID = g %v% 'gID',
                      gender = g %v% 'gender',
                      tst = g %v% 'tst',
                      esl = g %v% 'esl',
                      isop = g %v% 'isop',
                      ##age = g %v% 'age',
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

mmg <- as.network(sm) # for mega-matrix graph
## set.vertex.attribute(mmg, 'label', node_attr$node) ### Why is this missing? Investigate.
set.vertex.attribute(mmg, 'gID', node_attr$gID)
set.vertex.attribute(mmg, 'gender', node_attr$gender)
set.vertex.attribute(mmg, 'tst', node_attr$tst)
set.vertex.attribute(mmg, 'esl', node_attr$esl)
set.vertex.attribute(mmg, 'isop', node_attr$isop)
##set.vertex.attribute(mmg, 'age', node_attr$age)
set.vertex.attribute(mmg, 'intel', node_attr$intel)
set.vertex.attribute(mmg, 'gameknowl', node_attr$gameknowl)
set.vertex.attribute(mmg, 'inst', node_attr$inst)
set.vertex.attribute(mmg, 'sim', node_attr$sim)
set.vertex.attribute(mmg, 'consc', node_attr$consc)
set.vertex.attribute(mmg, 'agree', node_attr$agree)
set.vertex.attribute(mmg, 'neur', node_attr$neur)
set.vertex.attribute(mmg, 'open', node_attr$open)
set.vertex.attribute(mmg, 'extra', node_attr$extra)

## bigmodel <- ergm(
##     mmg ~ edges + mutual +
##         ##diff('tst', dir = 'h-t') +
##         nodeicov('tst') + #nodeicov('age') +
##         nodeicov('intel') +
##         nodeicov('gameknowl') +
##         nodeicov('consc') + nodeicov('agree') + nodeicov('neur') +
##         nodeicov('open') + nodeicov('extra') +
##         nodeifactor('gender') + nodeifactor('esl') + nodeifactor('isop') +
##         ##nodeifactor('inst') +
##         ##nodeifactor('sim') +
##         nodeocov('tst') + #nodeocov('age') +
##         nodeocov('intel') +
##         nodeocov('gameknowl') +
##         nodeocov('consc') + nodeocov('agree') + nodeocov('neur') +
##         nodeocov('open') + nodeocov('extra') +
##         nodeofactor('gender') + nodeofactor('esl') + nodeofactor('isop') +
##         ##nodeofactor('inst')## + nodeofactor('sim')
##         nodematch('gender') + nodematch('esl')
##    ,
##     constraints = ~ blockdiag('gID'))

bigmodel_notst <- ergm(
    mmg ~ edges + mutual + transitive + intransitive +
        ##diff('tst', dir = 'h-t') +
        ##nodeicov('tst') + #nodeicov('age') +
        nodeicov('intel') +
        nodeicov('gameknowl') +
        nodeicov('consc') + nodeicov('agree') + nodeicov('neur') +
        nodeicov('open') + nodeicov('extra') +
        nodeifactor('gender') + nodeifactor('esl') + nodeifactor('isop') +
        ##nodeifactor('inst') +
        ##nodeifactor('sim') +
        ##nodeocov('tst') + #nodeocov('age') +
        nodeocov('intel') +
        nodeocov('gameknowl') +
        nodeocov('consc') + nodeocov('agree') + nodeocov('neur') +
        nodeocov('open') + nodeocov('extra') +
        nodeofactor('gender') + nodeofactor('esl') + nodeofactor('isop') +
        ##nodeofactor('inst')## + nodeofactor('sim')
        nodematch('gender') + nodematch('esl')
   ,
    constraints = ~ blockdiag('gID'))
summary(bigmodel_notst)

bigmodel_triads <- ergm(
    mmg ~ edges + mutual + transitive + intransitive +
        ##diff('tst', dir = 'h-t') +
        nodeicov('tst') + #nodeicov('age') +
        nodeicov('intel') +
        nodeicov('gameknowl') +
        nodeicov('consc') + nodeicov('agree') + nodeicov('neur') +
        nodeicov('open') + nodeicov('extra') +
        nodeifactor('gender') + nodeifactor('esl') + nodeifactor('isop') +
        ##nodeifactor('inst') +
        ##nodeifactor('sim') +
        nodeocov('tst') + #nodeocov('age') +
        nodeocov('intel') +
        nodeocov('gameknowl') +
        nodeocov('consc') + nodeocov('agree') + nodeocov('neur') +
        nodeocov('open') + nodeocov('extra') +
        nodeofactor('gender') + nodeofactor('esl') + nodeofactor('isop') +
        ##nodeofactor('inst')## + nodeofactor('sim')
        nodematch('gender') + nodematch('esl')
   ,
    constraints = ~ blockdiag('gID'))
summary(bigmodel_triads)


## Insert AIC analysis here.

## m1 <- ergm(
##     mmg ~ edges + mutual,
##     constraints = ~ blockdiag('gID'))
## m2 <- ergm(
##     mmg ~ edges + mutual + transitive,
##     constraints = ~ blockdiag('gID'))
## m3 <- ergm(
##     mmg ~ edges + mutual + transitive + intransitive,
##     constraints = ~ blockdiag('gID'))
## m4 <- ergm(
##     mmg ~ edges + mutual + transitive + intransitive + nodeicov('tst'),
##     constraints = ~ blockdiag('gID'))
## m5 <- ergm(
##     mmg ~ edges + mutual + transitive + intransitive + nodeicov('tst') + diff('tst', dir = 'h-t'),
##     constraints = ~ blockdiag('gID'))
## m6 <- ergm(
##     mmg ~ edges + mutual + transitive + intransitive + nodeicov('tst') + diff('tst', dir = 'h-t') +
##         nodematch('gender'),
##     constraints = ~ blockdiag('gID'))
## m7 <- ergm(
##     mmg ~ edges + mutual + transitive + intransitive + nodeicov('tst') + diff('tst', dir = 'h-t') +
##         nodematch('gender') + nodeifactor('gender') + nodeofactor('gender'),
##     constraints = ~ blockdiag('gID'))
## notst <- ergm(
##     mmg ~ edges + mutual + transitive + intransitive +
##         nodematch('gender') + nodeifactor('gender') + nodeofactor('gender'),
##     constraints = ~ blockdiag('gID'))
## notriads <- ergm(
##     mmg ~ edges + mutual +
##         nodeicov('tst') + diff('tst', dir = 'h-t') +
##         nodematch('gender') + nodeifactor('gender') + nodeofactor('gender'),
##     constraints = ~ blockdiag('gID'))
## notriadsnotst <- ergm(
##     mmg ~ edges + mutual +
##         nodematch('gender') + nodeifactor('gender') + nodeofactor('gender'),
##     constraints = ~ blockdiag('gID'))
## models <- list(m1, m2, m3, m4, m5, m6, m7, notst, notriads, notriadsnotst)
## AIC <- sapply(models, function (x) summary(x)$aic)
## AIC <- AIC - min(AIC)
## names(AIC) <- c('m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'notst', 'notriads', 'notriadsnotst')
## AIC

## coefs_to_prob <- function(model) {
##     coefs <- summary(model)$coefs$Estimate
##     probs <- exp(coefs)/(1 + exp(coefs))
##     probs
## }

## bestmodel <- models[[which(AIC == min(AIC))]]
## summary(bestmodel)
## exp(summary(bestmodel)$coefs$Estimate)
## summary(notst)
## exp(summary(notst)$coefs$Estimate)
## summary(notriads)
## exp(summary(notriads)$coefs$Estimate)

## mmgm <- ergm(mmg ~
##                  edges +
##                  mutual +
##                  nodeifactor('gender') +
##                  nodeofactor('gender'),
##              constraints = ~ blockdiag('gID')
##              )
## summary(mmgm)
