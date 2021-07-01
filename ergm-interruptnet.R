##library(statnet)
library(ergm)
library(ergm.count)
set.seed(12345)

### Uncomment one
loc <- './data/networks-iss/'
## loc <- './data/networks-nss/'
## loc <- './data/networks-both/'

ext <- '.gml'
graph_files <- list.files(loc, pattern = paste0(ext, "$"))
##graph_files <- graph_files[grep(ext, graph_files)]
graphs <- list()
for (file in graph_files) {
    file_loc <- paste(loc, file, sep = "")
    g <- igraph::read_graph(file_loc, format = gsub("\\.", "", ext))
    g <- intergraph::asNetwork(g)
    gID <- gsub(ext, "", file)
    network::set.vertex.attribute(g, 'gID', gID)
    graphs[[gID]] <- g
}

                                        # By-hand correction
                                        # igraph will not import the try int(); except str() solution for "age" that worked in NetworkX
                                        # removed age as node attribute in generate-networks.py
##set.vertex.attribute(graphs[["XSP"]], "age", NA, 6) # "XSP07"
##delete.vertices(graphs[["XSP"]], 6) # to use age

### Isolated nodes appear to have an effect on ERGM results. To recover significance for the `mutual` coefficient, remove isolates:
## for (i in 1:length(graphs)) {
##     isolates <- get.vertex.attribute(graphs[[i]], "vertex.names")[which(has.edges(graphs[[i]]) == FALSE)]
##     delete.vertices(graphs[[i]], isolates)
## }

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

## the as.network function isn't working: not retaining edge weights. 
#mmg <- network(sm, matrix.type = 'adjacency', ignore.eval = FALSE, edge.check = TRUE) # for mega-matrix graph
smc <- sm
rownames(smc) <- seq(1:nrow(smc))
colnames(smc) <- seq(1:ncol(smc))
mmg <- as.network(smc, directed = TRUE, matrix.type = 'a', ignore.eval = FALSE, names.eval = 'weight')#, ignore.eval = FALSE, names.eval = 'weight')
set.vertex.attribute(mmg, 'label', node_attr$node)
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

############################
## Begin Test Models Here ##
############################

### transitive and intransitive are not implemented for weighted ERGM
## interesting options:
#### diff('tst', dir = 'h-t')
#### what's the difference between idegree() and nodeifactor()?
#### nodeicov('tst'), nodeocov('tst')
#### nodemix('gender') instead of nodematch('gender')

## Ok.
## Need to redo the AIC analysis, including the node attributes for ESL and ISOP.
## Justification for including these and no others are that these are obvious, salient features
## that may play into the probability that someone may interrupt someone else. 
## Intelligence is not obvious, although perhaps it too should be included.
## And the game knowledge score. Shit.

## bigmodel <- ergm(
##     mmg ~ sum + mutual +
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
##         nodematch('gender') + nodematch('esl') # diff() or absdiff() on 'intel' or 'gameknowl' etc.
##    ,
##     coef = -1, reference = ~ DiscUnif(0, 10),
##     response = 'weight', constraints = ~ blockdiag('gID'))
## summary(bigmodel)$aic # old -5389.614; with nodematch stats -5384.146
## coefs <- summary(bigmodel)$coefs
## coefs <- coefs[order(coefs[4]),]
## coefs$alpha <- 0.05/(nrow(coefs) + 1 - seq(1, nrow(coefs)))
## coefs$test <- coefs[4] < coefs[5]
## coefs[,4:6]

bigmodel_notst <- ergm(
    mmg ~ sum + mutual + transitiveweights +
        ####nodeicov('tst') + #nodeicov('age') +
        nodeicov('intel') +
        nodeicov('gameknowl') +
        nodeicov('consc') + nodeicov('agree') + nodeicov('neur') +
        nodeicov('open') + nodeicov('extra') +
        nodeifactor('gender') + nodeifactor('esl') + nodeifactor('isop') +
        ##nodeifactor('inst') + nodeifactor('sim') +
        ####nodeocov('tst') + #nodeocov('age') +
        nodeocov('intel') +
        nodeocov('gameknowl') +
        nodeocov('consc') + nodeocov('agree') + nodeocov('neur') +
        nodeocov('open') + nodeocov('extra') +
        nodeofactor('gender') + nodeofactor('esl') + nodeofactor('isop') +
        ##nodeofactor('inst') + nodeofactor('sim')
        nodematch('gender') + nodematch('esl')
   ,
    coef = -1, reference = ~ DiscUnif(0, 10),
    response = 'weight', constraints = ~ blockdiag('gID'))
summary(bigmodel_notst)

bigmodel_triads <- ergm(
    mmg ~ sum + mutual + transitiveweights +
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
        nodematch('gender') + nodematch('esl') # diff() or absdiff() on 'intel' or 'gameknowl' etc.
   ,
    coef = -1, reference = ~ DiscUnif(0, 10),
    response = 'weight', constraints = ~ blockdiag('gID'))
summary(bigmodel_triads) # 

## m1 <- ergm(
##     mmg ~ sum + mutual,
##     coef = -1, reference = ~ DiscUnif(0, 10),
##     response = 'weight', constraints = ~ blockdiag('gID'))

## m2 <- ergm(
##     mmg ~ sum + mutual +
##         diff('tst', dir = 'h-t'),
##     coef = -1, reference = ~ DiscUnif(0, 10),
##     response = 'weight', constraints = ~ blockdiag('gID'))

## m3 <- ergm(
##     mmg ~ sum + mutual +
##         diff('tst', dir = 'h-t') + nodeicov('tst'),
##     coef = -1, reference = ~ DiscUnif(0, 10),
##     response = 'weight', constraints = ~ blockdiag('gID'))

## m4 <- ergm(
##     mmg ~ sum + mutual +
##         diff('tst', dir = 'h-t') + nodeicov('tst') +
##         nodematch('gender'),
##     coef = -1, reference = ~ DiscUnif(0, 10),
##     response = 'weight', constraints = ~ blockdiag('gID'))

## m5 <- ergm(
##     mmg ~ sum + mutual +
##         diff('tst', dir = 'h-t') +
##         nodeicov('tst') +
##         nodematch('gender') +
##         nodeifactor('gender') + nodeofactor('gender'),
##     coef = -1, reference = ~ DiscUnif(0, 10),
##     response = 'weight', constraints = ~ blockdiag('gID'))

## m6 <- ergm(
##     mmg ~ sum + mutual +
##         ##diff('tst', dir = 'h-t') +
##         nodeicov('tst') +
##         nodeocov('tst') + 
##         nodematch('gender') +
##         nodeifactor('gender') + nodeofactor('gender'),
##     coef = -1, reference = ~ DiscUnif(0, 10),
##     response = 'weight', constraints = ~ blockdiag('gID'))
## salientmodel <- ergm(
##     mmg ~ sum + mutual +
        
## notst <- ergm(
##     mmg ~ sum + mutual +
##         nodematch('gender') +
##         nodeifactor('gender') + nodeofactor('gender'),
##     coef = -1, reference = ~ DiscUnif(0, 10),
##     response = 'weight', constraints = ~ blockdiag('gID'))
## models <- list(m1, m2, m3, m4, m5, notst)
## AIC <- sapply(models, function (x) summary(x)$aic)
## AIC <- AIC - min(AIC)
## names(AIC) <- c('m1', 'm2', 'm3', 'm4', 'm5', 'notst')
## AIC

## ## coefficient values from the best model
## bestmodel <- models[[which(AIC == min(AIC))]]
## summary(bestmodel)
## exp(summary(bestmodel)$coefs$Estimate)
## summary(notst)
## exp(summary(notst)$coefs$Estimate)


#####################################
## Distribution Check - DO NOT CUT ##
#####################################


## ## Check that DiscUnif approximation is appropriate
## actual <- mmg %e% 'weight'
## nszs <- # for non-structural zeroes
##     lapply(graphs, function (x) {
##         m <- as.sociomatrix(x, attrname = 'weight')
##         s <- sum(colSums(m == 0)) - nrow(m)
##         return(s)
##     })
## nsz <- rep(0, Reduce('+', nszs))
## actual <- c(nsz, sort(actual))
## nsims <- length(actual)
## y <- network.initialize(2, directed = TRUE)
## trgeo <- as.vector(simulate(y ~ sum, coef = -1.1, reference = ~ DiscUnif(0, max(actual)), response = "weight", output = "stats", nsim = nsims))
## geo <- as.vector(simulate(y ~ sum, coef = -1.25, reference = ~ Geometric, response = 'weight', output = 'stats', nsim = nsims))
## trgeo_df <- data.frame(x = trgeo, src = 'trgeo')
## geo_df <- data.frame(x = geo, src = 'geo')
## actual_df <- data.frame(x = actual, src = 'actual')
## plotting_data <- bind_rows(trgeo_df, geo_df, actual_df)
## ggplot(plotting_data, aes(x = x, group = src, fill = src)) + # , y = ..count../sum(..count..)
##     geom_histogram(binwidth = 1, alpha = .5, position = 'dodge') +
##     facet_wrap('src')








## xlims <- c(0, max(actual))
## colors = c(rgb(1, 0, 0, 1/4), rgb(0, 1, 0, 1/4), rgb(0, 0, 1, 1/4))
## breaks = 30

## trgeo_hist <- hist(trgeo,
##                    breaks = breaks,
##                    probability = TRUE)#,
##                    ##freq = FALSE,
##                    ##plot = FALSE)
## geo_hist <- hist(geo,
##                  breaks = breaks,
##                  density = TRUE)#,
##                  ##freq = FALSE,
##                  ##plot = FALSE)
## actual_hist <- hist(actual,
##                     breaks = breaks,
##                     density = TRUE)#,
##                     ##freq = FALSE,
##                     ##plot = FALSE)

## ## normalize <- function (h) {
## ##     h$counts <- h$counts/sum(h$counts)
## ## }

## plot(trgeo_hist, col = colors[1], xlim = xlims)
## plot(geo_hist, col = colors[2], xlim = xlims, add = TRUE)
## plot(actual_hist, col = colors[3], xlim = xlims, add = TRUE)
## legend('topright', pch = 22, col = colors, pt.bg = colors, legend = c('TrGeo', 'Geo', 'Actual'))


##############
## OLD CODE ##
##############

## ## mcmc.diagnostics(test_model)
## ## ## gof(ergmFitObject, GOF=~model) #### run help(gof)
## ## summary(test_model)

## ## next model needs to be a valued ergm
## mmgm_count <- ergm(mmg ~
##                        ## reviewer asked for: homophily, preferential attachment (?), transitivity, reciprocity (pretty sure that's mutual?)
##                        ## can look at nodematch() for homophily
##                        ## can look at intransitive, but doesn't seem to be implemented in valued ERGMs; there is a transitive statistic as well
##                        sum +
##                        mutual + #(form = 'geometric') +, larger AIC
##                        ##nodeicov('tst') +
##                        ##nodeocov('tst') +
##                        ##CMP +
##                        nodeifactor('gender') +
##                        nodeofactor('gender'),
##                    coef = -1,
##                    reference = ~ DiscUnif(0, 10), ## this works and covers over 95% of edge values. It is, however, not correct. It may do as an approximation for now. AIC are very high. Investigate model appropriateness, but probably not before the end of the month.
##                    ##reference = ~ Geometric, # Poisson crashes the sess, Geometric says it's not implemented
##                    response = 'weight',
##                    constraints = ~ blockdiag('gID')#,
##                    #MCMC.prop.weights = 'random'
##                    )
## summary(mmgm_count) 

