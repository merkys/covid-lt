#!/usr/bin/Rscript --vanilla

library( 'randomForest' )
library( 'BlandAltmanLeh' )

errors = list()

for (method in c( 'our', 'mutabind2' )) {
    D = read.csv( paste( 'train-dataset-', method, '.tab', sep='' ), sep="\t", quote='' )
    train_size = round( length(D[,1]) * 0.8 )

    set.seed( 1410 )
    
    is_train = logical( length(D[,1]) )
    is_train[1:train_size] = TRUE
    is_train = sample( is_train )

    train = D[ is_train, 2:length(D)]
    test  = D[!is_train, 2:length(D)]

    names(train) = names(D)[2:length(D)]
    names(test)  = names(D)[2:length(D)]

    load( paste( 'binding-evaluator-model-', method, '.RData', sep='' ) )
    errors[[method]] = data.frame( mutation = D[!is_train,1], ddG = predict(model, newdata=test) )
}

D = read.csv( paste( 'train-dataset-mutabind2.tab', sep='' ), sep="\t", quote='' )
D = data.frame( mutation = D$mutation, ddG = D$ddG )
merged = merge( merge( D, errors[['our']], by='mutation' ), errors[['mutabind2']], by='mutation' )
names(merged) <- c( 'mutation', 'true', 'our', 'mutabind2' )

svg( '/dev/stdout' )
BlandAltmanLeh::bland.altman.plot( merged$our, merged$mutabind2,
                                   xlab = 'average of estimated ddG',
                                   ylab = 'difference of estimated ddG',
                                   pch = 'x' )
