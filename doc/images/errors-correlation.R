#!/usr/bin/Rscript --vanilla

library( 'randomForest' )

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
    P = predict(model, newdata=test)
    errors[[method]] = data.frame( mutation = D[!is_train,1], error = P - test$ddG )
}

merged = merge( errors[['our']], errors[['mutabind2']], by='mutation' )
names(merged) <- c( 'mutation', 'our', 'mutabind2' )

cor.test(D[,2], D[,3])
