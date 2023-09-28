#!/usr/bin/Rscript --vanilla

library( 'randomForest' )

values = list()

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
    values[[method]] = data.frame( mutation = D[!is_train,1], value = P )
}

merged = merge( values[['our']], values[['mutabind2']], by='mutation' )
names(merged) <- c( 'mutation', 'our', 'mutabind2' )

svg( '/dev/stdout' )
plot( merged[,3], merged[,2], pch='x',
      xlab = 'ddG values of model trained on MutaBind2 terms',
      ylab = 'ddG values of model trained on our terms' )

x = c( 1, 2, 3 )
y = c( 1, 2, 3 )

diagonal = lm( y ~ x )
fit = lm( merged[,2] ~ merged[,3] )

abline( fit, col = 'red' )
abline( diagonal, col='black' )

legend( 'bottomright', c( 'y = x', 'trend' ), col = c( 'black', 'red' ), lty = c( 1, 1 ) )
