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

D = read.csv( 'train-dataset-mutabind2.tab', sep="\t", quote='' )
experimental = data.frame( mutation = D$mutation, ddG = D$ddG )

merged = merge( merge( values[['our']], values[['mutabind2']], by='mutation' ), experimental, by='mutation' )
names(merged) <- c( 'mutation', 'our', 'mutabind2', 'experimental' )

x = c( 1, 2, 3 )
y = c( 1, 2, 3 )

diagonal = lm( y ~ x )

svg( '/dev/stdout' )

par( mfrow = c( 3, 3 ) )
for (i in 2:4) {
    for (j in 2:4) {
        plot( merged[,j], merged[,i], xlab = names(merged)[j], ylab = names(merged)[i], pch = 'x' )

        fit = lm( merged[,i] ~ merged[,j] )
        abline( fit, col = 'red' )
        abline( diagonal, col='black' )
    }
}
