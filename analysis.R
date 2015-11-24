library(Matrix)
load("objects/cpM.rda")

library(lineprof)

newM <- cpM[1:3000,1:3000]
prof <- lineprof({
system.time({
    
    tmp=randomSVD(newM)
})


system.time({
    tmp <- svd(newM)
})



## crossprod(X,Y) = t(X) %*% Y

testA <- Matrix(rgamma(1000*300, shape = 3), ncol = 300, nrow = 1000)
sv.test <- svd(testA)
tmp = randomSVD(testA, method = "qr")

all.equal(sv.test$d[1:100], tmp$d)
all.equal(as.numeric(sv.test$u[, 1:100]), as.numeric(tmp$u))
all.equal(as.numeric(sv.test$v[, 1:100]), as.numeric(tmp$v))

debugonce(randomEigen)
tmp$u[1:10,1:10] - sv.test$u[1:10,1:10]
tmp$u[1:10,1] - sv.test$u[1:10,1]
tmp$u[1:10,2]
tmp$v[1:10,1:10] - sv.test$v[1:10,1:10]
