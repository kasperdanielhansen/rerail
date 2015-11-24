library(Matrix)
randomSVD <- function(A, p = 2, k = 100,
                      method = c("svd", "qr")) {
    method <- match.arg(method)
    m <- nrow(A)
    n <- ncol(A)
    l <- n+p
    G <- Matrix(rnorm(n*l), ncol = l, nrow = n)
    H <- A %*% (crossprod(A) %*% G)
    if(method == "svd") {
        sv0 <- svd(crossprod(H))
        omega <- diag(1/sqrt(sv0$d))
        Q <- H %*% sv0$u %*% omega
        T <- crossprod(Q, A)
    }
    if(method == "qr") {
        Q <- qr.Q(qr(H,0))
        T <- crossprod(Q, A)
    }
    sv <- svd(T)
    v <- sv$v[, 1:k]
    d <- sv$d[1:k]
    u <- Q %*% sv$u[,1:k]
    return(list(u = u, d = d, v = v))
}

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
