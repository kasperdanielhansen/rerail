library(Matrix)
randomSVD <- function(A, p = 2, k = 100,
                      method = c("svd", "qr")) {
    method <- match.arg(method)
    m <- nrow(A)
    n <- ncol(A)
    l <- n+p
    G <- Matrix(rnorm(n*l), ncol = l, nrow = n) # n-l
    H <- A %*% (crossprod(A) %*% G) # m - l 
    if(method == "svd") {
        sv0 <- svd(crossprod(H)) # l - l 
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

