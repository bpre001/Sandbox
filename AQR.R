library(pracma)

A <- matrix(1:9,
            nrow = 3,
            byrow = FALSE)

A
R0 <- rref(A)
R0

decomp <- qr(A)
Q <- qr.Q(decomp)
R <- qr.R(decomp)
Q
R

QQT <- Q %*% t(Q)
QTQ <- t(Q) %*% Q
QTQ

cat("A:\n"); print(A)
cat("\nR0:\n"); print(R0)

cat("\nQ:\n"); print(-Q)
cat("\nR:\n"); print(-R)


cat("\nVerify A = Q %*% R:\n"); print(round(Q %*% R, 10))
