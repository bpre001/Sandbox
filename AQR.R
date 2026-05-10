library(pracma)

A <- matrix(c(1, 2, 3,
              4, 5, 6), nrow = 2, byrow = FALSE)

R0 <- rref(A)

decomp <- qr(A)
Q <- qr.Q(decomp)
R <- qr.R(decomp)

QQT <- Q %*% t(Q)
QTQ <- t(Q) %*% Q
QTQ

cat("A:\n"); print(A)
cat("\nR0:\n"); print(R0)

cat("\nQ:\n"); print(Q)
cat("\nR:\n"); print(R)


cat("\nVerify A = Q %*% R:\n"); print(round(Q %*% R, 10))
