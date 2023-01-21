Rcpp::sourceCpp("ZIE.cpp")



M <- rbind(g1,g2,g3) |> scale()
group <- as.factor(c(rep("g1", 400), rep("g1", 400), rep("g2", 200))) |> as.integer()


RM <- matrix(c(rnorm(800), rnorm(200)+runif(200,5,6)), 500, 2) |> scale()
group <- c(rep(1, 400), rep(2, 100))



plot(M, col=group)
res = nearest_diff_class(M, group)
drop(res$first_diff_class)
drop(res$nn_idx+1)
mask <- drop(res$first_diff_class) == res$first_diff_class[res$nn_idx+1]
pre <- M[drop(res$first_diff_class == 20 & mask == TRUE & group == "g1"),]
M[drop(res$first_diff_class == 20 & mask == TRUE & group == "g1"),] |>
  points(col="blue", cex=1.4)


points(x$centroids, col="green", cex=1.5, pch=19)
x <- mean_shift(data, 20)
most_frequent_row(x)


table(res$first_diff_class[group==2])
M[drop(res$first_diff_class == 17 & group==1),] |> points(col="blue", cex=1.4)
ms <- mean_shift(x, 15)


pre <- RM[drop(res$first_diff_class == 36 & group==1),]
x <- mean_shift(pre, 15)
res4$data |> points(col="blue", cex=1.4)
most_frequent_row2(x$centroids)

res2 <- nearest_diff_class(res$data, res$labels)
res3 <- nearest_diff_class(res2$data, res2$labels)
res4 <- nearest_diff_class(res3$data, res3$labels)
res5 <- nearest_diff_class(res4$data, res4$labels)
res6 <- nearest_diff_class(res5$data, res5$labels)
res7 <- nearest_diff_class(res6$data, res6$labels)
res8 <- nearest_diff_class(res7$data, res7$labels)
res9 <- nearest_diff_class(res8$data, res8$labels)
res10 <- nearest_diff_class(res9$data, res9$labels)
res11 <- nearest_diff_class(res10$data, res10$labels)
res12 <- nearest_diff_class(res11$data, res11$labels)

plot(res12$data, col=res12$labels)
points(res9$data, col="blue", cex=1.4)

irism <- as.matrix(iris[,1:4])
irisl <- as.integer(as.factor(iris[,5]))

pca <- prcomp(irism)$x[,1:2]
plot(as.matrix(pca), col=irisl)
