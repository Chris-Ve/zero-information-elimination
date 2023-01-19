Rcpp::sourceCpp("ZIE.cpp")



M <- rbind(g1,g2,g3) |> scale()
group <- as.factor(c(rep("g1", 20), rep("g1", 20), rep("g2", 20))) |> as.integer()


RM <- matrix(c(rnorm(800), rnorm(200)+runif(200,5,6)), 500, 2) |> scale()
group <- c(rep(1, 400), rep(2, 100))



plot(RM, col=group)
res = nearest_diff_class(RM, group)
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
RM[drop(res$first_diff_class == 2 & group==2),] |> points(col="blue", cex=1.4)
ms <- mean_shift(x, 15)


pre <- RM[drop(res$first_diff_class == 36 & group==1),]
x <- mean_shift(pre, 15)
x$centroids |> points(col="green", pch=19)
most_frequent_row2(x$centroids)
