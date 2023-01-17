Rcpp::sourceCpp("ZIE.cpp")



M <- rbind(g1,g2,g3) |> scale()
group <- as.factor(c(rep("g1", 20), rep("g1", 20), rep("g2", 20)))







plot(M, col=group)
res = nearest_diff_class(M, group)
drop(res$first_diff_class)
drop(res$nn_idx+1)
mask <- drop(res$first_diff_class) == res$first_diff_class[res$nn_idx+1]
pre <- M[drop(res$first_diff_class == 20 & mask == TRUE & group == "g1"),]
M[drop(res$first_diff_class == 20 & mask == TRUE & group == "g1"),] |>
  points(col="blue", cex=1.4)


points(x$centroids, col="green", cex=1.5, pch=19)
most_frequent_row(x)


