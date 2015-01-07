data(subway)
transpose = as.data.frame(t(data[, 1 : ncol(data) - 1]))
transpose$labels = factor("subway")
transpose = scorr(transpose, alpha = 0.2, useDensity = F)
