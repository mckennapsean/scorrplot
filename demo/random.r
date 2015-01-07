# create a uniform random dataset
# 10,000 variables with 10 observations each
m = matrix(runif(10 * 10000), ncol = 10, nrow = 10000)
data = data.frame(m)

# add random labels (as factors!)
data$labels = as.factor(ceiling((0.1 + runif(10000)) * 4))
levels(data$labels) = c("A", "B", "C", "D", "E")
names = c();
for(i in 1 : 10000)
  names[i] = paste("name", i, sep = "-")

rownames(data) = names

# start up the s-CorrPlot interactive display
# disable density estimation since random data is spread out
scorr(data, alpha = 0.2, useDensity = F)

# get the number of variables or points on the screen
scorr.get.size()

# select a new variable of interest
scorr.set.primary(data[3, 1 : ncol(data) - 1]);

# query the top 100 correlations in our current projection
cor = scorr.get.cor(100)

# highlight our first ten variables in the plot
scorr.highlight.index(1 : 10)

# highlight some other variables by their names
scorr.highlight.name(c("name-45", "name-101"))

# create a static plot of our current view
# must have called scorr() first in order to use the get functions!
scorr.plot(data, scorr.get.primary(), scorr.get.secondary(), 0.2)
