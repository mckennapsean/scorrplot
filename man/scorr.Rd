\name{scorr}
\alias{scorr}
\alias{scorr.highlight.index}
\alias{scorr.highlight.name}
\alias{scorr.get.selected}
\alias{scorr.get.cor}
\alias{scorr.get.acor}
\alias{scorr.get.name}
\alias{scorr.close}
\title{s-CorrPlot: Orthogonal Projections for Visualizing Correlation}
\description{
This package provides an interactive visualization tool to explore correlation
structures in data sets.

The visualization is based on 2-dimensional orthogonal projection of the
z-scores of each sample (data point). This allows various interpretation of the
projected point. Most importantly for points directly on the unit circle the
correlations to all other is exactly visualized through distances. For details
of the orthogonal projection please have a look at the paper in the reference
section.



\bold{Display Description:}

The call \code{scorr} starts up a display window with a scatter plot of the
orthogonal projection of the z-scores. The horizontal lines correspond to 0.9,
0.8, ..., 0.1, 0, -0.1, ..., -0.9 correlation compared to the most left point on
the unit circle. For each point a label has to be assigned. The boxes on the
bottom show the color for each label. The bar plots on the left hand show the
principal component standard deviations for the complete data sets and for each
group of points with the same label. On the right is a parallel coordinate plot
that can show the profile of selected points.  


\bold{Mouse Interactions:}

By moving the mouse over data points a list of data point names appears on the
right hand side. These are the points that are in the vicinity of the cursor.  By
clicking with the mouse on a data point the projection plane is oriented such
that the point is in the most left point on the unit circle and the horizontal
lines correspond to correlations to that point. And the text list on the left
hand side is fixed when the cursor is not over any data points to the points in
the neighborhood of the last click. By clicking on a name from the text list the
plane is oriented through that data point. The projection plane can also be
oriented to any of the principal components by clicking on the corresponding bar
on the left hand side.

Since a plane through the origin is fully described through two (non collinear)
data points an additional point can be selected to lie on the unit circle of the
projection. This is achieved with right mouse clicks and can be done as for the
normal mouse click on data points, bars or the data point name in text list.

With shift click on a data point or a data point name the data point is
permanently display in the parallel coordinates plot and highlighted in
the projection plane. 

By clicking on the axis on the parallel coordinates plot the points are colored
according to the values in that coordinate. By clicking on multiple coordinates
an the average values of those coordinates, multiple coordinates can be added
and removed by clicking on them. To facilitate comparison between two sets of
coordinate colorings ctrl-click creates another coloring by averaging the
coordinate values selected through ctrl-click.  

Clicking on the label boxes will highlight or un-highlight the set of points of
that label.

\bold{Keyboard Interactions:}

The projection can be changed to Spearman's rank correlation by pressing 's' and
to Pearson's correlation by 'p' (default). 

By pressing 'c' removes all highlighting and data points in the parallel
coordinates plot.

Pressing the left and right arrows, decreases and increases transparency of the
data points, respectively.

By pressing 'b', the background color can be changed from the default dark to a
bright white background, for the purposes of taking screenshots from the program
for editing or placing as a figure into a paper.

By pressing 'a', the text labels for each item can be toggled on and off for the
color red. This allows for easier visual searching of non-red items.

The coloring scheme can be changed through pressing '1' fro label coloring, '2'
for coloring set one of average coordinate values (selected through clicking on
coordinates in parallel coordinate plot) and '3' for coloring set two of average
coordinate values (ctr-click in parallel coordinates plot). 



\bold{Queries through R console:}

To start up the display call \code{scorr}. s-CorrPlot expects a data frame
with each row a data point. Names for each data points are extract through the
\code{rownames} of the data frame and coordinate names through \code{colnames}. The last
column needs to be a factor which indicates for each data point to which labels
it belongs and the names of the labels correspond to the levels of that factor.
For details see the setup of the data frame in the example section. 

Once the s-CorrPlot is started there are several commands to interact with the display
through the R console. 

To retrieve the currently n top correlation or anti-correlated points (i.e. left
most and right most points, respectively).  call \code{scorr.get.cor(n)} or
\code{scorr.get.acor(n)} which will return a data frame with n rows with the
first column the index of the data point into the data frame passed to the s-CorrPlot,
the second row the correlation coefficients to vector corresponding to the left
most point of the unit circle and the last row the name of the data point.

To highlight a set of genes pass a list of indices (into the data frame) of the
data points to highlight through \code{scorr.highlight.index} or a list of
names through \code{scorr.highlight.name}

To set primary or secondary projection projection directy pass a vector of
lenngth equal the number of observation for each data point to
\code{scorr.set.primary} or \code{scorr.set.secondary}, respectively.

}


\usage{
  scorr(data, perms = 0, threshold = 0, coloring = 1:nrows(data), useDensity
      = T, showProfile = T, showPatch = F)
  scorrehighlight.index(indices)
  scorr.highlight.name(names)
  scorr.get.size()
  scorr.get.selected()
  scorr.get.density()
  scorr.get.corr(r)
  scorr.get.acorr(r)
  scorr.get.cor(n)
  scorr.get.acor(n)
  scorr.get.name(indices)
  scorr.set.primary(v),
  scorr.set.secondary(v),
  scorr.close()
}

\arguments{
  \item{data}{A data frame with n columns and m rows. Each row is a data point
(sample) and the first 1:(n-1) columns are observations and the nth column is a
label assignment in the form of a factor. Names for the data points should be
supplied in \code{rownames(data)} and names each column can be supplied in
\code{colnames(data)} }
  \item{perms}{Number of permutations for permutation test of projection. The
permutation test computes the probability that there is a projection with more
points with correlation r>a for a>=0 and r<a for a<=0 than the current shown
projection.}
  \item{n}{Number of top anti-/correlations to retrieve}
  \item{indices}{Data point indices to retrieve names for or to highlight in the s-CorrPlot
display} 
  \item{names}{Data point names to highlight in the s-CorrPlot display}
  \item{threshold}{Data points with varance below threhold are excluded}
  \item{coloring}{Alternate continuous values for rach data  point for coloring
    the data instead of lables (swicth by pressing c)}
  \item{useDensity}{Set transparancey based of plotted points based on a density
    estimate of the distribution on the correlation sphere}
   \item{showProfile}{Show parallel coordinates plot}
  \item{showPatch}{Show data point as a greysacle patch image} 
}
\value{
  For cor/acor the top n correlated/anti-correlated data points to the current
first projection direction (left most point of the unit circle in the display).
  All other return methods  nothing (NULL).
}

\author{
  Samuel Gerber
}
\seealso{
}
\references{
In preparation
}
\examples{
library(scorr)

####
#Create uniform random data sets with 10 dimensions and 10000 points
####

m <- matrix(runif(10*10000), ncol=10, nrow=10000);
data <- data.frame(m)

#random labels
#(labels need to be in the last column of the data frame)
data$labels = as.factor(ceiling( (0.1 + runif(10000)) * 4) )
#names for labels
levels(data$labels) = c("A", "B", "C", "D", "E");
#names for data points
names = c();
for(i in 1:10000){ names[i] = paste("name", i, sep="-")}
rownames(data) = names;


####
#use the s-CorrPlot to explore data
####

#start up interactive display
scorr(data)

#get the number of points on the screen
scorr.get.size()

#query current selected data point index
sel <- scorr.get.selected()

#query all points for their density estimation values
#must have used the density estimation for this command
dens <- scorr.get.density()

#query the positive correlations between 0.95 and 1 to the current projection
corr <- scorr.get.corr(0.95)

#query the negative correlations between -0.95 and -1 to the current projection
acorr <- scorr.get.acorr(0.95)

#query top 100 correlations to the current projection
cor <- scorr.get.cor(100)

#query top 100 correlations to the current projection
acor <- scorr.get.acor(100)


#highlight points 1 through 10
scorr.highlight.index(1:10);

#highlight by names
scorr.highlight.name(c("name-45", "name-101"))

#make sure to close the s-CorrPlot before closing R
#( you can also close the display by pressing 'q' in the s-CorrPlot
#  display )
scorr.close();


}
\keyword{nonparametric,models,nonlinear}
