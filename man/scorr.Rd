\name{scorr}
\alias{scorr}
\alias{scorr.set.primary}
\alias{scorr.set.secondary}
\alias{scorr.get.primary}
\alias{scorr.get.secondary}
\alias{scorr.get.size}
\alias{scorr.get.selected}
\alias{scorr.get.density}
\alias{scorr.get.corr}
\alias{scorr.get.acorr}
\alias{scorr.get.cor}
\alias{scorr.get.acor}
\alias{scorr.get.name}
\alias{scorr.highlight.index}
\alias{scorr.highlight.name}
\alias{scorr.close}
\alias{scorr.plot}
\title{s-CorrPlot: Visualizing Correlation}

\description{

The s-CorrPlot is an interactive scatterplot for visually exploring pairwise
correlation coefficients between variables in large datasets. Variables are
projected as points on a scatterplot with respect to some user-selected
variables of interest, driven by a geometric interpretation of correlation.
The correlation of all other variables to the selected one is indicated by
vertical gridlines in the plot. By selecting new variables of interest, a
user can create simple tours of the correlation space through animations
between different projections of the data.

For further details about the s-CorrPlot, please refer to our companion paper.



\bold{Tool Overview:}

The call \code{scorr(data)} will initialize the s-CorrPlot tool with your input
dataset. In this view, the default projection is based off of the first two
principal component vectors of the dataset, but you can select any other
variables in the plot to animate to a new projection of the data. These
principal components are illustrated as interactive bar-charts on the left of
the screen. A dataset must provide factor labels for the entire dataset, and
each unique label is given its own color as shown in the boxes along the bottom.
Lastly, on the right is both a parallel coordinate plot showing a profile of
the observations for each variable as well as a text display for all variables
near the cursor in the display.



\bold{Data Format:}

The s-CorrPlot expects a data frame where each row is your variable and its
columns correspond to different observations, where \code{rownames} and
\code{colnames} correspond to these respectively. Additionally, a final column
of nominal labels is expected (the factor data type), so that each variable
gets an associated color in the plot. Each row variable will become a point
in the final s-CorrPlot display, and the observations will be displayed for
these points using the parallel coordinates view. For further details, please
see the package's provided data \code{data(NAME)} as well as the data frame
constructed in the example.



\bold{Mouse Interactions:}

Hovering the mouse over the data points will update the parallel coordinate
plot and text display views by the variables near the cursor. Additionally,
you can left-click variables or the bar-charts to select them as the new primary
variable of interest, thus animating the plot accordingly and this point moves
to the right-most part of the circular plot. All correlation coefficients to
this primary variable can be read directly by vertical gridlines, where anti-
correlated variables are furthest away and zero correlation is near the middle.

Additionally, you can right-click the bars or points in order to select them as
the new secondary variable of interest, and then the s-CorrPlot will animate the
projection by moving all variables vertically to their new locations. It is also
possible to shift-click on variables in order to affix them in the parallel
coordinates plot for quick comparison against other variables.



\bold{Keyboard Interactions:}

To analyze correlation with Spearman's rank correlation, press 's'. You can
return to the default of Pearson's correlation by pressing 'p'.

To remove all highlighting and selected variables in the display, press 'c'.

To adjust transparency of data points, press the left and right arrow keys.

To adjust the size of data points, press the up and down arrow keys.

To adjust the kernel bandwidth for density estimation, press '<' or '>'.

To change the background color of the plot, press 'b'.



\bold{Queries through R console:}

Once you have initialized the s-CorrPlot tool, there are several commands to
interact with the display from the R console.

The primary or secondary projection can be set directly by passing in a vector
with the observations stored in the columns (or length equal to the number of
observations in the input data). These can be input using
\code{scorr.set.primary()} and \code{scorr.set.secondary()}. Additionally,
these same vectors can be output to R using the associated get functions.

The values used for density estimation in the current view of the tool can be
output back into R using the call \code{scorr.get.density()}.

To retrieve the currently n top correlation or anti-correlated points (i.e.
left-most and right-most points, respectively), call \code{scorr.get.cor(n)} or
\code{scorr.get.acor(n)}, which returns a data frame with n rows with columns corresponding to each variable's index in the input data, its correlation to the primary variable, and its rowname. Similarly, you can call either
\code{scorr.get.corr(r)} or \code{scorr.get.acorr(r)} in order to output all
variables with above or below the given correlation coefficient, respectively.

To highlight a set of variables, pass a list of indices to highlight into the
command \code{scorr.highlight.index()} or a list of names to
\code{scorr.highlight.name()}.

The call \code{scorr.close()} will close the interactive tool safely.

The function \code{scorr.plot()} takes data and two variables as input in order
to produce a static version of the s-CorrPlot, in case the interactive tool
is unable to build or run on a machine. This is useful for Windows computers,
since they are unable to compile the tool.

}


\usage{
  scorr(data,
        useDensity = T,
        showProfile = T,
        showPatch = F,
        threshold = 0,
        coloring = 1 : nrows(data),
        perms = 0
  )
  scorr.set.primary(v)
  scorr.set.secondary(v)
  scorr.get.primary()
  scorr.get.secondary()
  scorr.get.size()
  scorr.get.selected()
  scorr.get.density()
  scorr.get.corr(r)
  scorr.get.acorr(r)
  scorr.get.cor(n)
  scorr.get.acor(n)
  scorr.get.name(indices)
  scorr.highlight.index(indices)
  scorr.highlight.name(names)
  scorr.close()
  scorr.plot(data2, v2, v2)
}

\arguments{
  \item{data}{A data frame with d columns and n rows. Each row is a variable
and the first 1:(d-1) columns are observations and the nth column is a
label assignment in the form of a factor. Names for the variables should be
supplied in \code{rownames(data)}, and names for each column can be supplied in
\code{colnames(data)} }
  \item{useDensity}{Set transparancey of the plotted variables based on a
density estimate of the distribution along the correlation sphere.}
  \item{showProfile}{Show the parallel coordinates plot.}
  \item{showPatch}{Show each variable as a greysacle patch image.}
  \item{threshold}{Variables with variance below the threshold are excluded.}
  \item{coloring}{Alternate continuous values for each variable for coloring
the data instead of using labels (switch by pressing c).}
  \item{perms}{Number of permutations for permutation test of projection. The
permutation test computes the probability that there is a projection with more
points with correlation r>a for a>=0 and r<a for a<=0 than the current shown
projection.}
  \item{v}{Vector with d-1 columns in order to set a new projection.}
  \item{r}{Correlation coefficient to threshold which variables to output.}
  \item{n}{Number of top anti-/correlations to retrieve.}
  \item{indices}{Variable indices to retrieve the names of, likely in order
to highlight them in the s-CorrPlot display.}
  \item{names}{Data point names to highlight in the s-CorrPlot display}
  \item{data2}{A data frame where rows correspond to variables. May include the
final column as factor labels for automatic coloring but not required.}
  \item{v2}{Both a primary and secondary vector for defining a static s-CorrPlot
view. Can include labels from the rows or pipe in from the interactive display.}
}
\value{
  A vector in correlation space is returned for calling both primary/secondary.
  All combinations of corr/acorr and cor/acor return a data frame with indices, correlation coefficients, and variable names for valid variables.
  Calling size returns the number of variables or points on the display.
  Calling selected returns the index of the primary variable of interest.
  Calling density returns all the density estimate values used for the current
  calculation of transparency in the view.
  Calling name returns a list of names corresponding to the input indices.
  All other methods return nothing (NULL).
}

\author{
  Samuel Gerber & Sean McKenna
}
\seealso{
}
\references{
http://mckennapsean.com/projects/s-corrplot/
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
\keyword{multivariate,models,nonlinear,nonparametric,iplot}
