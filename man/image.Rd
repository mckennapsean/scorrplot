\docType{data}
\name{image}
\alias{image}
\title{Image Patches of a Full Image}
\format{A data frame with 131072 rows and 82 columns}
\usage{
data(image)
}
\description{
A dataset containing each image patch as a row or variable
with the columns corresponding to the pixels of the patch.
The entire dataset contains the values of the grayscale image
for that patch split up across the columns. There is also a
final column dictating which image that patch corresponds to;
it is possible to show multiple datasets simultaneously using
colors in the s-CorrPlot display. This is the format of a data frame
for input to the scorrplot package.
}
\source{http://graphics.cs.williams.edu/data/images.xml#2}
\keyword{datasets}

