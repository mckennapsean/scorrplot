[**s-CorrPlot**: Visualizing Correlation](http://mckennapsean.github.io/scorrplot)
==================================================================================

The **s-CorrPlot** is a new scatterplot for visually exploring pairwise correlation coefficients between all variables in large datasets.



About the **s-CorrPlot**
------------------------

The degree of correlation between variables is used in many data analysis applications as a key measure of similarity. The most common techniques for visualizing correlation, like scatterplot matrices and clustered heatmaps, however, do not scale well to large datasets, either computationally or visually. We present a new visualization that is capable of encoding pairwise correlation between hundreds of thousands variables, called the **s-CorrPlot**. The **s-CorrPlot** is based on a 2D scatterplot and exploits the geometric structure underlying Pearson’s correlation to derive a novel spatial encoding. The **s-CorrPlot** not only depicts a visually precise measure of correlation, but also supports visualizing metadata using encoding channels like color. We implemented the **s-CorrPlot** as an open-source proof of concept visualization in order to validate its effectiveness through a variety of methods including a case study with a biology collaborator.


This proof of concept employs simple multidimensional exploration techniques, to demonstrate how this visual encoding can employ the vast set of these exploration techniques for exploring correlation.

For further details, please read the description and derivation of the **s-CorrPlot** in our paper.



Code Dependencies
-----------------

The **s-CorrPlot** has been implemented within a proof of concept prototype: *Gyroscope*. The proof of concept is integrated within the R statistical framework for data input and output. The *Gyroscope* code is written in C++ and OpenGL.

The source code has been compiled and tested on both Mac OS X and Linux.

To compile *Gyroscope* through R, the following components must be installed:

-  **R** ([version 2.10 or newer](http://www.r-project.org/))
-  **Xcode** (for *Mac*, [install Xcode](https://developer.apple.com/xcode/) and [download command line tools](https://developer.apple.com/support/xcode/))
-  **Fortran compiler** (for *Mac*, [update your version](http://cran.r-project.org/bin/macosx/tools/))
-  **OpenGL** (for *Linux*, OpenGL library [like freeglut](http://freeglut.sourceforge.net/))



Installation
------------

Be sure to install all dependencies, as detailed above.

Then, simply run the install script as root:

    sudo ./install

If it prints out "Done", then *Gyroscope* has installed correctly inside R, as the library *"gyroscope"*.



Visualizing Correlation
-----------------------

If you do not need a fully interactive prototype, you can ignore the above dependencies and installation process. Instead, simply create static plots of the s-CorrPlot by using the provided **s-CorrPlot**.R function within the *code/* folder.

To explore correlation using the prototype, you can load the data from the paper in our proof of concept visualization.
Sample scripts load the data in R, and the data is provided in R format so that you can input your own data, too.

Each script loads a different dataset, corresponding to figures from the paper. The variables being correlated are listed first in the description:

-  **genes**.R
   -  genes in two brain regions, *Figure 6*
-  **geneDensity**.R
   -  genes in several regions, *Figure 3*
-  **imagePatches**.R
   -  subset of image patches of a full image, *Figure 1(c)*
-  **imagePatches-full**.R
   -  complete image patches for two image datasets, *Figure 1(b)*
-  **subway-stops**.R
   -  stations of subway ridership, *Figure 4(b)*
-  **subway-time**.R
   -  years of subway ridership, *Figure 4(a)*

Certain datasets have been anonymized in order to protect our biology collaborator's sensitive data.

These scripts can be run from terminal or loaded in R:

    ./genes.R

*or*

    R
    source("genes.R")

Please note that these scripts cannot be run from the R GUI program; they must be executed from the terminal.

For further instructions on how to use *Gyroscope*, please read the [code documentation](http://mckennapsean.github.io/scorrplot/documentation.html).



Uninstallation
--------------

Run the uninstall script as root or start R to remove the *Gyroscope* package:

    sudo ./uninstall

*or*

    R
    remove.packages("gyroscope")



License
-------

This project's codebase *(Gyroscope)* is licensed by GPLv2.



Authors
-------

[Sean McKenna](http://www.seanpmckenna.com/), [Miriah Meyer](http://www.cs.utah.edu/~miriah/), [Christopher Gregg](http://www.neuro.utah.edu/people/faculty/gregg.html), & [Samuel Gerber](http://www.math.duke.edu/~sgerber/)

The *Gyroscope* code was originally designed and developed by [Samuel Gerber](http://www.math.duke.edu/~sgerber/).



Contact
-------

If you have any difficulties or questions, please contact <sean@cs.utah.edu>.
