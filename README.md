# afblue

This is an R-package that makes it easy to compute the BLUE for estimating 
allele frequency amongst the founders given the pedigree connecting a sample.

It also provides an algorithm for determining optimal sibling removal schemes
and it makes it easy to compute effective sample sizes.

You can get the package thus:
```r
devtools::install_github("eriqande/afblue")
```

However, if you want to see how the package was used in 
my short piece with Robin Waples, then, after you have done the above line,
you should get the whole
repository which includes the scripts.  Do it like this:
```sh 
git clone https://github.com/eriqande/afblue
```
Then open up the Rstudio project in that repository, and run the code in 
`R-main/01-predict-ess.R` (with the working directory being the top level of the 
repository/Rstudio project).

Look at the code to see how I have used the different functions. Some
documentation is also available by doing `help(package = "afblue")`


