
# mapbayr-shiny

<!-- badges: start -->
<!-- badges: end -->

This is a collection of Shiny applications dedicated to model-based therapeutic drug monitoring, which are based on [mapbayr](https://github.com/FelicienLL/mapbayr) for the *maximum a posteriori* Bayesian estimation of pharmacokinetic parameters. 

Github let you the possibility to browse the code (see above), but not to see the applications in action, because a server is needed for running them. 

Fortunately, if R is installed on your computer, you can run them locally! 

First, make sure **shiny**, **mrgsolve** and **mapbayr** packages are installed on your computer.
Then, choose the application you want to run with `runGitHub()`. For instance: 
```r
shiny::runGitHub("mapbayr-shiny", "FelicienLL", subdir = "901-carboplatin")
```

You can also clone or download the repository on your machine.

### Warning

These application are provided **for illustrative purpose only**! Although they may be based on models built from clinical data, they should not be used to directly to individualize doses of patients. These applications are distributed **without any warantee**.