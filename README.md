# visualization_shiny

R Shiny App to visualize the expected cost

### Getting Started

-   Depedencies and environment: `renv`
    -   `renv` is a virtual environment for R. It is similar to `virtualenv` in Python.
    -   To install `rvenv`, run the following command: `install.packages("renv")`
    -   To install all the dependencies, run the following command: (it will install all the dependencies mentioned in the `renv.lock` file) `renv::restore()`
    -   (Optional) To create a new `renv` environment, run the following command: `renv::init()`
    -   (Optional) To activate the `renv` environment, run the following command: `renv::activate()`
    -   (Optional) To deactivate the `renv` environment, run the following command: `renv::deactivate()`
-   Start R shiny app
    -   To start the R shiny app, run the following command: `shiny::runApp()`
    -   The app will be available at `http://localhost:4367/`

### Parameters
- Detailed parameter and its description can be found in `parameters_lookup.csv` file

### Branch Logic
- The branch logic can be view in the `placeholder.png` file
