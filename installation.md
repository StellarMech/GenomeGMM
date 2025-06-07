# The package depends on R(>= 3.5.0)
--------------
## 1.Introduction
This package provides a function to estimate polyploid genome size using a Gaussian Mixture Model (GMM) based on k-mer frequency histograms.
## 2. Installation
- First of all, you need install dependencies required by  GenomeGMM under R, if you do not have them yet. Open a **R terminal**,
```R
install.package("devtools")
```
- Secondly, you have two ways to intsall  GenomeGMM
  - way 1 : install with a local copy
    - 1.1 download a zipped  GenomeGMM and create a tar.gz
      ```linux
      cd path/to/target/folder
      wget 
      unzip master.zip
      tar -czvf 
        ```
    - 1.2 install under R
      ```R
       install.packages("")
       q("no")
      ```
  - way 2 : directly install from github
    ```R
     install.packages("devtools")
     devtools::install_github("")
    ```
