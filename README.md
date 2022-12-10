# Linkage Analysis in R
<img src="https://camo.githubusercontent.com/f0c75c42d670c1153720d72688ab576936104b7b9a80ea336eeba472949394e6/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f522533452533442d332e322e342d3636363666662e737667" style="width:105px"> <img src="https://camo.githubusercontent.com/628aedf920d3cee6c4467d1f63915b015f131749861bcec4178bdeca3cf3810b/68747470733a2f2f7777772e722d706b672e6f72672f6261646765732f76657273696f6e2f75736574686973" style="width:115px"> <img src="https://camo.githubusercontent.com/7d329492423dc3eceea8e9b170b73d4f5e8d4b0c4878376bb0f94c7c91c2608e/68747470733a2f2f73766773686172652e636f6d2f692f5a68792e737667" style="width:80px"> <img src="https://img.shields.io/github/license/Naereen/StrapDown.js.svg" style="width:110px">

## Friendly Reminder

If you use or take inspiration from this repository please cite with this link: [santurini/DAG-Linkage-Analysis-in-R](https://github.com/santurini/DAG-Linkage-Analysis-in-R)

Your support will be truly appreciated and feel free to contact me at my following links or just send me an email:
- [Linkedin](https://www.linkedin.com/in/arturo-ghinassi-50b8a0219/)
- [Kaggle](https://www.kaggle.com/santurini)
- ghinassi.1863151@studenti.uniroma1.it

## Overview

In this work we implemented to universal hypothesis test to check whether a link or a path existed or not in a directed acyclic graph.
First we checked the properties of the tests and then we applied both on a real case study.

## Content

**scripts** this folder contains the scripts:
  - _functions.R_ which contains all the functions to run the scripts
  - _LinkageTest.R_ which contains the implementation of the test for a single link
  - _PathwayTest.R_ which contains the implementation of the test for a directed pathway
  - _Cytometry.R_ which contains the application of the tests on real data
  - _plotSparsity.R_ which contains the script fro the scatter plots of the test results with different sparsity values

**report** this folder contains the final report of the case study:
  - _Report.html_ is the knitted html file written in rmarkdown
  
**data** is a folder which contains:
  - _cytometry-data.xlsx_ the excel file with 9 different sheets

## Theory concepts

The aim of this repository is to implement the Universal Hypothesis Test that is explained in the following image:
<p align="center">
  <img src="https://user-images.githubusercontent.com/91251307/152517117-4ec8e32c-8dc5-498f-989d-416259ab667b.png" style="width:800px">
</p>

### Test of graph linkages

This is the definition we used to implement the graph linkage test:
<p align="center">
  <img src="https://user-images.githubusercontent.com/91251307/152516437-8c2ba0df-181f-49cd-afc2-08fad0ea7e5a.png" style="width:800px">
</p>


### Test of directed pathway
This is the definition we used to implement the graph pathway test:
<p align="center">
  <img src="https://user-images.githubusercontent.com/91251307/152516816-2b99b5ae-b8b6-4a54-9f59-3ce30df069c2.png" style="width:800px">
</p>

## Tools

To be able to implement it we used the **clrdag** package which contains the *MLEdag* function, the linck to the repository is the following:
  - [chunlinli/clrdag](https://github.com/chunlinli/clrdag)

## Data

Once implemented we tested them on real world [data](https://docs.google.com/spreadsheets/d/11Xv1YJHUshjz2kfA-11XqVzsQaRxTfSmoEMyM5rXypM/edit#gid=357906815) regarding [Cell Signaling](https://en.wikipedia.org/wiki/Flow_cytometry) represented as a DAG:
<p align="center">
  <img src="https://user-images.githubusercontent.com/91251307/152517293-f0d880e3-2e92-4015-8d10-52692a02a83e.png" style="width:500px" title="Cell Signaling">
  <img src="https://user-images.githubusercontent.com/91251307/152517581-0ade8d8d-0ede-428c-b6fc-d80dca78613c.png" style="width:500px" title="Protein DAG">
</p>


## Overview

  1) Estimating the variance and log-likelihood: 
      - Sigma estimate:
        <p align="center">
          <img src="https://user-images.githubusercontent.com/91251307/152529071-638701d9-d995-44da-a838-0025b74fe1e7.png" style="width:400px">
        </p>
      - log-Likelihood estimate:
        <p align="center">
          <img src="https://user-images.githubusercontent.com/91251307/152527515-798d5ae0-9431-430c-8d1a-80b25d1f347b.png" style="width:600px">
        </p>
        
  2) Implementing the LRT (Likelihood ratio test) and Crossfit test
  3) Implementing the Universal tests on random data to compute size and power
     - Formula to generate the random data:
        <p align="center">
          <img src="https://user-images.githubusercontent.com/91251307/152529840-011ffda1-abc6-43a5-8e48-3352969419a4.png" style="width:600px">
        </p>
  4) Applying the test on real world data to check linkages between proteins
  
  <br />
<br />
<p align="center">
    <img src="https://user-images.githubusercontent.com/50860347/147412786-183da6b0-990f-4016-9f2e-0719d8066f5b.png" style="width: 100%"/>
<p>

<br />
