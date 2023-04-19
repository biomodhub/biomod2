---
name: Bug report
about: Unexpected error in function or output
title: Error in BIOMOD_xxx - [short error summary here]
labels: bug
assignees: ''

---

_Please make sure to close the issue once you consider it as solved_
_Please use screenshots only when you cannot copy-paste the object, e.g. for figures or maps_

**Error and context**
A short description of the error and its context.

**Code used to get the error**
Please add the code used to reproduce your error, starting with `BIOMOD_FormatingData` up to the function that bugged. Please add as well the output of `show` for the different object used or generated. 

```
# Example for an error in BIOMOD_EnsembleForecasting
# make sure to add other relevant functions (e.g. bm_PseudoAbsences or bm_CrossValidation)
# make sure to remove functions you did not run

myBiomodData <- BIOMOD_FormatingData( ** write arguments here **)
show(myBiomodData)
# paste output here
show(myExpl) # if using an environment raster
# paste output here

myBiomodModelOut <- BIOMOD_Modeling( ** write arguments here **)
show(myBiomodModelOut)
# paste output here

myBiomodEM <- BIOMOD_EnsembleModeling( ** write arguments here **)
show(myBiomodEM)
# paste output here

myBiomodProj <- BIOMOD_Projection( ** write arguments here **)
show(myBiomodProj)
# paste output here
show(myExplFuture) # if projecting on a new environment raster
# paste output here

myBiomodEMProj <- BIOMOD_EnsembleForecasting( ** write arguments here **)
show(myBiomodEMProj)
# paste output here
show(myExplFuture) # if projecting on a new environment raster
# paste output here
```

**Environment Information**
Please paste the output of `sessionInfo()` in your current R session below.
```
# paste output of sessionInfo() here
```

**Additional information**
If you have any additional information or context you can add it here.
