---
name: Help Wanted
about: " How to use biomod2 functions and outputs"
title: Help with BIOMOD_xxx - [short question here]
labels: help wanted
assignees: ''

---

_Please make sure to close the issue once you consider it as solved_
_Please use screenshots only when you cannot copy-paste the object, e.g. for figures or maps_

**Context and question**
Please describe the context of the question and ask the question here.

**Code used**
Please add the code used here to illustrate your question. Start with `BIOMOD_FormatingData` up to the function on which you have questions. Please add as well the output of `show` for the different object used or generated. 

```
# Example for a question about BIOMOD_EnsembleForecasting argument or outputs
# make sure to add other relevant functions (e.g. bm_PseudoAbsences, bm_CrossValidation, get_evaluations ...)
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
# if question is on argument, stop here
show(myBiomodEMProj)
# paste output here
``` 

**Other outputs** _(optional)_
If you have questions about a figure or maps, you can add a screenshot here.

**Environment Information**
Please paste the output of `sessionInfo()` in your current R session below.
```
# paste output of sessionInfo() here
```
