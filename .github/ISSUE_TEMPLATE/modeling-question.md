---
name: Modeling Question
about: About modeling workflow and output
title: Question about [generic topic] - [short and concise question]
labels: modeling question
assignees: ''

---

_Please make sure to close the issue once you consider it as solved_
_Please use screenshots only when you cannot copy-paste the object, e.g. for figures or maps_

**Context and question**
Ask your question here and add the relevant context.

**Related code** _(optional)_
If your question is linked to some function outputs, please add the code used, starting with `BIOMOD_FormatingData` up to the output in question. Please add as well the output of `show` for the different object used or generated. 

```
# Example for a question on BIOMOD_EnsembleModeling output
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

# for a question on evaluation of ensemble models:
get_evaluations(myBiomodEM)
# paste output here
```

**Other outputs** _(optional)_
If your question is about a figure or a map, you can add a screenshot here.

**Environment Information** _(optional)_
If the question relates to model output, please paste the output of `sessionInfo()` in your current R session below.
```
# paste output of sessionInfo() here
```
