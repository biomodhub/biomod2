# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# biomod2 package object updating
# 13/06/18 - Damien G.
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

# Description =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# This file contains a set of function usefull 
# to update your biomod2 objects when you make
# a package updating and you want to keep working
# with simulation done with a previous version of 
# the package.
#
# The functions should be updated ecah time a slot
# is added to a biomod2 object.
#
# The main function is update.objects which takes
# a output of BIOMOD_FormatingData() (i.e 
# 'BIOMOD.formated.data' or 'BIOMOD.formated.data.PA')
# or a BIOMOD_Modeling() output. If you set recursive
# parameter to TRUE (default), all depending objects
# (e.g individual models) will be also updated.
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

# to do =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# Implement the same functions for BIOMOD_ModelingOptions(),
# BIOMOD_projection and BIOMOD_EnsembleModeling() outputs.
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

update_objects <- function(obj, recursive=TRUE){
  if(inherits(obj,'BIOMOD.formated.data') | inherits(obj,'BIOMOD.formated.data.PA')){
    cat("\n=-=-=- BIOMOD.formated.data conversion")
    update.objects_BIOMOD.formated.data(obj, recursive=recursive)
  }
  
  if(inherits(obj,'BIOMOD.models.out')){
    cat("\n=-=-=- BIOMOD.models.out conversion")
    update.objects_BIOMOD.models.out(obj, recursive=recursive)
  }
  
}

update.objects_BIOMOD.formated.data <- function(obj,recursive=TRUE){
  slots_to_add <- test_slots(obj)
  if(length(slots_to_add)){
    ref_obj <- new(class(obj))
    slot_to_add_str <- paste(slots_to_add,"=ref_obj@",slots_to_add, sep="", collapse=",")
    eval(parse(text= paste ( "obj <- new('",class(obj),"', obj,",slot_to_add_str,")",sep="")))
#     eval(parse(text= paste ( "current_obj <- new('",class(obj),"', obj,",slot_to_add_str,")",sep="")))
#     obj <- current_obj
  } else {
    cat("\tnoting to do!\n")
  }
  return(obj)
}

update.objects_biomod2_model <- function(obj,model.out=NULL){
  slots_to_add <- test_slots(obj)
  if(length(slots_to_add)){
    ref_obj <- new(class(obj))
    
    ## special args
    if(sum(c("expl_var_names","expl_var_type", "expl_var_range") %in% slots_to_add, na.rm=T) > 0){
      if(length(model.out)){
        data <- getModelsInputData(model.out,'expl.var')
        
        if("expl_var_names" %in% slots_to_add) ref_obj@expl_var_names <- names(data)
        if("expl_var_type" %in% slots_to_add) ref_obj@expl_var_type <- get_var_type(data)
        if("expl_var_range" %in% slots_to_add) ref_obj@expl_var_range <- get_var_range(data)
      }
    }

    slot_to_add_str <- paste(slots_to_add,"=ref_obj@",slots_to_add, sep="", collapse=",")
    eval(parse(text= paste ( "obj <- new('",class(obj),"', obj,",slot_to_add_str,")",sep="")))
#     eval(parse(text= paste ( "current_obj <- new('",class(obj),"', obj,",slot_to_add_str,")",sep="")))
#     obj <- current_obj
  } else {
    cat("\n\tnoting to do!\n")
  }
  return(obj)
}

update.objects_BIOMOD.models.out <- function(obj, recursive=T){
  slots_to_add <- test_slots(obj)
  if(length(slots_to_add)){
    ref_obj <- new("BIOMOD.models.out")
    slot_to_add_str <- paste(slots_to_add,"=ref_obj@",slots_to_add, sep="", collapse=",")
    eval(parse(text= paste ( "current_obj <- new('BIOMOD.models.out', obj,",slot_to_add_str,")",sep="")))
  } 
  
  if(recursive){
    cat("\n\tBIOMOD.formated.data checking...")
    data <- getModelsInputData(obj)
    data <- update.objects_BIOMOD.formated.data(data)
    save(data, file=obj@formated.input.data@link)
    
    cat("\n\tIndividual models checking...\n")
    mod_to_check <- BIOMOD_LoadModels(obj)  
    for(mtc in mod_to_check){
      cat("\t",mtc,"\n")
      assign(mtc,update.objects_biomod2_model(get(mtc),model.out=obj))
      save(list=mtc,file=file.path(obj@sp.name,"models",obj@modeling.id, mtc))
    }
  }
  
  obj_name <- tail(unlist(strsplit(obj@link,.Platform$file.sep,fixed=T)),1)
  assign(obj_name, obj)
  save( list=obj_name,file=obj@link)
  
}

test_slots <- function(obj){
  out <- NULL
  for(slot in slotNames(obj)){
    test <- try(getElement(obj,slot),silent=T)
    if(inherits(test,"try-error")){
      cat("\t\tredefining",slot,"slot\n")
      out <- c(out, slot)
    }
  }
  return(out)
}
