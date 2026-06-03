library(urca)
#
# ADF-test
#
urcadf <- function(){
  if (!checkActiveDataSet()) return()
  if (!checkNumeric(n=1)) return()
  initializeDialog(title="ADF Test")
  xBox <- variableListBox(top, Numeric(), title="Variable (pick one)")
  onOK <- function(){
    x <- getSelection(xBox)
    var <- paste(.activeDataSet, "$", x, sep="")
    ttype <- tclvalue(testtypeVariable)
    lags <- tclvalue(lagsVariable)
    if (length(x) == 0){
      errorCondition(recall=urcadf, message="You must select a variable.")
      return()
      }
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    doItAndPrint(paste("summary(ur.df(y= ", var, ", type = ", ttype, ", lags = ",  lags, "))", sep=""))
    tkdestroy(top)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="ur.df")
  radioButtons(name="testtype", buttons=c("none", "drift", "trend"), values=c("'none'", "'drift'", "'trend'" ), labels=c("no deterministic regressor", "drift only", "drift and trend"), title="Type of test")
  rightFrame <- tkframe(top)
  lagsFrame <- tkframe(rightFrame)
  lagsVariable <- tclVar("4")
  lagsField <- tkentry(lagsFrame, width="2", textvariable=lagsVariable)
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(tklabel(rightFrame, text=""))
  tkgrid(tklabel(lagsFrame, text="Maximum number of lags = ", fg="blue"), lagsField, sticky="w")
  tkgrid(lagsFrame, sticky="w")
  tkgrid(testtypeFrame, rightFrame, sticky="nw")
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  tkgrid.configure(lagsField, sticky="e")
  dialogSuffix(rows=4, columns=2)
}
#
# ERS-test
#
urcaers <- function(){
  if (!checkActiveDataSet()) return()
  if (!checkNumeric(n=1)) return()
  initializeDialog(title="ERS Test")
  xBox <- variableListBox(top, Numeric(), title="Variable (pick one)")
  onOK <- function(){
    x <- getSelection(xBox)
    var <- paste(.activeDataSet, "$", x, sep="")
    ttype <- tclvalue(testtypeVariable)
    modtype <- tclvalue(modelVariable)
    lags <- tclvalue(lagsVariable)
    if (length(x) == 0){
      errorCondition(recall=urcaers, message="You must select a variable.")
      return()
      }
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    doItAndPrint(paste("summary(ur.ers(y= ", var, ", type = ", ttype, ", model = ", modtype, ", lag.max = ",  lags, "))", sep=""))
    tkdestroy(top)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="ur.ers")
  radioButtons(name="testtype", buttons=c("DFGLS", "Ptest"), values=c("'DF-GLS'", "'P-test'"), 
               labels=c("DF-GLS statistic", "P-test statistic"), title="Type of test")
  radioButtons(name="model", buttons=c("const", "trend"), values=c("'constant'", "'trend'"), 
               labels=c("Include constant", "Include constant + trend"), title="Model type")
  rightFrame <- tkframe(top)
  lagsFrame <- tkframe(rightFrame)
  lagsVariable <- tclVar("4")
  lagsField <- tkentry(lagsFrame, width="2", textvariable=lagsVariable)
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(tklabel(rightFrame, text=""))
  tkgrid(tklabel(lagsFrame, text="Maximum number of lags = ", fg="blue"), lagsField, sticky="w")
  tkgrid(lagsFrame, sticky="w")
  tkgrid(testtypeFrame, rightFrame, sticky="nw")
  tkgrid(modelFrame, rightFrame, sticky="nw")
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  tkgrid.configure(lagsField, sticky="e")
  dialogSuffix(rows=4, columns=2)
}
#
# KPSS-test
#
urcakpss <- function(){
  if (!checkActiveDataSet()) return()
  if (!checkNumeric(n=1)) return()
  initializeDialog(title="KPSS Test")
  xBox <- variableListBox(top, Numeric(), title="Variable (pick one)")
  radioButtons(name="testtype", buttons=c("constant", "trend"), 
        values=c("'mu'", "'tau'"), initialValue="'mu'", 
        labels=c("Include constant", "Include trend"), title="Deterministic Part")
  radioButtons(name="lags", buttons=c("short", "long", "nil"), 
        values=c("'short'", "'long'", "'nil'"), initialValue="'short'",         labels=c("use short lags", "use long lags", "use no lags"), title="Lag Selection")
  onOK <- function(){
    x <- getSelection(xBox)
    if (length(x) == 0){
      errorCondition(recall=urcakpss, message="You must select a variable.")
      return()
    }
    var <- paste(.activeDataSet, "$", x, sep="")
    ttype <- tclvalue(testtypeVariable)
    ltype <- tclvalue(lagsVariable)
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    command <- paste("summary(ur.kpss(", var, ", type = ", ttype, ", lags = ", ltype, "))", sep="")
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="ur.kpss")
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(testtypeFrame, lagsFrame,  sticky="w")
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  dialogSuffix(rows=3, columns=2)
}
#
# Schmidt-Phillips test
#
urcasp <- function(){
  if (!checkActiveDataSet()) return()
  if (!checkNumeric(n=1)) return()
  initializeDialog(title="SP Test")
  xBox <- variableListBox(top, Numeric(), title="Variable (pick one)")
  onOK <- function(){
    x <- getSelection(xBox)
    var <- paste(.activeDataSet, "$", x, sep="")
    ttype <- tclvalue(testtypeVariable)
    pdtype <- tclvalue(poldegVariable)
    sltype <- tclvalue(signifVariable)
    if (length(x) == 0){
      errorCondition(recall=urcasp, message="You must select a variable.")
      return()
      }
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    doItAndPrint(paste("summary(ur.sp(y= ", var, ", type = ", ttype, ", pol.deg = ", pdtype, ", signif = ",  sltype, "))", sep=""))
    tkdestroy(top)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="ur.sp")
  radioButtons(name="testtype", buttons=c("tau", "rho"), values=c("'tau'", "'rho'"), 
               labels=c("tau statistic", "rho statistic"), title="Type of test")
  radioButtons(name="poldeg", buttons=c("pb1", "pb2", "pb3", "pb4"), values=c("1", "2", "3", "4"), 
               labels=c("1st", "2nd", "3rd", "4th"), title="Select pol. degree")
  radioButtons(name="signif", buttons=c("sb1", "sb5", "sb10"), values=c("0.01", "0.05", "0.1"), 
               labels=c("alpha=1%", "alpha=5%", "alpha=10%"), title="Sig. Level")
  tkgrid(getFrame(xBox), testtypeFrame, sticky="nw")
  tkgrid(poldegFrame, signifFrame, sticky="nw")
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  dialogSuffix(rows=4, columns=2)
}
#
# Phillips-Perron test
#
urcapp <- function(){
  if (!checkActiveDataSet()) return()
  if (!checkNumeric(n=1)) return()
  initializeDialog(title="PP Test")
  xBox <- variableListBox(top, Numeric(), title="Variable (pick one)")
  onOK <- function(){
    x <- getSelection(xBox)
    var <- paste(.activeDataSet, "$", x, sep="")
    ttype <- tclvalue(testtypeVariable)
    modtype <- tclvalue(modelVariable)
    lagtype <- tclvalue(lagsVariable)
    if (length(x) == 0){
      errorCondition(recall=urcapp, message="You must select a variable.")
      return()
      }
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    doItAndPrint(paste("summary(ur.pp(x= ", var, ", type = ", ttype, ", model = ", modtype, ", lags = ",  lagtype, "))", sep=""))
    tkdestroy(top)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="ur.pp")
  radioButtons(name="testtype", buttons=c("Zalpha", "Ztau"), values=c("'Z-alpha'", "'Z-tau'"), 
               labels=c("Z-alpha statistic", "Z-tau statistic"), title="Type of test")
  radioButtons(name="model", buttons=c("const", "trend"), values=c("'constant'", "'trend'"), 
               labels=c("Include constant", "include constant + trend"), title="Model type")
  radioButtons(name="lags", buttons=c("short", "long"), values=c("'short'", "'long'"), 
               labels=c("short lags", "long lags"), title="Lags for error correction")
  tkgrid(getFrame(xBox), testtypeFrame, sticky="nw")
  tkgrid(modelFrame, lagsFrame, sticky="nw")
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  dialogSuffix(rows=4, columns=2)
}
#
# Zivot-Andrews test
#
urcaza <- function(){
  if (!checkActiveDataSet()) return()
  if (!checkNumeric(n=1)) return()
  initializeDialog(title="Zivot & Andrews Test")
  xBox <- variableListBox(top, Numeric(), title="Variable (pick one)")
  onOK <- function(){
    x <- getSelection(xBox)
    var <- paste(.activeDataSet, "$", x, sep="")
    modtype <- tclvalue(modelVariable)
    lags <- tclvalue(lagsVariable)
    ptype <- if(tclvalue(plotVariable) == 0) "FALSE" else "TRUE"
    if (length(x) == 0){
      errorCondition(recall=urcaza, message="You must select a variable.")
      return()
      }
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    command <- paste("ur.za(y= ", var, ",  model = ", modtype, ", lag = ",  lags, ")", sep="")
    logger(paste("ZAstat <- ", command, sep=""))
    assign("ZAstat", justDoIt(command), envir=.GlobalEnv)
    doItAndPrint("summary(ZAstat)")
    if(ptype==TRUE) {
      justDoIt("x11()")
      justDoIt("plot(ZAstat)")
    }
    logger("remove(ZAstat)") 
    remove(ZAstat, envir=.GlobalEnv)       
    tkdestroy(top)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="ur.za")
  radioButtons(name="model", buttons=c("const", "trend", "both"), values=c("'intercept'", "'trend'", "'both'"), 
               labels=c("Include constant", "Include trend", "Include both"), title="Model type")
  checkBoxes(frame="plotFrame", boxes="plot", initialValues="0", labels="Plot path of Zivot & Andrews Statistic?")
  rightFrame <- tkframe(top)
  lagsFrame <- tkframe(rightFrame)
  lagsVariable <- tclVar("4")
  lagsField <- tkentry(lagsFrame, width="2", textvariable=lagsVariable)
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(tklabel(rightFrame, text=""))
  tkgrid(tklabel(lagsFrame, text="Maximum number of lags = ", fg="blue"), lagsField, sticky="w")
  tkgrid(lagsFrame, sticky="w")
  tkgrid(modelFrame, rightFrame, sticky="nw")
  tkgrid(plotFrame, sticky="nw")
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  tkgrid.configure(lagsField, sticky="e")
  dialogSuffix(rows=3, columns=2)
}
#
# Phillips-Ouliaris test
#
urcacapo <- function(){
  dataSets <- listDataSets()
  if (length(dataSets) == 0){
    tkmessageBox(message="There are no data sets from which to choose.", icon="error", type="ok")
    tkfocus(CommanderWindow())
    return()
  }
  initializeDialog(title="Phillips & Ouliaris Test")
  xBox <- variableListBox(top, dataSets, title="Data Sets (pick one)")
  onOK <- function(){
    x <- getSelection(xBox)
    if (length(x) == 0){
      errorCondition(recall=urcacapo, message="You must select a data set.")
      return()
    }
    meantype <- tclvalue(demeanVariable)
    lags <- tclvalue(lagsVariable)
    ttype <- tclvalue(typeVariable)
    tolvar <- tclvalue(tolVariable)
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    doItAndPrint(paste("summary(ca.po(z = ", x, ",  demean = ", meantype, ", lag = ",  lags, ", type = ", ttype, ", tol = ", tolvar, "))", sep=""))
    tkdestroy(top)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="ca.po")
  radioButtons(name="demean", buttons=c("none", "constant", "trend"), values=c("'none'", "'constant'", "'trend'"), 
               labels=c("None", "Include constant", "Include trend"), title="Demean?")
  radioButtons(name="lags", buttons=c("short", "long"), values=c("'short'", "'long'"), 
               labels=c("short lags", "long lags"), title="Lags for error correction")
  radioButtons(name="type", buttons=c("Pu", "Pz"), values=c("'Pu'", "'Pz'"), 
               labels=c("Pu statistic", "Pz statistic"), title="Type of test")
  rightFrame <- tkframe(top)
  tolFrame <- tkframe(rightFrame)
  tolVariable <- tclVar("NULL")
  tolField <- tkentry(tolFrame, width="8", textvariable=tolVariable)
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(tklabel(rightFrame, text=""))
  tkgrid(tklabel(tolFrame, text="Tolerance level =  ", fg="blue"), tolField, sticky="w")
  tkgrid(tolFrame, sticky="w")
  tkgrid(demeanFrame, rightFrame, sticky="nw")
  tkgrid(lagsFrame, typeFrame, sticky="nw")
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  tkgrid.configure(tolField, sticky="e")
  dialogSuffix(rows=5, columns=2)
}
#
# Johansen-Procedure
#
urcacajo <- function(){
  dataSets <- listDataSets()
  if (length(dataSets) == 0){
    tkmessageBox(message="There are no data sets from which to choose.", icon="error", type="ok")
    tkfocus(CommanderWindow())
    return()
  }
  initializeDialog(title="Johansen's Procedure")
  xBox <- variableListBox(top, dataSets, title="Data Sets (pick one)")
  assign("UpdateModelNumber", UpdateModelNumber() + 1, envir=.GlobalEnv)
  modelName <- tclVar(paste("VECMmodel.", UpdateModelNumber, sep=""))
  modelFrame <- tkframe(top)
  model <- tkentry(modelFrame, width="20", textvariable=modelName)
  onOK <- function(){
    x <- getSelection(xBox)
    if (length(x) == 0){
      errorCondition(recall=urcacajo, message="You must select a data set.")
      return()
    }
    ttype <- tclvalue(typeVariable)
    spect <- tclvalue(specVariable)
    ctype <- tclvalue(ecdetVariable)
    lags <- tclvalue(lagVariable)
    seas <- tclvalue(seasonVariable)
    dummy <- tclvalue(dumVariable)
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    modelValue <- tclvalue(modelName)
    if (!is.valid.name(modelValue)){
      assign("UpdateModelNumber", UpdateModelNumber() - 1, envir=.GlobalEnv)
      errorCondition(recall=urcacajo, message=paste('"', modelValue, '" is not a valid name.', sep=""))
      return()
    }
    if (is.element(modelValue, listVECMmodels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type="Model"))){
        assign("UpdateModelNumber", UpdateModelNumber() - 1, envir=.GlobalEnv)
        if (GrabFocus()) tkgrab.release(top)
        tkdestroy(top)
        urcacajo()
        return()
      }
    }
    command <- paste("ca.jo(x = ", x, ",  type = ", ttype, ", ecdet = ",  ctype, ", K = ",
                     lags, ", spec = ", spect, ", season = ", seas, ", dumvar = ", dummy , ")", sep="")
    logger(paste(modelValue, " <- ", command, sep=""))
    assign(modelValue, justDoIt(command), envir=.GlobalEnv)
    doItAndPrint(paste("summary(", modelValue, ")", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="ca.jo")
  radioButtons(name="type", buttons=c("eigen", "trace"), values=c("'eigen'", "'trace'"), 
               labels=c("Eigenvalue statistic", "Trace statistic"), title="Type of statistic")
  radioButtons(name="spec", buttons=c("long", "trans"), values=c("'longrun'", "'transitory'"), 
               labels=c("longrun specification", "transitory specification"), title="VECM specification")
  radioButtons(name="season", buttons=c("none", "monthly", "quarterly"), values=c("NULL", "12", "4"), 
               labels=c("None", "Monthly seasonality", "Quarterly seasonality"), title="Seasonality")
  radioButtons(name="ecdet", buttons=c("none", "const", "trend"), values=c("'none'", "'const'", "'trend'"), 
               labels=c("none", "constant", "trend"), title="Deterministic Variable in Cointegration") 
  rightFrame <- tkframe(top)
  lagFrame <- tkframe(rightFrame)
  lagVariable <- tclVar("2")
  lagField <- tkentry(lagFrame, width="2", textvariable=lagVariable)
  dumFrame <- tkframe(rightFrame)
  dumVariable <- tclVar("NULL")
  dumField <- tkentry(dumFrame, width="8", textvariable=dumVariable)
  tkgrid(tklabel(modelFrame, text="Enter name for model:"), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(tklabel(rightFrame, text=""), sticky="w")
  tkgrid(tklabel(lagFrame, text="Lag order =  ", fg="blue"), lagField, sticky="w")
  tkgrid(lagFrame, sticky="w")
  tkgrid(tklabel(dumFrame, text="Matrix of dummy variables =  ", fg="blue"), dumField, sticky="w")
  tkgrid(dumFrame, sticky="w")
  tkgrid(typeFrame, rightFrame, sticky="nw")
  tkgrid(specFrame, seasonFrame, sticky="nw")
  tkgrid(ecdetFrame, rightFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  tkgrid.configure(lagField, sticky="e")
  tkgrid.configure(dumField, sticky="e")
  dialogSuffix(rows=8, columns=2)
}
#
# Linear Trend test
#
urcalttest <- function(){
  VECMmodels <- listVECMmodels()
  if (length(VECMmodels) == 0){
    tkmessageBox(message="There are no VECM models from which to choose.", icon="error", type="ok")
    tkfocus(CommanderWindow())
    return()
  }
  initializeDialog(title="Linear Trend Test")
  xBox <- variableListBox(top, VECMmodels, title="VECM models (pick one)")
  onOK <- function(){
    x <- getSelection(xBox)
    if (length(x) == 0){
      errorCondition(recall=urcalttest, message="You must select a VECM model.")
      return()
    }
    rint <- tclvalue(rankVariable)
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    doItAndPrint(paste("lttest(z = ", x, ", r = ", rint, ")", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="lttest")
  rankFrame <- tkframe(top)
  rankVariable <- tclVar("1")
  rankField <- tkentry(rankFrame, width="2", textvariable=rankVariable)
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(rankFrame, sticky="w")
  tkgrid(tklabel(rankFrame, text="Number of cointegrating relationships =  ", fg="blue"), rankField, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  tkgrid.configure(rankField, sticky="e")
  dialogSuffix(rows=3, columns=2)
}
#
# Restrictions on Loading matrix
#
urcaalrtest <- function(){
  VECMmodels <- listVECMmodels()
  matrices <- listMatrix()
  if (length(VECMmodels) == 0){
    tkmessageBox(message="There are no VECM models from which to choose.", icon="error", type="ok")
    tkfocus(CommanderWindow())
    return()
  }
  if (length(matrices) == 0){
    tkmessageBox(message="There are no restriction matrices defined from which to choose.", icon="error", type="ok")
    tkfocus(CommanderWindow())
    return()
  }
  initializeDialog(title="Test Restrictions on Loading Vectors")
  xBox <- variableListBox(top, VECMmodels, title="VECM models (pick one)")
  yBox <- variableListBox(top, matrices, title="Restriction matrices (pick one)")
  onOK <- function(){
    x <- getSelection(xBox)
    if (length(x) == 0){
      errorCondition(recall=urcaalrtest, message="You must select a VECM model.")
      return()
    }
    y <- getSelection(yBox)
    if (length(y) == 0){
      errorCondition(recall=urcaalrtest, message="You must select a restriction matrix.")
      return()
    }
    rint <- tclvalue(rankVariable)
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    doItAndPrint(paste("summary(alrtest(z = ", x, ", A = ", y , ", r = ", rint, "))", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="alrtest")
  rankFrame <- tkframe(top)
  rankVariable <- tclVar("1")
  rankField <- tkentry(rankFrame, width="2", textvariable=rankVariable)
  tkgrid(getFrame(xBox), getFrame(yBox), sticky="nw")
  tkgrid(rankFrame, sticky="w")
  tkgrid(tklabel(rankFrame, text="Number of cointegrating relationships =  ", fg="blue"), rankField, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  tkgrid.configure(rankField, sticky="e")
  dialogSuffix(rows=3, columns=1)
}
#
# Restrictions on cointegration matrix
#
urcablrtest <- function(){
  VECMmodels <- listVECMmodels()
  matrices <- listMatrix()
  if (length(VECMmodels) == 0){
    tkmessageBox(message="There are no VECM models from which to choose.", icon="error", type="ok")
    tkfocus(CommanderWindow())
    return()
  }
  if (length(matrices) == 0){
    tkmessageBox(message="There are no restriction matrices defined from which to choose.", icon="error", type="ok")
    tkfocus(CommanderWindow())
    return()
  }
  initializeDialog(title="Test Restrictions on Cointegration Vectors")
  xBox <- variableListBox(top, VECMmodels, title="VECM models (pick one)")
  yBox <- variableListBox(top, matrices, title="Restriction matrices (pick one)")
  onOK <- function(){
    x <- getSelection(xBox)
    if (length(x) == 0){
      errorCondition(recall=urcaablrtest, message="You must select a VECM model.")
      return()
    }
    y <- getSelection(yBox)
    if (length(y) == 0){
      errorCondition(recall=urcaalrtest, message="You must select a restriction matrix.")
      return()
    }
    rint <- tclvalue(rankVariable)
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    doItAndPrint(paste("summary(blrtest(z = ", x, ", H = ", y , ", r = ", rint, "))", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="blrtest")
  rankFrame <- tkframe(top)
  rankVariable <- tclVar("1")
  rankField <- tkentry(rankFrame, width="2", textvariable=rankVariable)
  tkgrid(getFrame(xBox), getFrame(yBox), sticky="nw")
  tkgrid(rankFrame, sticky="w")
  tkgrid(tklabel(rankFrame, text="Number of cointegrating relationships =  ", fg="blue"), rankField, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  tkgrid.configure(rankField, sticky="e")
  dialogSuffix(rows=3, columns=1)
}
#
# Restrictions for partly known cointegrating vectors
#
urcabh5lrtest <- function(){
  VECMmodels <- listVECMmodels()
  matrices <- listMatrix()
  numerics <- listNumeric()
  elements <- c(matrices, numerics)
  if (length(VECMmodels) == 0){
    tkmessageBox(message="There are no VECM models from which to choose.", icon="error", type="ok")
    tkfocus(CommanderWindow())
    return()
  }
  if (length(matrices) == 0){
    tkmessageBox(message="There are no restriction matrices defined from which to choose.", icon="error", type="ok")
    tkfocus(CommanderWindow())
    return()
  }
  initializeDialog(title="Test validity of partly known cointegrating Vectors")
  xBox <- variableListBox(top, VECMmodels, title="VECM models (pick one)")
  yBox <- variableListBox(top, elements, title="Restriction matrices (pick one)")
  onOK <- function(){
    x <- getSelection(xBox)
    if (length(x) == 0){
      errorCondition(recall=urcabh5lrtest, message="You must select a VECM model.")
      return()
    }
    y <- getSelection(yBox)
    if (length(y) == 0){
      errorCondition(recall=urcabh5lrtest, message="You must select a restriction matrix.")
      return()
    }
    rint <- tclvalue(rankVariable)
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    doItAndPrint(paste("summary(bh5lrtest(z = ", x, ", H = ", y , ", r = ", rint, "))", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="bh5lrtest")
  rankFrame <- tkframe(top)
  rankVariable <- tclVar("2")
  rankField <- tkentry(rankFrame, width="2", textvariable=rankVariable)
  tkgrid(getFrame(xBox), getFrame(yBox), sticky="nw")
  tkgrid(rankFrame, sticky="w")
  tkgrid(tklabel(rankFrame, text="Number of cointegrating relationships =  ", fg="blue"), rankField, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  tkgrid.configure(rankField, sticky="e")
  dialogSuffix(rows=3, columns=1)
}
#
# Restrictions for partly known cointegrating vectors
#
urcabh6lrtest <- function(){
  VECMmodels <- listVECMmodels()
  matrices <- listMatrix()
  numerics <- listNumeric()
  elements <- c(matrices, numerics)
  if (length(VECMmodels) == 0){
    tkmessageBox(message="There are no VECM models from which to choose.", icon="error", type="ok")
    tkfocus(CommanderWindow())
    return()
  }
  if (length(matrices) == 0){
    tkmessageBox(message="There are no restriction matrices defined from which to choose.", icon="error", type="ok")
    tkfocus(CommanderWindow())
    return()
  }
  initializeDialog(title="Test restrictions of partly known cointegrating Vectors")
  xBox <- variableListBox(top, VECMmodels, title="VECM models (pick one)")
  yBox <- variableListBox(top, elements, title="Restriction matrices (pick one)")
  onOK <- function(){
    x <- getSelection(xBox)
    if (length(x) == 0){
      errorCondition(recall=urcabh6lrtest, message="You must select a VECM model.")
      return()
    }
    y <- getSelection(yBox)
    if (length(y) == 0){
      errorCondition(recall=urcaalrtest, message="You must select a restriction matrix.")
      return()
    }
    rint <- tclvalue(rankVariable)
    r1int <- tclvalue(r1Variable)
    maxiter <- tclvalue(maxiterVariable)
    conv <- tclvalue(conVariable)
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    doItAndPrint(paste("summary(bh6lrtest(z = ", x, ", H = ", y , ", r = ", rint, ", r1 =", r1int, ", conv.val =", conv, ", max.iter =", maxiter, "))", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="bh6lrtest")
  rankFrame <- tkframe(top)
  r1Frame <- tkframe(top)
  rankVariable <- tclVar("2")
  r1Variable <- tclVar("1")
  rightFrame <- tkframe(top)
  maxiterFrame <- tkframe(rightFrame)
  maxiterVariable <- tclVar("50")
  maxiterField <- tkentry(maxiterFrame, width="4", textvariable=maxiterVariable)
  conFrame <- tkframe(rightFrame)
  conVariable <- tclVar("0.0001")
  conField <- tkentry(conFrame, width="8", textvariable=conVariable)
  tkgrid(tklabel(rightFrame, text=""))
  rankField <- tkentry(rankFrame, width="2", textvariable=rankVariable)
  r1Field <- tkentry(r1Frame, width="2", textvariable=r1Variable)
  tkgrid(getFrame(xBox), getFrame(yBox), sticky="nw")
  tkgrid(rankFrame, rightFrame, sticky="w")
  tkgrid(tklabel(rankFrame, text="Number of cointegrating relationships =  ", fg="blue"), rankField, sticky="w")
  tkgrid(r1Frame, rightFrame, sticky="w")
  tkgrid(tklabel(r1Frame, text="Number of restricted ci relationships =  ", fg="blue"), r1Field, sticky="w")
  tkgrid(maxiterFrame, sticky="w")
  tkgrid(tklabel(maxiterFrame, text="Maximum number of iterations = ", fg="blue"), maxiterField, sticky="w")
  tkgrid(conFrame, sticky="w")
  tkgrid(tklabel(conFrame, text="Convergence criteria = ", fg="blue"), conField, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  tkgrid.configure(rankField, sticky="e")
  tkgrid.configure(r1Field, sticky="e")
  tkgrid.configure(maxiterField, sticky="e")
  tkgrid.configure(conField, sticky="e")
  dialogSuffix(rows=4, columns=2)
}
#
# Restrictions for loading and cointegration matrix
#
urcaablrtest <- function(){
  VECMmodels <- listVECMmodels()
  matrices <- listMatrix()
  if (length(VECMmodels) == 0){
    tkmessageBox(message="There are no VECM models from which to choose.", icon="error", type="ok")
    tkfocus(CommanderWindow())
    return()
  }
  if (length(matrices) == 0){
    tkmessageBox(message="There are no restriction matrices defined from which to choose.", icon="error", type="ok")
    tkfocus(CommanderWindow())
    return()
  }
  initializeDialog(title="Test Restrictions on Loading & Cointegration vectors")
  xBox <- variableListBox(top, VECMmodels, title="VECM models (pick one)")
  yBox <- variableListBox(top, matrices, title="Restriction matrices for cointegration vectors (pick one)")
  lBox <- variableListBox(top, matrices, title="Restriction matrices for loading vectors (pick one)")
  onOK <- function(){
    x <- getSelection(xBox)
    if (length(x) == 0){
      errorCondition(recall=urcaablrtest, message="You must select a VECM model.")
      return()
    }
    y <- getSelection(yBox)
    if (length(y) == 0){
      errorCondition(recall=urcaablrtest, message="You must select a cointegration restriction matrix.")
      return()
    }
    l <- getSelection(lBox)
    if (length(l) == 0){
      errorCondition(recall=urcaablrtest, message="You must select a loading restriction matrix.")
      return()
    }
    rint <- tclvalue(rankVariable)
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    doItAndPrint(paste("summary(ablrtest(z = ", x, ", H = ", y, ", A = ", l,  ", r = ", rint, "))", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="ablrtest")
  rankFrame <- tkframe(top)
  rankVariable <- tclVar("1")
  rankField <- tkentry(rankFrame, width="2", textvariable=rankVariable)
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(getFrame(yBox), getFrame(lBox), sticky="nw")
  tkgrid(rankFrame, sticky="w")
  tkgrid(tklabel(rankFrame, text="Number of cointegrating relationships =  ", fg="blue"), rankField, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  tkgrid.configure(rankField, sticky="e")
  dialogSuffix(rows=3, columns=2)
}
#
# VECM with a structural shift
#
urcacajolst <- function(){
  dataSets <- listDataSets()
  if (length(dataSets) == 0){
    tkmessageBox(message="There are no data sets from which to choose.", icon="error", type="ok")
    tkfocus(CommanderWindow())
    return()
  }
  initializeDialog(title="VECM with level shift")
  xBox <- variableListBox(top, dataSets, title="Data Sets (pick one)")
  assign("UpdateModelNumber", UpdateModelNumber + 1, envir=.GlobalEnv)
  modelName <- tclVar(paste("VECMmodel.", UpdateModelNumber, sep=""))
  modelFrame <- tkframe(top)
  model <- tkentry(modelFrame, width="20", textvariable=modelName)
  onOK <- function(){
    x <- getSelection(xBox)
    if (length(x) == 0){
      errorCondition(recall=urcacajo, message="You must select a data set.")
      return()
    }
    ttype <- if(tclvalue(trendVariable) == 0) "FALSE" else "TRUE"
    lags <- tclvalue(lagVariable)
    seas <- tclvalue(seasonVariable)
    if (GrabFocus()) tkgrab.release(top)
    tkdestroy(top)
    modelValue <- tclvalue(modelName)
    if (!is.valid.name(modelValue)){
      assign("UpdateModelNumber", UpdateModelNumber - 1, envir=.GlobalEnv)
      errorCondition(recall=urcacajolst, message=paste('"', modelValue, '" is not a valid name.', sep=""))
      return()
    }
    if (is.element(modelValue, listVECMmodels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type="Model"))){
        assign("UpdateModelNumber", UpdateModelNumber - 1, envir=.GlobalEnv)
        if (GrabFocus()) tkgrab.release(top)
        tkdestroy(top)
        urcacajolst()
        return()
      }
    }
    command <- paste("cajolst(x = ", x, ",  trend = ", ttype, ", K = ", lags, ", season = ", seas, ")", sep="")
    logger(paste(modelValue, " <- ", command, sep=""))
    assign(modelValue, justDoIt(command), envir=.GlobalEnv)
    doItAndPrint(paste("summary(", modelValue, ")", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="cajolst")
  radioButtons(name="season", buttons=c("none", "monthly", "quarterly"), values=c("NULL", "12", "4"), 
               labels=c("None", "Monthly seasonality", "Quarterly seasonality"), title="Seasonality")
  checkBoxes(frame="trendFrame", boxes="trend", initialValues="1", labels="Include linear trend in the auxiliary regressions?")
  lagFrame <- tkframe(top)
  lagVariable <- tclVar("2")
  lagField <- tkentry(lagFrame, width="2", textvariable=lagVariable)
  tkgrid(tklabel(modelFrame, text="Enter name for model:"), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(tklabel(lagFrame, text="Lag order =  ", fg="blue"), lagField, sticky="w")
  tkgrid(lagFrame, sticky="w")
  tkgrid(trendFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  tkgrid.configure(lagField, sticky="e")
  dialogSuffix(rows=5, columns=1)
}
#
# Utility Functions
#
listVECMmodels <- function(envir=.GlobalEnv, ...) {
  objects <- ls(envir=envir, ...)
  if (length(objects) == 0) NULL
  else objects[sapply(objects, function(.x) "ca.jo" == (class(eval(parse(text=.x), envir=envir))[1]))]
}

listMatrix <- function(envir=.GlobalEnv, ...) {
  objects <- ls(envir=envir, ...)
  if (length(objects) == 0) NULL
  else objects[sapply(objects, function(.x) "matrix" == (class(eval(parse(text=.x), envir=envir))[1]))]
}

listNumeric <- function(envir=.GlobalEnv, ...) {
  objects <- ls(envir=envir, ...)
  if (length(objects) == 0) NULL
  else objects[sapply(objects, function(.x) "numeric" == (class(eval(parse(text=.x), envir=envir))[1]))]
}
