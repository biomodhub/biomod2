##
##  p r a c m a . R  'pracma' Demos
##

require("pracma")

# Remember graphics settings
opar <- par(ask = dev.interactive(orNone = TRUE))

if (interactive())
    readline("Press ENTER to continue: ")
cat("Some demos will follow shortly. BYE.\n")

# Reset graphics parameters
par(opar)
