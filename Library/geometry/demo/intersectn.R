library(geometry)

## 2D example
ps1 <- rbind(c(0,   sqrt(3)),
             c(3/2, -sqrt(3)/2),
             c(-3/2, -sqrt(3)/2))
ps2 <- ps1
ps2[,2] <- -ps2[,2]

is <-  intersectn(ps1, ps2)
plot(is, asp=1)

## 3D example
ps1a <- rbox(2, C=0.5)
dt1a <- delaunayn(ps1a)
ps1b <- rbox(2, C=0.5) + 2
dt1b <- delaunayn(ps1b)
ps1 <- rbind(ps1a, ps1b)
dt1 <- rbind(dt1a, dt1b + nrow(ps1a))
tetramesh(dt1, ps1, alpha=0.5, col="yellow")

ps2 <-  rbox(2, C=0.5) + 0.5
dt2 <- delaunayn(ps2)

tetramesh(dt2, ps2, alpha=0.5, col="red", clear=FALSE)

vol <- 0
for (i in 1:nrow(dt1)) {
  for (j in 1:nrow(dt2)) {
    is <- intersectn(ps1[dt1[i,],], ps2[dt2[j,],])
    vol <- vol + is$ch$vol
  }
}
message(paste("Volume of overlap should be 0.125. It is:", vol))

