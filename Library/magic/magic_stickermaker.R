## Notes:  does not seem to work on linux; knotR.png produced on windows
## Function magicplot2() is a bespoke version of magicplot()

library("magic")
library("hexSticker")

`magicplot2` <-
function (m) 
{
    par(pch = 16)
    n <- nrow(m)
    jj <- sort(t(m[n:1, ]), index.return = TRUE)$ix
    x <- process(jj, n)
    y <- (jj - 1)%/%n
    par(pty = "s", xaxt = "n", yaxt = "n")
    plot(x, y, type = "n", asp = 1, xlab = "", ylab = "", frame = FALSE)
    points(x, y, type="l",lwd=14)
}


png(file="magic_icon.png",width=1000,height=1000,bg="transparent")
magicplot2(magic(4))
dev.off()


sticker("magic_icon.png", package="magic", p_size=24, s_x=0.975, s_y=1.0,
s_width=0.83,asp=sqrt(3)/2, white_around_sticker=TRUE, h_fill="#7733FF",
h_color="#000000", filename="magic.png")

