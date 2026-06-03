## ----text-entry, fig.align='center', fig.width=6, out.width='75%'-------------
myTree <- ape::read.tree(text = "((A, B), ((C, D), (E, F)));")
plot(myTree)

## ----double-bracket, fig.align='center', fig.width=6, out.width='75%'---------
badTree <- ape::read.tree(text = "((A, B), (((C, D), ((E), F))));")
plot(badTree)
ape::nodelabels(bg = c(3, 3, 2, 2, 3, 3, 2))

## ----branch-lengths, fig.align='center', fig.width=6, out.width='75%'---------
myTree <- ape::read.tree(
  text = "((A:1, B:1):2, ((C:1, D:1):2, (E:1, F:1):2):4);"
)
plot(myTree)

