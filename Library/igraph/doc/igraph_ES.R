## ----echo = TRUE, eval = FALSE------------------------------------------------
# install.packages("igraph")

## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(fig.width=6, fig.height=6)

## ----setup--------------------------------------------------------------------
library("igraph")

## -----------------------------------------------------------------------------
g <- make_empty_graph()

## -----------------------------------------------------------------------------
g <- make_graph(edges = c(1,2, 1,5), n=10, directed = FALSE)

## ----echo = TRUE--------------------------------------------------------------
g <- make_graph(~ 1--2, 1--5, 3, 4, 5, 6, 7, 8, 9, 10)

## ----echo = TRUE--------------------------------------------------------------
g

## ----echo = TRUE--------------------------------------------------------------
summary(g)

## ----echo = TRUE--------------------------------------------------------------
g <- make_graph("Zachary")

## -----------------------------------------------------------------------------
plot(g)

## -----------------------------------------------------------------------------
g <- add_vertices(g, 3)

## -----------------------------------------------------------------------------
g <- add_edges(g, edges = c(1,35, 1,36, 34,37))

## ----echo = TRUE, eval=FALSE--------------------------------------------------
# g <- g + edges(c(1,35, 1,36, 34,37))

## ----echo = TRUE, error = TRUE------------------------------------------------
try({
g <- add_edges(g, edges = c(38, 37))
})

## ----echo = TRUE--------------------------------------------------------------
g <- g %>%
  add_edges(edges = c(1, 34)) %>%
  add_vertices(3) %>%
  add_edges(edges = c(38, 39, 39, 40, 40, 38, 40, 37))
g

## ----echo = TRUE--------------------------------------------------------------
edge_id_para_borrar <- get_edge_ids(g, c(1,34))
edge_id_para_borrar

## -----------------------------------------------------------------------------
g <- delete_edges(g, edge_id_para_borrar)

## ----echo = TRUE--------------------------------------------------------------
g <- make_ring(10) %>% delete_edges("10|1")
plot(g)

## ----echo = TRUE--------------------------------------------------------------
g <- make_ring(5)
g <- delete_edges(g, get_edge_ids(g, c(1,5, 4,5)))
plot(g)

## -----------------------------------------------------------------------------
g1 <- graph_from_literal(
  A-B:C:I,
  B-A:C:D,
  C-A:B:E:H,
  D-B:E:F,
  E-C:D:F:H,
  F-D:E:G,
  G-F:H,
  H-C:E:G:I,
  I-A:H
)
plot(g1)

## ----echo = TRUE--------------------------------------------------------------
is_chordal(g1, fillin=TRUE)

## ----echo = TRUE--------------------------------------------------------------
chordal_graph <- add_edges(g1, is_chordal(g1, fillin=TRUE)$fillin)
plot(chordal_graph)

## ----echo = TRUE--------------------------------------------------------------
graph1 <- make_tree(127, 2, mode = "undirected")
summary(g)

## -----------------------------------------------------------------------------
graph2 <- make_tree(127, 2, mode = "undirected")

## ----echo = TRUE--------------------------------------------------------------
identical_graphs(graph1, graph2)

## ----echo = TRUE--------------------------------------------------------------
graph1 <- sample_grg(100, 0.2)
summary(graph1)

## ----echo = TRUE--------------------------------------------------------------
graph2 <- sample_grg(100, 0.2)
identical_graphs(graph1, graph2)

## ----echo = TRUE--------------------------------------------------------------
isomorphic(graph1, graph2)

## -----------------------------------------------------------------------------
g <- make_graph(
  ~ Alice-Boris:Himari:Moshe,
  Himari-Alice:Nang:Moshe:Samira,
  Ibrahim-Nang:Moshe,
  Nang-Samira
)

## ----echo = TRUE--------------------------------------------------------------
V(g)$age <- c(25, 31, 18, 23, 47, 22, 50)
V(g)$gender <- c("f", "m", "f", "m", "m", "f", "m")
E(g)$is_formal <- c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE)
summary(g)

## ----echo = TRUE, eval=FALSE--------------------------------------------------
# g <- make_graph(
#   ~ Alice-Boris:Himari:Moshe,
#   Himari-Alice:Nang:Moshe:Samira,
#   Ibrahim-Nang:Moshe,
#   Nang-Samira
# ) %>%
#   set_vertex_attr("age", value = c(25, 31, 18, 23, 47, 22, 50)) %>%
#   set_vertex_attr("gender", value = c("f", "m", "f", "m", "m", "f", "m")) %>%
#   set_edge_attr("is_formal", value = c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE))
# summary(g)

## ----echo = TRUE--------------------------------------------------------------
E(g)$is_formal
E(g)$is_formal[1] <- TRUE
E(g)$is_formal

## ----echo = TRUE--------------------------------------------------------------
g$date <- c("2022-02-11")
graph_attr(g, "date")

## ----echo = TRUE--------------------------------------------------------------
match(c("Ibrahim"), V(g)$name)

## ----echo = TRUE--------------------------------------------------------------
V(g)$name[1:3] <- c("Alejandra", "Bruno", "Carmina")
V(g)

## ----echo = TRUE--------------------------------------------------------------
g <- delete_vertex_attr(g, "gender")
V(g)$gender

## ----echo = TRUE--------------------------------------------------------------
degree(g)

## ----echo = TRUE--------------------------------------------------------------
degree(g, 7)

## ----echo = TRUE--------------------------------------------------------------
degree(g, v = c(3,4,5))

## ----echo = TRUE--------------------------------------------------------------
degree(g, v = c("Carmina", "Moshe", "Nang"))

## ----echo = TRUE--------------------------------------------------------------
degree(g, "Bruno")

## ----echo = TRUE--------------------------------------------------------------
edge_betweenness(g)

## ----echo = TRUE--------------------------------------------------------------
ebs <- edge_betweenness(g)
as_edgelist(g)[ebs == max(ebs), ]

## ----echo = TRUE--------------------------------------------------------------
which.max(degree(g))

## ----echo = TRUE--------------------------------------------------------------
graph <- graph.full(n=10)
only_odd_vertices <- which(V(graph)%%2==1)
length(only_odd_vertices)

## ----echo = TRUE--------------------------------------------------------------
seq <- V(graph)[2, 3, 7]
seq

## ----echo = TRUE--------------------------------------------------------------
seq <- seq[1, 3]    # filtrar un conjunto de vértices existente
seq

## ----echo = TRUE, error = TRUE------------------------------------------------
try({
seq <- V(graph)[2, 3, 7, "foo", 3.5]
})

## ----echo = TRUE--------------------------------------------------------------
V(g)[age < 30]$name

## ----echo = TRUE--------------------------------------------------------------
`%notin%` <- Negate(`%in%`)

## ----echo = TRUE--------------------------------------------------------------
V(g)$degree <- c("A", "B", "B+", "A+", "C", "A", "B")
V(g)$degree[degree(g) == 3]

## ----echo = TRUE--------------------------------------------------------------
V(g)$name[degree(g) == 3]

## ----echo = TRUE, warning = FALSE---------------------------------------------
E(g)[.from(3)]

## ----echo = TRUE, warning = FALSE---------------------------------------------
E(g)[.from("Carmina")]

## ----echo = TRUE--------------------------------------------------------------
E(g) [ 3:5 %--% 5:6 ]

## -----------------------------------------------------------------------------
V(g)$gender <- c("f", "m", "f", "m", "m", "f", "m")

## ----echo = TRUE--------------------------------------------------------------
men <- V(g)[gender == "m"]$name
men

## ----echo = TRUE--------------------------------------------------------------
women <- V(g)[gender == "f"]$name
women

## ----echo = TRUE--------------------------------------------------------------
E(g)[men %--% women]

## ----echo = TRUE--------------------------------------------------------------
as_adjacency_matrix(g)

## -----------------------------------------------------------------------------
layout <- layout_with_kk(g)

## -----------------------------------------------------------------------------
layout <- layout_as_tree(g, root = 2)

## -----------------------------------------------------------------------------
layout <- layout_with_kk(g)

## -----------------------------------------------------------------------------
plot(g, layout = layout, main = "Red social con el algoritmo de diseño Kamada-Kawai")

## -----------------------------------------------------------------------------
plot(
  g,
  layout = layout_with_fr,
  main = "Red social con el algoritmo de disposición Fruchterman-Reingold"
)

## -----------------------------------------------------------------------------
V(g)$color <- ifelse(V(g)$gender == "m", "yellow", "red")
plot(
  g,
  layout = layout,
  vertex.label.dist = 3.5,
  main = "Red social - con los géneros como colores"
)

## -----------------------------------------------------------------------------
plot(
  g,
  layout = layout,
  vertex.label.dist = 3.5,
  vertex.color = as.factor(V(g)$gender))

## -----------------------------------------------------------------------------
plot(
  g,
  layout = layout,
  vertex.label.dist = 3.5,
  vertex.size = 20,
  vertex.color = ifelse(V(g)$gender == "m", "yellow", "red"),
  edge.width = ifelse(E(g)$is_formal, 5, 1)
)

## ----session-info-------------------------------------------------------------
sessionInfo()

