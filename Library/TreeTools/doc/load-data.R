## ----Contrast matrix----------------------------------------------------------
contrast.matrix <- matrix(data = c(
# 0 1 -  # Each column corresponds to a character-state
  1, 0, 0, # Each row corresponds to a token, here 0, denoting the 
           # character-state set {0} 
  0, 1, 0, # 1 | {1}
  0, 0, 1, # - | {-}
  1, 1, 0, # A | {01}
  1, 1, 0, # + | {01}
  1, 1, 1  # ? | {01-}
), ncol = 3, # ncol corresponds to the number of columns in the matrix
byrow = TRUE)
dimnames(contrast.matrix) <- list(
  c(0, 1, "-", "A", "+", "?"), # A list of the tokens corresponding to each row
                               # in the contrast matrix
  c(0, 1, "-") # A list of the character-states corresponding to the columns 
               # in the contrast matrix
)

contrast.matrix

