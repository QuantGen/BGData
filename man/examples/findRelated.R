# Load example data
bg <- BGData:::loadExample()

G <- getG(bg@geno)
findRelated(G)
