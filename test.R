set_A <- subset(set_A, orig.ident == c("FVB11", "FVB12"))
levels(set_A$orig.ident)
set_A$orig.ident <- droplevels(set_A$orig.ident)
