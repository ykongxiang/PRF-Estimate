detext_Slipseq <- function(frameStart,mrnaSeq){
    pattern <- "(A{3})(?=A{3})|(T{3})(?=T{3})|(C{3})(?=C{3})|(G{3})(?=G{3})"
    exists <- grepl(pattern, mrnaSeq[1,frameStart])
        return(exists)
}