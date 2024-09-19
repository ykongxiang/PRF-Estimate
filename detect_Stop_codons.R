detect_Stop_codons <- function(mrnaSeq,posit，new_qframe){
    x <- floor(length(mrnaSeq)-/3)
    stop_codons <- c("UAG", "UGA", "UAA")
    stop_codon_positions <- 0
    for (i in 1:x) {
        codon <- substr(mrna_seq[3*i-2,3*i])
        if (codon %in% stop_codons) {
            stop_codon_positions <- c(stop_codon_positions, 3i)
        }
        if (stop_codon_positions>posits){
            return(TRUE)
        }else{
            return(FALSE)
        }
    }
}