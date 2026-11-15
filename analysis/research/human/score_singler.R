args <- commandArgs(trailingOnly=TRUE)
out <- args[[1L]]
stopifnot(getRversion() == "4.6.2", packageVersion("SingleR") == "2.14.2")
options(digits=17)
read_matrix <- function(file) {
    x <- as.matrix(read.delim(file.path(out, file), row.names=1L, check.names=FALSE))
    storage.mode(x) <- "double"
    x
}
reference <- read_matrix("reference.tsv")
meta <- read.delim(file.path(out, "reference_labels.tsv"), colClasses="character")
labels <- meta$label
stopifnot(identical(meta$ref_id, colnames(reference)))
trained <- SingleR::trainSingleR(ref=reference, labels=labels,
    test.genes=rownames(reference), genes="de", de.method="classic", de.n=NULL,
    aggr.ref=FALSE, restrict=rownames(reference), check.missing=TRUE,
    approximate=FALSE, num.threads=1L)
markers <- trained$markers$unique
stopifnot(length(markers) >= 2L)
writeLines(markers, file.path(out, "marker_genes.txt"))
for (cohort in c("GTEx", "HPA")) {
    query <- read_matrix(paste0(cohort, "_query.tsv"))
    stopifnot(identical(rownames(query), rownames(reference)))
    finite <- colSums(!is.finite(query)) == 0L
    diverse <- apply(query[markers, , drop=FALSE], 2L,
        function(x) all(is.finite(x)) && max(x) > min(x))
    valid <- finite & diverse
    native <- data.frame(sample_id=colnames(query),
        reason=ifelse(finite, ifelse(diverse, "", "CONSTANT_ON_MARKERS"), "NONFINITE_QUERY"),
        native_label="")
    scores <- matrix(NA_real_, ncol(query), length(trained$labels$unique),
        dimnames=list(colnames(query), trained$labels$unique))
    if (any(valid)) {
        result <- SingleR::classifySingleR(test=query[, valid, drop=FALSE], trained=trained,
            quantile=0.8, fine.tune=FALSE, prune=FALSE, num.threads=1L)
        stopifnot(identical(rownames(result), colnames(query)[valid]))
        scores[valid, colnames(result$scores)] <- result$scores
        native$native_label[valid] <- as.character(result$labels)
        native$reason[valid & rowSums(!is.finite(scores)) != 0L] <- "NONFINITE_NATIVE_SCORE"
    }
    score_text <- as.data.frame(scores, check.names=FALSE)
    for (label in colnames(scores))
        score_text[[label]] <- ifelse(is.na(scores[, label]), NA_character_,
            sprintf("%.17g", scores[, label]))
    write.table(cbind(native, score_text), file.path(out, paste0(cohort, "_native.tsv")),
        sep="\t", quote=FALSE, row.names=FALSE)
}
capture.output(sessionInfo(), file=file.path(out, "R_SESSION.txt"))
