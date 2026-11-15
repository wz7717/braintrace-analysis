# Native SingleR 2.14.2 implementation; explicit real/synthetic invocation only.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) %% 2L || !length(args)) stop("Expected option/value pairs")
opt <- as.list(args[seq.int(2L,length(args),2L)])
names(opt) <- args[seq.int(1L,length(args),2L)]
mode <- opt[["--mode"]]
if (!identical(mode,"synthetic") && !identical(mode,"real")) {
    cat("Explicit real or synthetic mode required\n",file=stderr())
    quit(save="no",status=77L)
}
is_real <- identical(mode,"real")
execution_scope <- if (is_real) "REAL_FINAL_EXTERNAL" else "SYNTHETIC_ONLY"
action <- opt[["--action"]]
output <- opt[["--out"]]
dir.create(output,recursive=TRUE,showWarnings=FALSE)
abort_contract <- function(code,detail="") {
    utils::write.table(data.frame(status="METHOD_CONTRACT_FAILURE",code=code,detail=detail),
        file.path(output,"contract_status.tsv"),sep="\t",quote=TRUE,row.names=FALSE,fileEncoding="UTF-8")
    cat(paste0(code,": ",detail,"\n"),file=stderr())
    quit(save="no",status=3L)
}
write_contract_ok <- function() utils::write.table(data.frame(status="COMPLETED",code="OK",detail=""),
    file.path(output,"contract_status.tsv"),sep="\t",quote=TRUE,row.names=FALSE,fileEncoding="UTF-8")
writeLines(paste0(execution_scope,"; real_data_execution_enabled=",if (is_real) "true" else "false"),file.path(output,"scope.txt"))
tryCatch({
    library_dir <- normalizePath(Sys.getenv("BRAINTRACE_B17_R_LIBRARY"),winslash="/",mustWork=TRUE)
    .libPaths(c(library_dir,.Library))
    if (getRversion() != "4.6.1") abort_contract("R_VERSION_MISMATCH")
    suppressPackageStartupMessages(library(SingleR))
    if (as.character(utils::packageVersion("SingleR")) != "2.14.2") abort_contract("SINGLER_VERSION_MISMATCH")
    package_path <- normalizePath(find.package("SingleR"),winslash="/",mustWork=TRUE)
    if (dirname(package_path) != library_dir) abort_contract("SINGLER_LIBRARY_OUTSIDE_EXISTING_ISOLATION")
    thread_env <- c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS")
    if (any(Sys.getenv(thread_env,unset="") != "1")) abort_contract("CHILD_THREAD_ENV_NOT_LOCKED")
    options(digits=17)
    utils::capture.output(list(scope=execution_scope,mode=mode,synthetic=!is_real,action=action,R_home=R.home(),libPaths=.libPaths(),
        SingleR_path=package_path,SingleR_version=utils::packageVersion("SingleR"),sessionInfo=utils::sessionInfo(),
        loaded_DLLs=lapply(getLoadedDLLs(),function(x)x[["path"]]),thread_environment=Sys.getenv(thread_env)),
        file=file.path(output,"runtime.dput"))
    if (identical(action,"build")) {
        input <- opt[["--reference-dir"]]
        genes <- readLines(file.path(input,"genes.txt"),warn=FALSE,encoding="UTF-8")
        supported <- c("4009","4132","4180","4219","4275","4391","4696","9001","9512")
        # Shared reader emits unquoted TAB-separated canonical JSON ref IDs.
        reference <- as.matrix(utils::read.delim(file.path(input,"reference.tsv"),row.names=1L,
            check.names=FALSE,quote="",comment.char="",stringsAsFactors=FALSE))
        storage.mode(reference) <- "double"
        meta <- utils::read.delim(file.path(input,"reference_meta.tsv"),check.names=FALSE,
            quote="",comment.char="",stringsAsFactors=FALSE)
        labels <- as.character(meta$label)
        if (length(genes)<2L || anyDuplicated(genes) || !identical(rownames(reference),genes) ||
            anyDuplicated(colnames(reference)) || !identical(colnames(reference),as.character(meta$ref_id)) ||
            nrow(meta)!=ncol(reference) || anyNA(labels) || any(!is.finite(reference)) ||
            !setequal(unique(labels),supported)) abort_contract("INVALID_FIXED_INPUT_AXES_OR_REFERENCE")
        utils::capture.output(list(train=formals(SingleR::trainSingleR),classify=formals(SingleR::classifySingleR)),
            file=file.path(output,"native_formals.dput"))
        markers_pre <- SingleR::getClassicMarkers(reference,labels=labels,
            assay.type="logcounts",check.missing=TRUE,de.n=NULL,num.threads=1L)
        marker_union_pre <- unique(unlist(markers_pre,use.names=FALSE))
        if (length(marker_union_pre)<2L || anyDuplicated(marker_union_pre) || !all(marker_union_pre %in% genes))
            abort_contract("MARKER_UNION_LT2_OR_INVALID")
        marker_indices <- match(marker_union_pre,genes)
        for (i in seq_len(ncol(reference))) {
            values <- reference[marker_indices,i]
            if (any(!is.finite(values)) || max(values)<=min(values))
                abort_contract("INVALID_TRAINING_REFERENCE_MARKER_NORM",colnames(reference)[i])
        }
        trained <- SingleR::trainSingleR(ref=reference,labels=labels,
            test.genes=genes,genes="de",de.method="classic",de.n=NULL,de.args=list(),
            aggr.ref=FALSE,restrict=genes,assay.type="logcounts",check.missing=TRUE,
            approximate=FALSE,num.threads=1L)
        actual_markers <- as.character(trained$markers$unique)
        if (length(actual_markers)<2L || anyDuplicated(actual_markers) || !setequal(actual_markers,marker_union_pre) ||
            !all(actual_markers %in% genes) || !setequal(trained$labels$unique,supported) ||
            !identical(trained$options$test.genes,genes)) abort_contract("TRAINED_MARKER_OR_AXIS_STRUCTURE_MISMATCH")
        utils::write.table(data.frame(gene_id=actual_markers),file.path(output,"marker_genes.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE,fileEncoding="UTF-8")
        utils::capture.output(list(scope=execution_scope,mode=mode,synthetic=!is_real,training_reference_columns=colnames(reference),
            reference_labels=labels,unique_labels=trained$labels$unique,genes=genes,
            markers_unique=actual_markers,markers_full=trained$markers$full,preflight_markers_full=markers_pre,
            options=trained$options,trained_top_level_types=lapply(trained,typeof)),file=file.path(output,"trained_structure.dput"))
        utils::capture.output(list(train=list(genes="de",de.method="classic",de.n=NULL,de.args=list(),
            aggr.ref=FALSE,restrict=genes,test.genes=genes,check.missing=TRUE,approximate=FALSE,num.threads=1L),
            classify=list(quantile=0.8,fine.tune=FALSE,prune=FALSE,num.threads=1L),
            omitted_BNPARAM_BPPARAM="Omitted as in B8; pinned release formals NULL; no backend replacement",
            reference_scale="Admitted final PLM values unchanged; ordinary-matrix logcounts declaration does not transform"),
            file=file.path(output,"actual_call_parameters.dput"))
        saveRDS(list(trained=trained,genes=genes,supported=supported,actual_markers=actual_markers,
                    reference_columns=colnames(reference),scope=execution_scope,mode=mode,synthetic=!is_real),file.path(output,"trained.rds"))
        write_contract_ok()
        cat("COMPLETED\tbuild\n")
    } else if (identical(action,"score")) {
        cached <- readRDS(file.path(opt[["--cache-dir"]],"trained.rds"))
        if (!identical(cached$scope,execution_scope) || (is_real &&
            (!identical(cached$mode,"real") || !identical(cached$synthetic,FALSE))))
            abort_contract("CACHE_SCOPE_MODE_MISMATCH")
        genes <- cached$genes; supported <- cached$supported; trained <- cached$trained
        query <- as.matrix(utils::read.delim(opt[["--query"]],row.names=1L,check.names=FALSE,
            quote="",comment.char="",stringsAsFactors=FALSE,na.strings=c("NaN","NA")))
        storage.mode(query) <- "double"
        if (ncol(query)!=1L || !identical(rownames(query),genes)) abort_contract("QUERY_FIXED_G_AXIS_MISMATCH")
        marker_indices <- match(cached$actual_markers,genes)
        status <- data.frame(record_id=colnames(query),status="NOT_RUN",reason="",native_label="",stringsAsFactors=FALSE)
        native_scores <- NULL
        if (any(!is.finite(query))) {
            status$status <- "QUERY_FAILURE"; status$reason <- "QUERY_NONFINITE_FIXED_G"
        } else if (max(query[marker_indices,1L])<=min(query[marker_indices,1L])) {
            status$status <- "QUERY_FAILURE"; status$reason <- "QUERY_CONSTANT_ON_NATIVE_MARKERS"
        } else {
            result <- tryCatch(SingleR::classifySingleR(test=query,trained=trained,quantile=0.8,
                fine.tune=FALSE,prune=FALSE,num.threads=1L),error=function(e)e)
            if (inherits(result,"error")) {
                status$status <- "QUERY_FAILURE"; status$reason <- paste0("QUERY_NATIVE_ERROR:",conditionMessage(result))
            } else {
                native_scores <- as.matrix(result$scores)
                if (!identical(rownames(result),colnames(query)) || nrow(native_scores)!=1L ||
                    !identical(colnames(native_scores),as.character(trained$labels$unique)) ||
                    !setequal(colnames(native_scores),supported)) abort_contract("NATIVE_SCORE_AXIS_MISMATCH")
                if ("pruned.labels" %in% colnames(result)) abort_contract("UNEXPECTED_PRUNED_LABELS")
                if (!setequal(S4Vectors::metadata(result)$common.genes,cached$actual_markers)) abort_contract("NATIVE_SCORE_MARKER_AXIS_MISMATCH")
                if (any(!is.finite(native_scores))) {
                    status$status <- "QUERY_FAILURE"; status$reason <- "QUERY_UNDEFINED_NATIVE_SCORE"; native_scores <- NULL
                } else {
                    status$status <- "OK"; status$native_label <- as.character(result$labels[1L])
                }
            }
        }
        if (!is.null(native_scores)) {
            native_table <- data.frame(record_id=colnames(query),stringsAsFactors=FALSE)
            for (label in colnames(native_scores)) native_table[[label]] <- sprintf("%.17g",native_scores[,label])
            utils::write.table(native_table,file.path(output,"native_scores.tsv"),sep="\t",quote=FALSE,row.names=FALSE,fileEncoding="UTF-8")
        }
        utils::write.table(status,file.path(output,"native_query_status.tsv"),sep="\t",quote=TRUE,row.names=FALSE,fileEncoding="UTF-8")
        write_contract_ok()
        cat("COMPLETED\tscore\n")
    } else abort_contract("UNKNOWN_ACTION")
},error=function(e) abort_contract("NATIVE_TRAINING_OR_STRUCTURE_ERROR",conditionMessage(e)))
