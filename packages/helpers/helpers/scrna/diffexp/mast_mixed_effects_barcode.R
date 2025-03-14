function(adata_, group, reference){
    # create a MAST object
    # use a quarter of all cores
    options(mc.cores=as.integer(parallel::detectCores() * .25))

    print("generating SingleCellAssay")
    sca <- SceToSingleCellAssay(adata_, class = "SingleCellAssay", check_sanity = FALSE)

    # add a column to the data which contains scaled number of genes that are expressed in each cell
    colData(sca)$ngeneson <- scale(colSums(assay(sca)>0))

    # store the columns that we are interested in as factors
    label <- factor(colData(sca)$label)
    barcode <- factor(colData(sca)$barcode)
    # set the reference level
    colData(sca)$label <- relevel(label, reference)
    colData(sca)$barcode <- barcode

    # define and fit the model

    print("generating model")
    zlmCond <- zlm(
       formula = ~ngeneson + label + (1 | barcode),
       sca = sca,
       method = 'glmer',
       ebayes = FALSE,
       strictConvergence=FALSE,
       parallel=TRUE,
       fitArgsD=list(nAGQ = 0)
    )
    print("performing likelihood test")
    # perform likelihood-ratio test for the condition that we are interested in
    summaryCond <- summary(zlmCond, doLRT=group)

    # get the table with log-fold changes and p-values
    summaryDt <- summaryCond$datatable
    result <- merge(
        summaryDt[contrast==group & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
        summaryDt[contrast==group & component=='logFC', .(primerid, coef)], # logFC coefficients
        by='primerid'
    )

    # MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
    result[,coef:=result[,coef]/log(2)]

    # do multiple testing correction
    result[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]

    result <- stats::na.omit(as.data.frame(result))
    return(result)
}