# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

sambada <- function(nomFichierParam, nomFichierEnv, nomFichierMarqueurs) {
    .Call(`_test2_sambada`, nomFichierParam, nomFichierEnv, nomFichierMarqueurs)
}

supervision <- function(nomFichier, numBlock, blockSize, maxDimension, selScore, scoreThreshold, sortScore, wordDelim) {
    .Call(`_test2_supervision`, nomFichier, numBlock, blockSize, maxDimension, selScore, scoreThreshold, sortScore, wordDelim)
}

recodePlink <- function(nbEch, nbSNP, nomFichierPlink, nomFichierSortie) {
    invisible(.Call(`_test2_recodePlink`, nbEch, nbSNP, nomFichierPlink, nomFichierSortie))
}

