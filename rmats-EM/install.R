## Installs rMATSISO package if it is not already installed

if (!require("rMATSISO")) {
    for (pkg in c("nloptr", "doParallel", "foreach")) {
        if (!require(pkg, character.only=T)) {
            install.packages(pkg)
        }
    }
    install.packages("./rMATSISO_1.0.0.tar.gz", type="source", repos=NULL)
}
