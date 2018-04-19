library(arules)

############################### FUNCTIONS ###############################

allE <- list(0.1, 0.05, 0.01)
allD <- list(0.01, 0.001)
allDat <- list(
  "chess.dat",
  "mushroom.dat",
  "T10I4D100K.txt",
  "retail.dat",
  "connect.dat",
  "pumsb_star.dat",
  "pumsb.dat",
  "T40I10D100K.dat",
  "kosarak.dat",
  "accidents.dat"
)

# Load
loadExperiment <- function(db) {
  for (e in allE) { for (d in allD) { load(sprintf("output/%s/(ε=%g)(δ=%g)", db, e, d), envir = globalenv()) }}
}

# Print C-like formatted strings
printf <- function(...) cat(sprintf(...))

# Convert a configuration to a filename
configToPath <- function(c) {
  sprintf("output/%s/(ε=%g)(δ=%g)", c$dataset, c$e, c$delta)
}
configToStr <- function(c) {
  sprintf("(ε=%g)(δ=%g)", c$e, c$delta)
}

# Compute VC dimension (i.e. d-bound)
d.bound <- function(db, D) {
  printf("Calculation VC dimension... ")
  xs <- lapply(as(db, "list"), function (x) {length(x)})
  rest <- 0
  for (q in D:1) {
    occur <- sum(xs == q)
    if (occur + rest >= q) {
      printf("VC: %d\n", q)
      return(q)
    }
    else {rest <- rest + occur}
  }
  stop("No VC!")
}

# Mine FIs
mine <- function(filename, support) {
  printf("\nMining... ")
  apriori( filename
         , parameter = list
             ( target = "frequent itemsets"
             , support = support
             #, maxlen = 100
             )
         , control = list(verbose = FALSE)
         )
}

# Calculate absolute frequencies
freqs <- function(db, itemset) {
  mapply(function (x, y) abs(x - y),
         itemset@quality$support,
         support(itemset, db))
}

# Min size
minSize <- function(m1, m2) {
  if (m2 > m1) {
    warning("sample size > actual size")
    m1
  }
  else m2
}

run <- function(conf) {
############################### MAIN ###############################

printf("\n\n==================================================================================\n")
DATASET <- conf$dataset; c <- conf$c; delta <- conf$delta; e <- conf$e; mu <- conf$mu
printf("Dataset: %s\n", DATASET)
printf("c: %g, δ:%g, ε:%g, μ:%g\n", c, delta, e, mu)

# Read dataset
printf("\nReading dataset...\n")
dataset <- read.transactions(DATASET, format="basket")
D <- as(summary(dataset)@lengthSummary[6][[1]], "integer")
L <- length(dataset)
I <- dim(dataset)[2]
printf("|D|: %d\n", L)
printf("|I|: %d\n", I)
printf("Δ: %d\n", D)

# Configure (automatically pick θ)
is <- itemFrequency(dataset)
theta <- median(is[is > mean(is)])
printf("θ: %f {mean:%f ~ max:%f}\n", theta, mean(is), max(is))

# Mine FIs
fi <- mine(dataset, theta)
fi.length <- length(fi)
printf("|FI|: %d\n", fi.length)

## Toivonen ##
printf("\n-------- TOIVONEN --------\n")
# Config
T.sampleSize <- minSize(L, floor((1/(2*e^2))*log(2/delta)))
T.theta <- max(0.01, theta - sqrt((1/(2*T.sampleSize))*log(1/mu)))
printf("T.θ: %f\n", T.theta)
printf("T.sampleSize: %d\n", T.sampleSize)
# Sample
T.sample <- sample(dataset, T.sampleSize, replace = TRUE)
# Mine
T.fi <- mine(T.sample, T.theta)
T.fi.length <- length(T.fi)
printf("|T.FI|: %d\n", T.fi.length)
# Metrics
T.true.positives <- length(intersect(fi, T.fi))
T.false.negatives <- length(setdiff(fi, T.fi))
T.false.positives <- length(setdiff(T.fi, fi))
printf("\nT.+: %d, T.--: %d, T.-+: %d\n", T.true.positives, T.false.negatives, T.false.positives)
T.precision <- T.true.positives / (T.true.positives + T.false.positives)
T.recall <- T.true.positives / (T.true.positives + T.false.negatives)
printf("T.precision: %f\nT.recall: %f\n", T.precision, T.recall)
T.freq <- freqs(dataset, head(T.fi, n = fi.length))
T.avgFreq <- mean(T.freq)
T.maxFreq <- max(T.freq)
printf("T.Freq: {mean: %f} {max: %f}\n", T.avgFreq, T.maxFreq)

## Riondato ##
printf("\n-------- RIONDATO --------\n")
# Config
VC <- d.bound(dataset, D)
R.sampleSize <- minSize(L, floor(((4*c)/(e^2))*(VC + log(1/delta))))
printf("R.sampleSize: %d\n", R.sampleSize)
# Sample
R.sample <- sample(dataset, R.sampleSize, replace = TRUE)
# Mine
R.fi <- mine(R.sample, theta)
R.fi.length <- length(R.fi)
printf("|R.FI|: %d\n", R.fi.length)
# Metrics
R.true.positives <- length(intersect(fi, R.fi))
R.false.negatives <- length(setdiff(fi, R.fi))
R.false.positives <- length(setdiff(R.fi, fi))
printf("\nR.+: %d, R.--: %d, R.-+: %d\n", R.true.positives, R.false.negatives, R.false.positives)
R.precision <- R.true.positives / (R.true.positives + R.false.positives)
R.recall <- R.true.positives / (R.true.positives + R.false.negatives)
printf("R.precision: %f\nR.recall: %f\n", R.precision, R.recall)
R.freq <- freqs(dataset, head(R.fi, n = fi.length))
R.avgFreq <- mean(R.freq)
R.maxFreq <- max(R.freq)
printf("R.Freq: {mean: %f} {max: %f}\n", R.avgFreq, R.maxFreq)

############################### EXPORT ###############################

# Export variable
exportName = configToStr(config)
exportPath = configToPath(config)
assign(
  exportName,
  list(
    config=conf,
    "|Δ|"=D,
    "|D|"=L,
    "|I|"=I,
    VC=VC,
    "θ"=theta,
    "Φ"=fi.length,
    T=list("θ"=T.theta,
           "Σ"=T.sampleSize,
           "Φ"=T.fi.length,
           "+"=T.true.positives,
           "--"=T.false.negatives,
           "-+"=T.false.positives,
           precision=T.precision,
           recall=T.recall,
           avgFreq=T.avgFreq,
           maxFreq=T.maxFreq
           ),
    R=list("Σ"=R.sampleSize,
           "Φ"=R.fi.length,
           "+"=R.true.positives,
           "--"=R.false.negatives,
           "-+"=R.false.positives,
           precision=R.precision,
           recall=R.recall,
           avgFreq=R.avgFreq,
           maxFreq=R.maxFreq
    )
  )
)

# Save results to disk
save(list=c(exportName), file=exportPath)


printf("==================================================================================\n\n")
}

############################### CONFIG ###############################
for (dataset in allDat) {
  for (e in allE) {
    for (delta in allD) {
      mu <- delta
      config <- list(dataset = dataset,
                     e = e,
                     delta = delta,
                     mu = mu,
                     c = 0.5)
      run(config)
    }
  }
}
