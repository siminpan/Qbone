library(devtools)

document()

load_all()

load_all(recompile = T)

check()

devtools::check(manual=TRUE)

tinytex::install_tinytex()
tinytex::tlmgr_install("makeindex")
tinytex::tlmgr_update()
tinytex::reinstall_tinytex()
library("tinytex")
devtools::build_manual(pkg = ".")
# readQbone ----
q0 = new("Qbone")
q0 = new("Qbone", version = packageVersion(pkg = "Qbone"))
?dorem

list1 = list(c("123", "345"), "234")
meta.data = data.frame(name = c("a_1", "b_1"), names = c("c", "b"))
rownames(meta.data) = meta.data[,1]

qbonedata = createQboneData(list1, meta.data, sampleid = 2)

qbone1 = createQboneObject(qbonedata, meta.data = meta.data)

qbone2 = createQboneObject(list1, meta.data = meta.data, sampleid = 2)

qbone <- addMetaData(object = qbone, metadata = meta.data)
qbone[["name"]] <- meta.data

cqo2 = addMetaData(qbone, c("c", "d", "e"), col.name = "names")

data.frame(row.names = qbonedata@data)

idents(q2) <- c("c", "d")
idents(q2)

# lassolist ----
data.dir = "/home/span/Documents/MOSJ-3DCT/data/csv.test"
q1 = readQbone(data.dir, groupbyfolder = T)

qa1 = getQboneData(q1, slot = 'data', assay = defaultAssay(q1))


q2 = thinData(q1,prop=0.0001)
q3 = lassolist2(q2)
q4 = lassolist2(q2)
q5 = lassolist2(q2)

document()
q6 = lassolist2(q2)
all.equal(q6,q4)

## old ----
raw.dataset <- getQboneData(q2, slot = 'data', assay = defaultAssay(q2))
a1 = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1))
a2 = c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1))
list1 = runlassolist(raw.dataset[[1]])

## assay.seed ----
q7 = lassolist(q2)
q8 = lassolist(q2, assay.seed = q7@assays[["Lasso"]]@scale.data[["lassolist"]])
all.equal(q8@assays[["Lasso"]]@data,
          q7@assays[["Lasso"]]@data)

set.seed(12345)
q9 = lassolist(q2)
set.seed(12345)
q10 = lassolist(q2)
all.equal(q9@assays[["Lasso"]]@data,
          q10@assays[["Lasso"]]@data)

# For this example, set the random seed
set.seed(423)
runif(3)
#> [1] 0.1089715 0.5973455 0.9726307

# Save the seed
set.seed(423)
oldseed <- .Random.seed
oldseed1 <- .Random.seed
.Random.seed <- oldseed1
oldseed2 <- .Random.seed

all.equal(oldseed,oldseed1)
all.equal(oldseed2,oldseed1)

runif(3)
.Random.seed <- q7@assays[["Lasso"]]@scale.data[["lassolist"]]

# Get the same random numbers as before, after saving the seed
runif(3)



## performance ----
library("microbenchmark")
document()
q7 = lassolist(q2)
q7.2 = lassolist2(q2, verbose=F)
microbenchmark(lassolist(q2), lassolist2(q2), times = 10L)

# eigen mm ----
one <- rep(1, x.n)
one1 = matrix(rep(1, x.n), ncol = x.n)

x3 = one %*% object.x
x4 = eigenmapmm(one1, object.x)

# idents ----
# sloop::ftype(packageVersion(x = "0.0.0.9000"))
class(package_version(x = '2.99.0'))
idents(cqo2) <- c("a1", "b1")
samples(cqo2)

idents <- factor(x = unlist(x = lapply(
  # X = colnames(x = meta.data),
  X = names(data),
  FUN = ExtractField,
  field = field,
  delim = delim
)))


# other ----
object <- new(
  Class = 'Qbone',
  assays = assay.list,
  meta.data =  data.frame(row.names = names(data@data)),
  active.assay = assay,
  active.ident = idents,
  project.name = project,
  version = packageVersion(pkg = 'Qbone'))


mate = c("c", "d")
names(mate) = c("a_1", "b_1")
cqo2 = addMetaData(qbone, mate, col.name = "names")
assays(qbonedata)
assays <- FilterObjects(object = qbone, classes.keep = 'QboneData')
names(getQboneData(qbonedata))

foo1 <- function(verbose =T ){
  for (i in 1:100){
    if (verbose) {
      message("Calculating cluster ", i)
    }
    Sys.sleep(1)
  }

}
for (i in 1:100){
    message("Calculating cluster ", i)
  }

call1 = function(x,y){print(x)}
