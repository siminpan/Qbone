library(devtools)

document()

load_all()

load_all(recompile = T)

check()

tinytex::tlmgr_install("makeindex")
devtools::build_manual(pkg = "/home/span/Documents/Qbone")

q0 = new("Qbone")
q0 = new("Qbone", version = packageVersion(pkg = "Qbone"))
?dorem

list1 = list(c("123", "345"), "234")
meta.data = data.frame(name = c("a_1", "b_1"))
rownames(meta.data) = meta.data[,1]

qbonedata = createQboneData(list1, meta.data)

qbone = createQboneObject(qbonedata, meta.data = meta.data)

qbone <- addMetaData(object = qbone, metadata = meta.data)
qbone[["name"]] <- meta.data

cqo2 = addMetaData(qbone, c("c", "d", "e"), col.name = "names")

data.frame(row.names = qbonedata@data)
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
