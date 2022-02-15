library(devtools)

document()

load_all()

check()

q0 = new("Qbone")
q0 = new("Qbone", version = packageVersion(pkg = "Qbone"))
?dorem

list1 = list(c("123", "345"), "234")
meta.data = data.frame(name = c("a_1", "b_1"))

qbonedata = CreateQboneData(list1, meta.data)

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
