library(devtools)

document()

load_all()

check()

q0 = new("Qbone")
q0 = new("Qbone", version = packageVersion(pkg = "Qbone"))


# sloop::ftype(packageVersion(x = "0.0.0.9000"))
class(package_version(x = '2.99.0'))
