sourceDir <- function(dir) {
  for (name in list.files(dir)) {
    source(file.path(dir, name))
  }
}