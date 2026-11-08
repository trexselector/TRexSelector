#!/usr/bin/env Rscript
# Builds the search index of the vignettes, the page titles and the headings with their anchors, and
# writes it into vignettes/vignette_nav.html between the two marker lines, where the search box of the
# shared header reads it. Run it after a heading or a page changes:
#
#   Rscript R/TRexSelector/tools/build_vignette_search_index.R
#
# The anchors follow pandoc's rule for automatic identifiers: lowercase, non-alphanumeric characters
# other than underscore, hyphen and period removed, runs of spaces to one hyphen, leading non-letters
# removed, "-1", "-2" appended to repeats within a page. Math in a heading counts with its TeX source,
# so "$p > n$" contributes "p-n".

vignettes <- file.path(dirname(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(),
                                                                                  value = TRUE))))),
                       "vignettes")
nav_file <- file.path(vignettes, "vignette_nav.html")

pandoc_id <- function(heading) {
  x <- tolower(heading)
  x <- gsub("[^[:alnum:]_. -]", "", x)
  x <- gsub("[[:space:]]+", "-", trimws(x))
  x <- sub("^[^[:alpha:]]+", "", x)
  if (x == "") "section" else x
}

# The visible text of a heading: inline code marks and trailing attributes removed.
heading_text <- function(line) {
  x <- sub("^#{1,3} +", "", line)
  x <- sub(" *\\{[^}]*\\} *$", "", x)
  x <- gsub("`", "", x)
  trimws(x)
}

page_entry <- function(file, number) {
  lines <- readLines(file, warn = FALSE)
  title <- sub("^title: *\"?(.*?)\"? *$", "\\1", grep("^title:", lines, value = TRUE)[1])
  in_code <- FALSE
  heads <- character(0)
  for (line in lines) {
    if (grepl("^```", line)) {
      in_code <- !in_code
      next
    }
    if (!in_code && grepl("^#{1,3} ", line)) heads <- c(heads, heading_text(line))
  }
  ids <- vapply(heads, pandoc_id, character(1), USE.NAMES = FALSE)
  seen <- character(0)
  for (i in seq_along(ids)) {
    base <- ids[i]
    k <- 0
    while (ids[i] %in% seen) {
      k <- k + 1
      ids[i] <- paste0(base, "-", k)
    }
    seen <- c(seen, ids[i])
  }
  list(page = sub("\\.Rmd$", ".html", basename(file)), number = number, title = title,
       heads = heads, ids = ids)
}

js_string <- function(x) paste0("\"", gsub("\"", "\\\\\"", gsub("\\\\", "\\\\\\\\", x)), "\"")

files <- c(file.path(vignettes, "TRexSelector.Rmd"),
           sort(list.files(vignettes, pattern = "^v[0-9]{2}-.*\\.Rmd$", full.names = TRUE)))
entries <- lapply(seq_along(files), function(i) page_entry(files[i], i - 1))
# one heading per line, so the lines of the header file stay within the ruler
rows <- vapply(entries, function(e) {
  heads <- paste0("    [", js_string(e$heads), ",", js_string(e$ids), "]", collapse = ",\n")
  sprintf("  {p:%s,n:%d,t:%s,h:[\n%s\n  ]}", js_string(e$page), e$number, js_string(e$title), heads)
}, character(1))
block <- c("<script>", "var trexSearchIndex = [", paste0(rows, collapse = ",\n"), "];", "</script>")

nav <- readLines(nav_file, warn = FALSE)
begin <- grep("^<!-- search index begin -->$", nav)
end <- grep("^<!-- search index end -->$", nav)
stopifnot(length(begin) == 1, length(end) == 1, begin < end)
writeLines(c(nav[seq_len(begin)], block, nav[end:length(nav)]), nav_file)
cat(sprintf("search index: %d pages, %d headings, written into %s\n", length(entries),
            sum(vapply(entries, function(e) length(e$heads), integer(1))), basename(nav_file)))
