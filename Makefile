.PHONY: build vignettes docs

version := $(shell less DESCRIPTION | grep 'Version' | sed 's/Version: \(.*\)$$/\1/')
pkg := build/slendr_$(version).tar.gz
logo := man/figures/logo.png

docs:
	rm -rf docs/reference
	R -e 'devtools::document()'
	R -e 'pkgdown::build_reference()'
	R -e 'pkgdown::build_reference_index()'

website:
	rm -rf docs/
	R -e 'devtools::install(upgrade = "never")'
	R -e 'source("logo.R")'
	R -e 'knitr::knit("README.Rmd", output = "README.md")'
	R -e 'devtools::document()'
	R -e 'pkgdown::build_site()'
	git restore docs/pkgdown.yml
	git restore docs/CNAME

build: $(pkg)

check: $(pkg)
	unset R_HAS_GGTREE; cd build; R CMD check --as-cran $(notdir $<)

winrel: README.md
	unset R_HAS_GGTREE; R -e 'devtools::check_win_release()'

windev: README.md
	unset R_HAS_GGTREE; R -e 'devtools::check_win_devel()'

winold: README.md
	unset R_HAS_GGTREE; R -e 'devtools::check_win_oldrelease()'

rhub: README.md
	unset R_HAS_GGTREE; R -e 'rhub::check_for_cran(platforms = c("macos-highsierra-release", "macos-highsierra-release-cran", "macos-m1-bigsur-release"))'

$(pkg): README.md
	R -e 'devtools::document()'
	unset R_HAS_GGTREE; mkdir -p build; cd build; R CMD build --log ../../slendr

README.md: README.Rmd $(logo)
	R -e 'devtools::install(upgrade = "never")'
	R -e 'knitr::knit("README.Rmd", output = "README.md")'

$(logo): logo.R
	Rscript logo.R

example_data:
	rm -rf inst/extdata/models/; mkdir -p inst/extdata/models/
	Rscript generate_examples.R

clean:
	rm -rf build
