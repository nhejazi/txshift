md:
	Rscript -e "rmarkdown::render('README.Rmd', output_file = 'README.md')"

site:
	Rscript -e "pkgdown::build_site()"

fastcheck:
	Rscript -e "devtools::check(build_args = '--no-build-vignettes')"

check:
	Rscript -e "devtools::check()"

test:
	Rscript -e "devtools::test()"

doc:
	Rscript -e "devtools::document()"

build:
	Rscript -e "devtools::build()"

style:
	Rscript -e "styler::style_pkg()"

